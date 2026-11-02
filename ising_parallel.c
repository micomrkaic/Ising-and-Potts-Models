#define _GNU_SOURCE
/*
 * ising_parallel.c — 2-D Ising and q-state Potts models, multithreaded
 *                    checkerboard Monte Carlo with SDL2 visualisation
 *
 * Decomposition
 *   A sweep is split into two phases by sublattice colour c = (i+j) mod 2.
 *   In a phase only sites of one colour are updated; all their nearest
 *   neighbours have the other colour and are frozen, so the ΔE reads are
 *   race-free without locks.  Within a phase the lattice is cut into
 *   horizontal row bands, one per thread.  Any thread count works.
 *
 * Algorithms
 *   Metropolis single-site update (parallel), Wolff single-cluster update
 *   (serial; the pixel repaint and the measurement are still parallel).
 *
 * Build
 *   make            (see Makefile)
 *
 * Headless physics check
 *   ./ising_parallel --batch -L 64 -T 2.0 -n 20000
 *   ./ising_parallel --batch -m potts -q 3 --scan 0.9:1.1:9 --csv
 */

#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <getopt.h>
#include <math.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define PROG        "ising_parallel"
#define VERSION     "2.1.1"

/* ── compile-time limits ─────────────────────────────────────────────────── */

#define MAXQ        16
#define MAXTHREADS  64
#define STATS_BUF   2000          /* GUI history window (frames)             */
#define NBLOCKS     20            /* jackknife blocks                        */

#define PANEL_W     340
#define SIM_W       800
#define WIN_W       (SIM_W + PANEL_W)
#define WIN_H       800
#define PLOT_H_PX   68
#define PLOT_W_PX   (PANEL_W - 24)

/* ── physics ─────────────────────────────────────────────────────────────── */

#define J           1.0
#define TC_ISING    2.269185314213022   /* 2/ln(1+√2)                        */

typedef enum { MODEL_ISING = 0, MODEL_POTTS = 1 } Model;

static double model_tc(Model m, int q) {
    return (m == MODEL_ISING) ? TC_ISING : 1.0 / log(1.0 + sqrt((double)q));
}
static const char *model_name(Model m) {
    return (m == MODEL_ISING) ? "ising" : "potts";
}

/* Complete elliptic integral of the first kind K(k), via AGM.             */
static double ellipK(double k) {
    double a = 1.0, b = sqrt(1.0 - k * k);
    if (b < 1e-300) return INFINITY;
    for (int it = 0; it < 60 && fabs(a - b) > 1e-16 * a; it++) {
        double an = 0.5 * (a + b);
        b = sqrt(a * b);
        a = an;
    }
    return M_PI / (2.0 * a);
}

/* Onsager: exact energy per spin of the infinite 2-D Ising lattice.       */
static double onsager_energy(double T) {
    double K = J / T;
    double s = sinh(2 * K), c = cosh(2 * K), t = tanh(2 * K);
    double kappa  = 2 * s / (c * c);
    double kappap = 2 * t * t - 1;
    if (fabs(kappap) < 1e-9) return -J * (c / s);     /* T = Tc: K diverges */
    return -J * (c / s) * (1.0 + (2.0 / M_PI) * kappap * ellipK(kappa));
}

/* Onsager–Yang: exact spontaneous magnetisation, zero for T ≥ Tc.         */
static double onsager_mag(double T) {
    if (T >= TC_ISING) return 0.0;
    double x = sinh(2 * J / T);
    return pow(1.0 - 1.0 / (x * x * x * x), 0.125);
}

/* ── lattice ─────────────────────────────────────────────────────────────── */

static int       L;                 /* side                                 */
static long      N;                 /* L*L                                  */
static int8_t   *spin;              /* Ising: ±1   Potts: 0..q-1            */
static uint32_t *pixels;            /* BGRA, L*L                            */
static Model     model = MODEL_ISING;
static int       q     = 3;

#define IDX(i, j) ((long)(i) * L + (j))

/* ── RNG: xoshiro256+ ────────────────────────────────────────────────────── */

typedef struct { uint64_t s[4]; } Rng;

static inline uint64_t rotl64(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

static inline uint64_t rng_next(Rng *r) {
    uint64_t res = rotl64(r->s[0] + r->s[3], 23) + r->s[0];
    uint64_t t   = r->s[1] << 17;
    r->s[2] ^= r->s[0]; r->s[3] ^= r->s[1];
    r->s[1] ^= r->s[2]; r->s[0] ^= r->s[3];
    r->s[2] ^= t;
    r->s[3]  = rotl64(r->s[3], 45);
    return res;
}
static inline double rng_double(Rng *r) {
    return (double)(rng_next(r) >> 11) * (1.0 / 9007199254740992.0);
}
static void rng_seed(Rng *r, uint64_t seed) {
    for (int i = 0; i < 4; i++) {                  /* splitmix64 */
        seed += 0x9e3779b97f4a7c15ULL;
        uint64_t z = seed;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        r->s[i] = z ^ (z >> 31);
    }
}

static Rng main_rng;                /* Wolff, initial configurations        */

/* ── Boltzmann tables ────────────────────────────────────────────────────── */
/* Ising: ΔE = 2J·s·Σnn, s·Σnn ∈ {-4..4} even → index s·Σnn/2+2 ∈ 0..4
   Potts: ΔE = J·(n_old − n_new) ∈ {-4..4}      → index ΔE+4     ∈ 0..8  */

static double boltz_i[5];
static double boltz_p[9];

static void build_tables(double T) {
    for (int k = 0; k < 5; k++) {
        double dE = (k - 2) * 4.0 * J;
        boltz_i[k] = (dE <= 0) ? 1.0 : exp(-dE / T);
    }
    for (int k = 0; k < 9; k++) {
        double dE = (k - 4) * J;
        boltz_p[k] = (dE <= 0) ? 1.0 : exp(-dE / T);
    }
}

/* ── barrier: spin briefly, then sleep ───────────────────────────────────── */

#define SPIN_LIMIT 4000

static inline void cpu_relax(void) {
#if defined(__x86_64__) || defined(__i386__)
    __builtin_ia32_pause();
#elif defined(__aarch64__)
    __asm__ __volatile__("yield");
#endif
}

typedef struct {
    pthread_mutex_t mtx;
    pthread_cond_t  cv;
    atomic_int      count;
    atomic_int      gen;
    int             parties;
} Barrier;

static void barrier_init(Barrier *b, int parties) {
    pthread_mutex_init(&b->mtx, NULL);
    pthread_cond_init(&b->cv, NULL);
    atomic_store(&b->count, 0);
    atomic_store(&b->gen, 0);
    b->parties = parties;
}

static void barrier_wait(Barrier *b) {
    int gen = atomic_load(&b->gen);
    if (atomic_fetch_add(&b->count, 1) + 1 == b->parties) {
        atomic_store(&b->count, 0);
        pthread_mutex_lock(&b->mtx);
        atomic_store(&b->gen, gen + 1);
        pthread_cond_broadcast(&b->cv);
        pthread_mutex_unlock(&b->mtx);
        return;
    }
    for (int i = 0; i < SPIN_LIMIT; i++) {
        if (atomic_load(&b->gen) != gen) return;
        cpu_relax();
    }
    pthread_mutex_lock(&b->mtx);
    while (atomic_load(&b->gen) == gen)
        pthread_cond_wait(&b->cv, &b->mtx);
    pthread_mutex_unlock(&b->mtx);
}

/* ── thread pool ─────────────────────────────────────────────────────────── */

typedef enum { CMD_SWEEP, CMD_PAINT, CMD_MEASURE, CMD_QUIT } Cmd;
typedef enum { COLOR_SPIN = 0, COLOR_THREAD = 1, COLOR_ENERGY = 2, COLOR_COUNT } ColorMode;

typedef struct {
    int   tid, r0, r1;              /* row band [r0, r1)                    */
    Rng   rng;
    long  sum_s;                    /* Ising: Σσ                            */
    long  cnt[MAXQ];                /* Potts: state counts                  */
    long  bonds;                    /* Σ (right+down) bond terms            */
} __attribute__((aligned(64))) Ctx;

static int        nthreads;
static Ctx        ctx[MAXTHREADS];
static pthread_t  tids[MAXTHREADS];
static Barrier    bar;
static Cmd        g_cmd;
static int        g_phase;
static ColorMode  g_color = COLOR_SPIN;

static const uint32_t thread_palette[8] = {    /* BGRA */
    0xFFD48B3Bu, 0xFF3C59E8u, 0xFF759E1Du, 0xFF17A7BAu,
    0xFFDD777Fu, 0xFF7E53D4u, 0xFF4A4AE2u, 0xFF888888u,
};
static uint32_t potts_palette[MAXQ];

static uint32_t hsv_bgra(double h, double s, double v) {
    int    hi = (int)floor(h / 60.0) % 6;
    double f  = h / 60.0 - floor(h / 60.0);
    double p = v * (1 - s), t = v * (1 - (1 - f) * s), u = v * (1 - f * s);
    double r, g, b;
    switch (hi) {
        case 0: r = v; g = t; b = p; break;
        case 1: r = u; g = v; b = p; break;
        case 2: r = p; g = v; b = t; break;
        case 3: r = p; g = u; b = v; break;
        case 4: r = t; g = p; b = v; break;
        default: r = v; g = p; b = u; break;
    }
    return 0xFF000000u | ((uint32_t)(r * 255) << 16) | ((uint32_t)(g * 255) << 8) | (uint32_t)(b * 255);
}
static void build_potts_palette(void) {
    for (int k = 0; k < MAXQ; k++)
        potts_palette[k] = hsv_bgra(360.0 * k / q, 0.78, (k & 1) ? 0.95 : 0.70);
}

/* ── Metropolis: one colour phase over this thread's band ────────────────── */

static void sweep_ising(Ctx *c, int phase) {
    Rng *rng = &c->rng;
    for (int i = c->r0; i < c->r1; i++) {
        const int8_t *up  = spin + IDX(i == 0 ? L - 1 : i - 1, 0);
        const int8_t *dn  = spin + IDX(i == L - 1 ? 0 : i + 1, 0);
        int8_t       *row = spin + IDX(i, 0);
        for (int j = (i + phase) & 1; j < L; j += 2) {
            int jl = j ? j - 1 : L - 1;
            int jr = (j == L - 1) ? 0 : j + 1;
            int s  = row[j];
            int nn = up[j] + dn[j] + row[jl] + row[jr];
            int sn = s * nn;                              /* even, -4..4 */
            if (sn <= 0 || rng_double(rng) < boltz_i[sn / 2 + 2])
                row[j] = (int8_t)-s;
        }
    }
}

static void sweep_potts(Ctx *c, int phase) {
    Rng *rng = &c->rng;
    const int qm1 = q - 1;
    for (int i = c->r0; i < c->r1; i++) {
        const int8_t *up  = spin + IDX(i == 0 ? L - 1 : i - 1, 0);
        const int8_t *dn  = spin + IDX(i == L - 1 ? 0 : i + 1, 0);
        int8_t       *row = spin + IDX(i, 0);
        for (int j = (i + phase) & 1; j < L; j += 2) {
            int jl = j ? j - 1 : L - 1;
            int jr = (j == L - 1) ? 0 : j + 1;
            int s  = row[j];
            int r  = s + 1 + (int)(rng_next(rng) % (uint64_t)qm1);
            if (r >= q) r -= q;                           /* uniform ≠ s  */
            int a = up[j], b = dn[j], cc = row[jl], d = row[jr];
            int nold = (a == s) + (b == s) + (cc == s) + (d == s);
            int nnew = (a == r) + (b == r) + (cc == r) + (d == r);
            int dE   = nold - nnew;                       /* units of J   */
            if (dE <= 0 || rng_double(rng) < boltz_p[dE + 4])
                row[j] = (int8_t)r;
        }
    }
}

/* ── measurement: partial sums over this thread's band ───────────────────── */

static void measure_band(Ctx *c) {
    long sum_s = 0, bonds = 0;
    long cnt[MAXQ] = {0};
    for (int i = c->r0; i < c->r1; i++) {
        const int8_t *row = spin + IDX(i, 0);
        const int8_t *dn  = spin + IDX(i == L - 1 ? 0 : i + 1, 0);
        if (model == MODEL_ISING) {
            for (int j = 0; j < L; j++) {
                int s = row[j];
                sum_s += s;
                bonds += s * (row[j == L - 1 ? 0 : j + 1] + dn[j]);
            }
        } else {
            for (int j = 0; j < L; j++) {
                int s = row[j];
                cnt[s]++;
                bonds += (row[j == L - 1 ? 0 : j + 1] == s) + (dn[j] == s);
            }
        }
    }
    c->sum_s = sum_s;
    c->bonds = bonds;
    memcpy(c->cnt, cnt, sizeof cnt);
}

/* ── paint: this thread's rows of the pixel buffer ───────────────────────── */

static inline uint32_t energy_color(int unsat, int max) {   /* 0..max */
    double t = (double)unsat / max;
    uint8_t r = (uint8_t)(t * 210), b = (uint8_t)((1 - t) * 210);
    return 0xFF000000u | ((uint32_t)r << 16) | (uint32_t)b;
}

static void paint_band(Ctx *c) {
    for (int i = c->r0; i < c->r1; i++) {
        const int8_t *row = spin + IDX(i, 0);
        uint32_t     *out = pixels + IDX(i, 0);
        switch (g_color) {
        case COLOR_SPIN:
            if (model == MODEL_ISING)
                for (int j = 0; j < L; j++) out[j] = row[j] > 0 ? 0xFFF5F0FFu : 0xFF0A0510u;
            else
                for (int j = 0; j < L; j++) out[j] = potts_palette[row[j]];
            break;
        case COLOR_THREAD: {
            /* spins drawn in the owning thread's hue: bright = up / high
               state, dark = down / state 0, so the row bands of the
               decomposition and the domain structure are both visible */
            uint32_t col  = thread_palette[c->tid & 7];
            uint32_t dark = (col >> 2) & 0x3F3F3Fu;      /* ~25 % brightness */
            if (model == MODEL_ISING) {
                for (int j = 0; j < L; j++)
                    out[j] = 0xFF000000u | (row[j] > 0 ? (col & 0xFFFFFFu) : dark);
            } else {
                for (int j = 0; j < L; j++) {
                    double   t = 0.25 + 0.75 * row[j] / (double)(q - 1);
                    uint32_t r = (uint32_t)(((col >> 16) & 0xFF) * t);
                    uint32_t g = (uint32_t)(((col >>  8) & 0xFF) * t);
                    uint32_t b = (uint32_t)(( col        & 0xFF) * t);
                    out[j] = 0xFF000000u | (r << 16) | (g << 8) | b;
                }
            }
            break;
        }
        case COLOR_ENERGY: {
            const int8_t *up = spin + IDX(i == 0 ? L - 1 : i - 1, 0);
            const int8_t *dn = spin + IDX(i == L - 1 ? 0 : i + 1, 0);
            for (int j = 0; j < L; j++) {
                int jl = j ? j - 1 : L - 1, jr = (j == L - 1) ? 0 : j + 1;
                int s = row[j], a = up[j], b = dn[j], cc = row[jl], d = row[jr];
                int unsat = (model == MODEL_ISING)
                          ? (4 - s * (a + b + cc + d)) / 2
                          : 4 - ((a == s) + (b == s) + (cc == s) + (d == s));
                out[j] = energy_color(unsat, 4);
            }
            break;
        }
        default: break;
        }
    }
}

/* ── worker loop and dispatch ────────────────────────────────────────────── */

static void run_cmd(Ctx *c, Cmd cmd, int phase) {
    switch (cmd) {
    case CMD_SWEEP:
        if (model == MODEL_ISING) sweep_ising(c, phase);
        else                      sweep_potts(c, phase);
        break;
    case CMD_PAINT:   paint_band(c);   break;
    case CMD_MEASURE: measure_band(c); break;
    default: break;
    }
}

/* Band 0 belongs to the main thread; bands 1..nthreads-1 to workers.  The
   main thread never idles while others work, so the machine is never
   oversubscribed by a spinning coordinator.                                */
static void *worker(void *arg) {
    Ctx *c = arg;
    for (;;) {
        barrier_wait(&bar);                 /* wait for a command            */
        Cmd cmd = g_cmd;
        if (cmd == CMD_QUIT) return NULL;
        run_cmd(c, cmd, g_phase);
        barrier_wait(&bar);                 /* report completion             */
    }
}

static void dispatch(Cmd cmd, int phase) {
    g_cmd   = cmd;
    g_phase = phase;
    if (nthreads > 1) barrier_wait(&bar);   /* release the workers           */
    if (cmd == CMD_QUIT) return;
    run_cmd(&ctx[0], cmd, phase);
    if (nthreads > 1) barrier_wait(&bar);   /* wait for the workers          */
}

static void pool_start(int nt, uint64_t seed) {
    nthreads = nt;
    barrier_init(&bar, nt);
    for (int t = 0; t < nt; t++) {
        ctx[t].tid = t;
        ctx[t].r0  = (int)((long)L * t / nt);
        ctx[t].r1  = (int)((long)L * (t + 1) / nt);
        rng_seed(&ctx[t].rng, seed ^ (0x9e3779b97f4a7c15ULL * (uint64_t)(t + 1)));
        if (t > 0) pthread_create(&tids[t], NULL, worker, &ctx[t]);
    }
}

static void pool_stop(void) {
    dispatch(CMD_QUIT, 0);
    for (int t = 1; t < nthreads; t++) pthread_join(tids[t], NULL);
}

/* ── full sweeps ─────────────────────────────────────────────────────────── */

static void metropolis_sweep(void) {
    dispatch(CMD_SWEEP, 0);
    dispatch(CMD_SWEEP, 1);
}

/* ── Wolff single-cluster update (serial) ────────────────────────────────── */

static long     *wstack;
static uint32_t *wmark;
static uint32_t  wgen = 0;
static double    p_add;             /* bond activation probability          */
static long      wolff_steps = 1;   /* cluster flips per "sweep"            */
static double    wolff_avg   = 0;   /* running mean cluster size            */

static void wolff_set_T(double T) {
    p_add = (model == MODEL_ISING) ? 1.0 - exp(-2.0 * J / T)
                                   : 1.0 - exp(-J / T);
    wolff_avg = 0;                  /* re-calibrate the sweep length        */
}

/* Returns the cluster size. */
static long wolff_step(void) {
    if (++wgen == 0) { memset(wmark, 0, (size_t)N * sizeof *wmark); wgen = 1; }

    long seed = (long)(rng_double(&main_rng) * N);
    int  s0   = spin[seed];
    int  s1;
    if (model == MODEL_ISING) s1 = -s0;
    else { s1 = s0 + 1 + (int)(rng_next(&main_rng) % (uint64_t)(q - 1)); if (s1 >= q) s1 -= q; }

    long top = 0, size = 1;
    wstack[top++] = seed;
    wmark[seed]   = wgen;
    spin[seed]    = (int8_t)s1;

    while (top > 0) {
        long cur = wstack[--top];
        int  i = (int)(cur / L), j = (int)(cur % L);
        long nb[4] = {
            IDX(i == 0 ? L - 1 : i - 1, j), IDX(i == L - 1 ? 0 : i + 1, j),
            IDX(i, j ? j - 1 : L - 1),      IDX(i, j == L - 1 ? 0 : j + 1),
        };
        for (int k = 0; k < 4; k++) {
            long n = nb[k];
            if (wmark[n] != wgen && spin[n] == s0 && rng_double(&main_rng) < p_add) {
                wmark[n] = wgen;
                spin[n]  = (int8_t)s1;
                wstack[top++] = n;
                size++;
            }
        }
    }
    return size;
}

/* Set the number of cluster flips per sweep to about N/⟨cluster⟩, using
   only data from earlier sweeps.  The count must be fixed *before* a sweep
   whose configuration is measured: stopping "once N spins have flipped"
   is a path-dependent stopping time and biases the sample towards the
   ordered states that follow a large cluster flip (observed as U_L 0.632
   instead of 0.611 at T_c).                                                 */
static void wolff_calibrate(void) {
    if (wolff_avg <= 0) {                       /* first call at this T     */
        long f = 0;
        for (int k = 0; k < 10; k++) f += wolff_step();
        wolff_avg = f / 10.0;
    }
    wolff_steps = lround((double)N / wolff_avg);
    if (wolff_steps < 1) wolff_steps = 1;
}

/* One Wolff "sweep": exactly wolff_steps cluster flips.  Returns the number
   of spins flipped; *avg receives this sweep's mean cluster size.           */
static long wolff_sweep(double *avg) {
    long flipped = 0;
    for (long k = 0; k < wolff_steps; k++) flipped += wolff_step();
    double a  = (double)flipped / wolff_steps;
    wolff_avg = 0.9 * wolff_avg + 0.1 * a;      /* EMA for later calibration */
    if (avg) *avg = a;
    return flipped;
}

/* ── observables ─────────────────────────────────────────────────────────── */

typedef struct { double m, e; } Obs;     /* order parameter, energy/spin   */

static Obs measure(void) {
    dispatch(CMD_MEASURE, 0);
    long bonds = 0;
    Obs o;
    if (model == MODEL_ISING) {
        long sum_s = 0;
        for (int t = 0; t < nthreads; t++) { sum_s += ctx[t].sum_s; bonds += ctx[t].bonds; }
        o.m = fabs((double)sum_s / N);
    } else {
        long cnt[MAXQ] = {0}, mx = 0;
        for (int t = 0; t < nthreads; t++) {
            bonds += ctx[t].bonds;
            for (int k = 0; k < q; k++) cnt[k] += ctx[t].cnt[k];
        }
        for (int k = 0; k < q; k++) if (cnt[k] > mx) mx = cnt[k];
        o.m = (q * (double)mx / N - 1.0) / (q - 1);
    }
    o.e = -J * (double)bonds / N;
    return o;
}

/* Jackknife over NBLOCKS blocks of a time series of (|m|, e).
   binder = 1 − ⟨m⁴⟩ / (3⟨m²⟩²), the Binder cumulant U_L.  Curves of U_L
   against T for different L cross at T_c; for the 2-D Ising model with
   periodic boundaries the crossing value is U* ≈ 0.6107.                  */
typedef struct {
    double m, m_err, e, e_err, chi, chi_err, cv, cv_err, binder, binder_err, tau_m;
    long   n;
} Stats;

static Stats analyse(const double *m, const double *e, long n, double T) {
    Stats st = {0};
    st.n = n;
    if (n < 2) { st.m = n ? m[0] : NAN; st.e = n ? e[0] : NAN; st.binder = NAN; return st; }

    int  B  = (n >= 2 * NBLOCKS) ? NBLOCKS : (int)(n / 2);
    long bs = n / B;
    n       = bs * B;                                /* drop remainder */

    double Sm[NBLOCKS], Sm2[NBLOCKS], Sm4[NBLOCKS], Se[NBLOCKS], Se2[NBLOCKS];
    double Tm = 0, Tm2 = 0, Tm4 = 0, Te = 0, Te2 = 0;
    for (int b = 0; b < B; b++) {
        Sm[b] = Sm2[b] = Sm4[b] = Se[b] = Se2[b] = 0;
        for (long k = b * bs; k < (b + 1) * bs; k++) {
            double m2 = m[k] * m[k];
            Sm[b] += m[k]; Sm2[b] += m2; Sm4[b] += m2 * m2;
            Se[b] += e[k]; Se2[b] += e[k] * e[k];
        }
        Tm += Sm[b]; Tm2 += Sm2[b]; Tm4 += Sm4[b]; Te += Se[b]; Te2 += Se2[b];
    }
    double inv = 1.0 / n;
    st.m   = Tm * inv;
    st.e   = Te * inv;
    st.chi = N * (Tm2 * inv - st.m * st.m) / T;
    st.cv  = N * (Te2 * inv - st.e * st.e) / (T * T);
    st.binder = 1.0 - (Tm4 * inv) / (3.0 * (Tm2 * inv) * (Tm2 * inv));

    double jm = 0, je = 0, jchi = 0, jcv = 0, ju = 0;
    double invj = 1.0 / (n - bs);
    for (int b = 0; b < B; b++) {
        double mm = (Tm - Sm[b]) * invj, mm2 = (Tm2 - Sm2[b]) * invj, mm4 = (Tm4 - Sm4[b]) * invj;
        double ee = (Te - Se[b]) * invj, ee2 = (Te2 - Se2[b]) * invj;
        double chi_b = N * (mm2 - mm * mm) / T;
        double cv_b  = N * (ee2 - ee * ee) / (T * T);
        double u_b   = 1.0 - mm4 / (3.0 * mm2 * mm2);
        jm += (mm - st.m) * (mm - st.m);
        je += (ee - st.e) * (ee - st.e);
        jchi += (chi_b - st.chi) * (chi_b - st.chi);
        jcv  += (cv_b - st.cv) * (cv_b - st.cv);
        ju   += (u_b - st.binder) * (u_b - st.binder);
    }
    double f = (double)(B - 1) / B;
    st.m_err = sqrt(f * jm);  st.e_err = sqrt(f * je);
    st.chi_err = sqrt(f * jchi); st.cv_err = sqrt(f * jcv);
    st.binder_err = sqrt(f * ju);

    /* integrated autocorrelation time of |m| from the blocking ratio       */
    double var = Tm2 * inv - st.m * st.m;
    st.tau_m = (var > 0) ? 0.5 * ((double)n * st.m_err * st.m_err / var - 1.0) : 0.0;
    if (st.tau_m < 0) st.tau_m = 0;
    return st;
}

/* ── lattice initialisation ──────────────────────────────────────────────── */

static void lattice_hot(void) {
    if (model == MODEL_ISING)
        for (long k = 0; k < N; k++) spin[k] = (rng_next(&main_rng) & 1) ? 1 : -1;
    else
        for (long k = 0; k < N; k++) spin[k] = (int8_t)(rng_next(&main_rng) % (uint64_t)q);
}
static void lattice_cold(void) {
    memset(spin, model == MODEL_ISING ? 1 : 0, (size_t)N);
}

static void lattice_alloc(int side) {
    L = side;
    N = (long)L * L;
    spin   = malloc((size_t)N);
    pixels = malloc((size_t)N * sizeof *pixels);
    wstack = malloc((size_t)N * sizeof *wstack);
    wmark  = calloc((size_t)N, sizeof *wmark);
    if (!spin || !pixels || !wstack || !wmark) { fprintf(stderr, "out of memory\n"); exit(1); }
}

/* ── batch mode ──────────────────────────────────────────────────────────── */

typedef struct {
    int    L, threads, q, use_wolff, batch, csv, frames;
    Model  model;
    double T, scan_lo, scan_hi;
    int    scan_n;
    long   sweeps, therm;
    uint64_t seed;
    bool   T_set;
} Opts;

static double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static void batch_header(const Opts *o) {
    if (o->csv)
        puts("model,q,L,threads,alg,T,Tc,sweeps,m,m_err,e,e_err,chi,chi_err,cv,cv_err,"
             "U,U_err,tau_m,m_exact,e_exact,avg_cluster,secs,sweeps_per_s");
    else
        printf("%-6s %2s %5s %3s %-5s %8s %8s %7s  %-18s %-18s %-20s %-20s %-18s %7s %8s %8s %8s %8s\n",
               "model", "q", "L", "thr", "alg", "T", "T/Tc", "sweeps", "m", "e",
               "chi", "Cv", "U_L", "tau_m", "m_exact", "e_exact", "secs", "sw/s");
}

static void batch_point(const Opts *o, double T) {
    build_tables(T);
    wolff_set_T(T);
    double avgc = NAN, acc = 0; long nc = 0;

    /* thermalise; the Wolff sweep length adapts here and is then frozen   */
    for (long k = 0; k < o->therm; k++) {
        if (o->use_wolff) { wolff_calibrate(); wolff_sweep(NULL); }
        else metropolis_sweep();
    }
    if (o->use_wolff) wolff_calibrate();
    double *ms = malloc((size_t)o->sweeps * sizeof *ms);
    double *es = malloc((size_t)o->sweeps * sizeof *es);
    if (!ms || !es) { fprintf(stderr, "out of memory\n"); exit(1); }

    double t0 = now_s();
    for (long k = 0; k < o->sweeps; k++) {
        if (o->use_wolff) { double a; wolff_sweep(&a); acc += a; nc++; }
        else metropolis_sweep();
        Obs ob = measure();
        ms[k] = ob.m; es[k] = ob.e;
    }
    double secs = now_s() - t0;
    if (nc) avgc = acc / nc;

    Stats st = analyse(ms, es, o->sweeps, T);
    double tc = model_tc(model, q);
    double mex = (model == MODEL_ISING) ? onsager_mag(T) : NAN;
    double eex = (model == MODEL_ISING) ? onsager_energy(T) : NAN;
    const char *alg = o->use_wolff ? "wolff" : "metro";

    if (o->csv)
        printf("%s,%d,%d,%d,%s,%.6f,%.6f,%ld,%.6f,%.6f,%.6f,%.6f,%.4f,%.4f,%.4f,%.4f,"
               "%.6f,%.6f,%.3f,%.6f,%.6f,%.2f,%.3f,%.1f\n",
               model_name(model), model == MODEL_ISING ? 2 : q, L, nthreads, alg, T, tc,
               o->sweeps, st.m, st.m_err, st.e, st.e_err, st.chi, st.chi_err,
               st.cv, st.cv_err, st.binder, st.binder_err, st.tau_m, mex, eex, avgc,
               secs, o->sweeps / secs);
    else
        printf("%-6s %2d %5d %3d %-5s %8.5f %8.5f %7ld  %8.5f ±%8.5f %8.5f ±%8.5f "
               "%9.3f ±%9.3f %9.4f ±%9.4f %8.5f ±%8.5f %7.2f %8.5f %8.5f %8.2f %8.0f\n",
               model_name(model), model == MODEL_ISING ? 2 : q, L, nthreads, alg, T, T / tc,
               o->sweeps, st.m, st.m_err, st.e, st.e_err, st.chi, st.chi_err,
               st.cv, st.cv_err, st.binder, st.binder_err, st.tau_m, mex, eex, secs,
               o->sweeps / secs);
    fflush(stdout);
    free(ms); free(es);
}

static int run_batch(const Opts *o) {
    batch_header(o);
    if (o->scan_n > 0) {
        for (int k = 0; k < o->scan_n; k++) {
            double T = (o->scan_n == 1) ? o->scan_lo
                     : o->scan_lo + (o->scan_hi - o->scan_lo) * k / (o->scan_n - 1);
            batch_point(o, T);
        }
    } else {
        batch_point(o, o->T);
    }
    return 0;
}

/* ── SDL helpers ─────────────────────────────────────────────────────────── */

static void set_color(SDL_Renderer *r, uint8_t R, uint8_t G, uint8_t B, uint8_t A) {
    SDL_SetRenderDrawColor(r, R, G, B, A);
}
static void fill_rect(SDL_Renderer *r, int x, int y, int w, int h) {
    SDL_Rect rc = {x, y, w, h};
    SDL_RenderFillRect(r, &rc);
}
static void hline(SDL_Renderer *r, int x1, int x2, int y) { SDL_RenderDrawLine(r, x1, y, x2, y); }

static void draw_text(SDL_Renderer *rend, TTF_Font *font, const char *text,
                      int x, int y, uint8_t R, uint8_t G, uint8_t B) {
    if (!font) return;
    SDL_Color col = {R, G, B, 255};
    SDL_Surface *surf = TTF_RenderUTF8_Blended(font, text, col);
    if (!surf) return;
    SDL_Texture *tex = SDL_CreateTextureFromSurface(rend, surf);
    if (tex) {
        SDL_Rect dst = {x, y, surf->w, surf->h};
        SDL_RenderCopy(rend, tex, NULL, &dst);
        SDL_DestroyTexture(tex);
    }
    SDL_FreeSurface(surf);
}

static TTF_Font *open_font(int pt) {
    const char *env = getenv("ISING_FONT");
    if (env) return TTF_OpenFont(env, pt);
    static const char *paths[] = {
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
        "/usr/share/fonts/truetype/ubuntu/Ubuntu-R.ttf",
        "/usr/share/fonts/TTF/DejaVuSans.ttf",
        "/usr/share/fonts/dejavu/DejaVuSans.ttf",
        "/opt/homebrew/share/fonts/dejavu/DejaVuSans.ttf",
        "/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
        "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/Library/Fonts/Arial.ttf",
        NULL
    };
    for (int i = 0; paths[i]; i++) {
        TTF_Font *f = TTF_OpenFont(paths[i], pt);
        if (f) return f;
    }
    return NULL;
}

/* ── GUI state ───────────────────────────────────────────────────────────── */

typedef struct {
    double T;
    Obs    cur;
    float  hist_m[STATS_BUF], hist_e[STATS_BUF], hist_chi[STATS_BUF], hist_cv[STATS_BUF];
    int    head, count;
    Stats  st;
    long   total_sweeps;
    int    use_wolff, sweeps_per_frame, paused;
    double avg_cluster;
} Sim;

static void sim_reset_hist(Sim *s) { s->head = s->count = 0; memset(&s->st, 0, sizeof s->st); }

static void sim_set_T(Sim *s, double T) {
    double tc = model_tc(model, q);
    if (T < 0.05 * tc) T = 0.05 * tc;
    if (T > 5.0 * tc)  T = 5.0 * tc;
    s->T = T;
    build_tables(T);
    wolff_set_T(T);
    sim_reset_hist(s);
}

static void sim_record(Sim *s) {
    s->hist_m[s->head] = (float)s->cur.m;
    s->hist_e[s->head] = (float)s->cur.e;
    /* stats over the window */
    static double mb[STATS_BUF], eb[STATS_BUF];
    int n = s->count < STATS_BUF ? s->count + 1 : STATS_BUF;
    int start = (s->head + 1 - n + STATS_BUF) % STATS_BUF;
    for (int k = 0; k < n; k++) {
        int ri = (start + k) % STATS_BUF;
        mb[k] = s->hist_m[ri]; eb[k] = s->hist_e[ri];
    }
    s->st = analyse(mb, eb, n, s->T);
    s->hist_chi[s->head] = (float)s->st.chi;
    s->hist_cv[s->head]  = (float)s->st.cv;
    s->head = (s->head + 1) % STATS_BUF;
    if (s->count < STATS_BUF) s->count++;
}

static void draw_plot(SDL_Renderer *rend, TTF_Font *font, const float *buf, int head, int count,
                      int x0, int y0, int pw, int ph, float vmin, float vmax,
                      uint8_t R, uint8_t G, uint8_t B, const char *label, float ref) {
    set_color(rend, 12, 12, 22, 255);
    fill_rect(rend, x0, y0, pw, ph);
    set_color(rend, 50, 50, 70, 255);
    SDL_Rect border = {x0, y0, pw, ph};
    SDL_RenderDrawRect(rend, &border);
    if (!isnan(ref) && ref > vmin && ref < vmax) {
        float fy = y0 + ph - (ref - vmin) / (vmax - vmin) * ph;
        set_color(rend, 200, 160, 60, 160);
        hline(rend, x0, x0 + pw, (int)fy);
    }
    if (count >= 2) {
        set_color(rend, R, G, B, 220);
        float px = 0, py = 0;
        int   n  = count < pw ? count : pw;
        for (int k = 0; k < n; k++) {
            int   ri = (head - n + k + STATS_BUF) % STATS_BUF;
            float v  = buf[ri];
            float fx = x0 + (float)k * pw / n;
            float t  = (vmax > vmin) ? (v - vmin) / (vmax - vmin) : 0.5f;
            t = t < 0 ? 0 : t > 1 ? 1 : t;
            float fy = y0 + ph - t * ph;
            if (k) SDL_RenderDrawLineF(rend, px, py, fx, fy);
            px = fx; py = fy;
        }
    }
    draw_text(rend, font, label, x0 + 3, y0 + 2, R, G, B);
}

static const char *color_names[COLOR_COUNT] = {"spin", "thread", "energy"};

static void draw_panel(SDL_Renderer *rend, TTF_Font *fl, TTF_Font *fs, const Sim *s,
                       int recording, int rec_frames) {
    int  px0 = SIM_W + 8, y = 6;
    char buf[192];
    double tc = model_tc(model, q);

    set_color(rend, 12, 12, 22, 255);
    fill_rect(rend, SIM_W, 0, PANEL_W, WIN_H);
    set_color(rend, 55, 55, 80, 255);
    SDL_RenderDrawLine(rend, SIM_W, 0, SIM_W, WIN_H);

    if (model == MODEL_ISING) snprintf(buf, sizeof buf, "2D Ising model");
    else                      snprintf(buf, sizeof buf, "2D Potts model, q = %d", q);
    draw_text(rend, fl, buf, px0, y, 200, 180, 255); y += 24;

    snprintf(buf, sizeof buf, "T = %.4f J/kB", s->T);
    draw_text(rend, fs, buf, px0, y, 255, 255, 200); y += 15;
    snprintf(buf, sizeof buf, "T\xe2\x82\x9c = %.4f  (T/T\xe2\x82\x9c = %.3f)", tc, s->T / tc);
    bool near = fabs(s->T - tc) < 0.02 * tc;
    draw_text(rend, fs, buf, px0, y, near ? 255 : 180, near ? 100 : 180, near ? 50 : 100); y += 15;
    const char *phase = s->T < tc ? "Ordered" : "Disordered";
    const char *kind  = (model == MODEL_POTTS && q > 4) ? "first-order transition"
                                                        : "continuous transition";
    snprintf(buf, sizeof buf, "%s phase  (%s)", phase, kind);
    draw_text(rend, fs, buf, px0, y, s->T < tc ? 100 : 180, s->T < tc ? 220 : 140, 110); y += 16;

    /* temperature bar: 0 .. 2.5 Tc */
    double Tmax = 2.5 * tc;
    int bar_w = PANEL_W - 20, bar_x = SIM_W + 10;
    set_color(rend, 40, 40, 60, 255); fill_rect(rend, bar_x, y, bar_w, 8);
    int fill = (int)(s->T / Tmax * bar_w); if (fill > bar_w) fill = bar_w;
    set_color(rend, 100, 180, 255, 255); fill_rect(rend, bar_x, y, fill, 8);
    int tcx = bar_x + (int)(tc / Tmax * bar_w);
    set_color(rend, 255, 150, 50, 255); SDL_RenderDrawLine(rend, tcx, y - 2, tcx, y + 10);
    y += 16;

    snprintf(buf, sizeof buf, "%s   %d threads   colour: %s%s",
             s->use_wolff ? "Wolff" : "Metropolis", nthreads, color_names[g_color],
             s->paused ? "   PAUSED" : "");
    draw_text(rend, fs, buf, px0, y, 160, 200, 160); y += 14;
    snprintf(buf, sizeof buf, "L = %d   sweeps/frame %d   total %ld", L, s->sweeps_per_frame, s->total_sweeps);
    draw_text(rend, fs, buf, px0, y, 120, 120, 140); y += 14;
    if (s->use_wolff) snprintf(buf, sizeof buf, "window %d frames   \xe2\x9f\xa8""cluster\xe2\x9f\xa9 %.1f", s->count, s->avg_cluster);
    else              snprintf(buf, sizeof buf, "window %d frames", s->count);
    draw_text(rend, fs, buf, px0, y, 120, 120, 140); y += 16;

    set_color(rend, 50, 50, 75, 255); hline(rend, SIM_W + 4, WIN_W - 4, y); y += 7;

    const Stats *st = &s->st;
    struct { const char *nm; double v, err; uint8_t r, g, b; } obs[] = {
        {"m",                   s->cur.m, NAN,        220, 220, 255},
        {"\xe2\x9f\xa8m\xe2\x9f\xa9",   st->m, st->m_err,   180, 180, 220},
        {"e",                   s->cur.e, NAN,        255, 200, 150},
        {"\xe2\x9f\xa8""e\xe2\x9f\xa9", st->e, st->e_err,   220, 170, 120},
        {"\xcf\x87",            st->chi, st->chi_err, 150, 255, 200},
        {"Cv",                  st->cv,  st->cv_err,  255, 150, 200},
    };
    for (int k = 0; k < 6; k++) {
        if (isnan(obs[k].err)) snprintf(buf, sizeof buf, "%-6s %10.5f", obs[k].nm, obs[k].v);
        else snprintf(buf, sizeof buf, "%-6s %10.5f \xc2\xb1 %.5f", obs[k].nm, obs[k].v, obs[k].err);
        draw_text(rend, fs, buf, px0, y, obs[k].r, obs[k].g, obs[k].b); y += 14;
    }
    snprintf(buf, sizeof buf, "U\xe2\x82\x97 %10.5f \xc2\xb1 %.5f   (U* \xe2\x89\x88 0.611 at T\xe2\x82\x9c)",
             st->binder, st->binder_err);
    draw_text(rend, fs, buf, px0, y, 200, 200, 120); y += 14;
    snprintf(buf, sizeof buf, "\xcf\x84\xe2\x82\x98 (frames)  %8.2f", st->tau_m);
    draw_text(rend, fs, buf, px0, y, 200, 200, 120); y += 14;
    if (model == MODEL_ISING) {
        snprintf(buf, sizeof buf, "exact  m\xe2\x88\x9e = %.5f   e\xe2\x88\x9e = %.5f",
                 onsager_mag(s->T), onsager_energy(s->T));
        draw_text(rend, fs, buf, px0, y, 200, 160, 60); y += 14;
    } else {
        snprintf(buf, sizeof buf, "exact  e(T\xe2\x82\x9c) = %.5f", -J * (1.0 + 1.0 / sqrt((double)q)));
        draw_text(rend, fs, buf, px0, y, 200, 160, 60); y += 14;
    }
    y += 3;
    set_color(rend, 50, 50, 75, 255); hline(rend, SIM_W + 4, WIN_W - 4, y); y += 6;

    int pw = PLOT_W_PX, ph = PLOT_H_PX, px = SIM_W + 8;
    float mref = (model == MODEL_ISING) ? (float)onsager_mag(s->T) : NAN;
    float eref = (model == MODEL_ISING) ? (float)onsager_energy(s->T) : NAN;
    draw_plot(rend, fs, s->hist_m, s->head, s->count, px, y, pw, ph, 0.f, 1.f, 180, 180, 255, "m", mref); y += ph + 4;
    draw_plot(rend, fs, s->hist_e, s->head, s->count, px, y, pw, ph, -2.f, 0.f, 255, 180, 100, "e", eref); y += ph + 4;
    draw_plot(rend, fs, s->hist_chi, s->head, s->count, px, y, pw, ph, 0.f, fmaxf(1.f, (float)(st->chi * 1.5)), 100, 255, 180, "\xcf\x87", NAN); y += ph + 4;
    draw_plot(rend, fs, s->hist_cv, s->head, s->count, px, y, pw, ph, 0.f, fmaxf(1.f, (float)(st->cv * 1.5)), 255, 120, 180, "Cv", NAN); y += ph + 6;

    set_color(rend, 50, 50, 75, 255); hline(rend, SIM_W + 4, WIN_W - 4, y); y += 5;
    static const char *keys[] = {
        "scroll / \xe2\x86\x91\xe2\x86\x93   temperature (shift: fine)",
        "+ / -        sweeps per frame",
        "W            Wolff / Metropolis",
        "M            Ising / Potts      [ ]  Potts q",
        "C            colour mode        space  pause",
        "H / K        hot / cold start   R  reset stats",
        "P            screenshot (BMP)   V  record (ffmpeg)",
        "Q / Esc      quit",
    };
    draw_text(rend, fs, "Controls", px0, y, 160, 160, 100); y += 13;
    for (size_t k = 0; k < sizeof keys / sizeof *keys; k++) {
        draw_text(rend, fs, keys[k], px0, y, 120, 120, 140); y += 12;
    }

    if (recording) {
        if ((SDL_GetTicks() / 500) % 2 == 0) {
            set_color(rend, 220, 30, 30, 255);
            SDL_Rect dot = {WIN_W - 22, 8, 14, 14};
            SDL_RenderFillRect(rend, &dot);
        }
        snprintf(buf, sizeof buf, "REC %ds", rec_frames / 60);
        draw_text(rend, fs, buf, WIN_W - 80, 8, 220, 30, 30);
    }
}

/* ── GUI main loop ───────────────────────────────────────────────────────── */

static void switch_model(Sim *s, Model m, int newq) {
    double ratio = s->T / model_tc(model, q);
    model = m;
    q     = newq;
    build_potts_palette();
    lattice_hot();
    sim_set_T(s, ratio * model_tc(model, q));
    s->total_sweeps = 0;
}

static int run_gui(const Opts *o) {
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL_Init: %s\n(no display? try --batch)\n", SDL_GetError());
        return 1;
    }
    if (TTF_Init() < 0) { fprintf(stderr, "TTF_Init: %s\n", TTF_GetError()); return 1; }

    SDL_Window *win = SDL_CreateWindow("Ising / Potts - parallel Monte Carlo",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIN_W, WIN_H,
        SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);
    SDL_Renderer *rend = win ? SDL_CreateRenderer(win, -1,
        SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC) : NULL;
    if (!win || !rend) { fprintf(stderr, "Window/Renderer: %s\n", SDL_GetError()); return 1; }
    SDL_RenderSetLogicalSize(rend, WIN_W, WIN_H);
    SDL_SetRenderDrawBlendMode(rend, SDL_BLENDMODE_BLEND);

    SDL_Texture *tex = SDL_CreateTexture(rend, SDL_PIXELFORMAT_BGRA32,
                                         SDL_TEXTUREACCESS_STREAMING, L, L);
    if (!tex) { fprintf(stderr, "Texture: %s\n", SDL_GetError()); return 1; }

    TTF_Font *fl = open_font(15), *fs = open_font(12);
    if (!fl || !fs)
        fprintf(stderr, "warning: no TTF font found; panel text disabled "
                        "(set ISING_FONT=/path/to/font.ttf)\n");

    static Sim sim;
    memset(&sim, 0, sizeof sim);
    sim.sweeps_per_frame = 2;
    sim.use_wolff = o->use_wolff;
    sim_set_T(&sim, o->T_set ? o->T : 1.5 * model_tc(model, q));

    FILE *ffmpeg = NULL; uint8_t *cap = NULL; int recording = 0, rec_frames = 0;
    int running = 1, frames = 0;
    Uint32 last = SDL_GetTicks();

    while (running) {
        SDL_Event ev;
        while (SDL_PollEvent(&ev)) {
            if (ev.type == SDL_QUIT) running = 0;
            else if (ev.type == SDL_MOUSEWHEEL) {
                double step = model_tc(model, q) / 50.0;
                if (SDL_GetModState() & KMOD_SHIFT) step /= 10.0;
                sim_set_T(&sim, sim.T + (ev.wheel.y > 0 ? step : -step));
            } else if (ev.type == SDL_KEYDOWN) {
                double step = model_tc(model, q) / 50.0;
                if (ev.key.keysym.mod & KMOD_SHIFT) step /= 10.0;
                switch (ev.key.keysym.sym) {
                case SDLK_ESCAPE: case SDLK_q: running = 0; break;
                case SDLK_UP:   sim_set_T(&sim, sim.T + step); break;
                case SDLK_DOWN: sim_set_T(&sim, sim.T - step); break;
                case SDLK_PLUS: case SDLK_EQUALS: case SDLK_KP_PLUS:
                    if (sim.sweeps_per_frame < 200) sim.sweeps_per_frame++;
                    break;
                case SDLK_MINUS: case SDLK_KP_MINUS:
                    if (sim.sweeps_per_frame > 1) sim.sweeps_per_frame--;
                    break;
                case SDLK_SPACE: sim.paused ^= 1; break;
                case SDLK_w: sim.use_wolff ^= 1; sim_reset_hist(&sim); break;
                case SDLK_c: g_color = (ColorMode)((g_color + 1) % COLOR_COUNT); break;
                case SDLK_m: switch_model(&sim, model == MODEL_ISING ? MODEL_POTTS : MODEL_ISING, q); break;
                case SDLK_LEFTBRACKET:  if (model == MODEL_POTTS && q > 2)    switch_model(&sim, model, q - 1); break;
                case SDLK_RIGHTBRACKET: if (model == MODEL_POTTS && q < MAXQ) switch_model(&sim, model, q + 1); break;
                case SDLK_h: lattice_hot();  sim_reset_hist(&sim); break;
                case SDLK_k: lattice_cold(); sim_reset_hist(&sim); break;
                case SDLK_r: sim_reset_hist(&sim); break;
                case SDLK_p: {
                    char fn[64]; time_t t0 = time(NULL);
                    strftime(fn, sizeof fn, "ising_%Y%m%d_%H%M%S.bmp", localtime(&t0));
                    SDL_Surface *sh = SDL_CreateRGBSurfaceWithFormat(0, WIN_W, WIN_H, 32, SDL_PIXELFORMAT_BGRA32);
                    if (sh && SDL_RenderReadPixels(rend, NULL, SDL_PIXELFORMAT_BGRA32, sh->pixels, sh->pitch) == 0) {
                        SDL_SaveBMP(sh, fn); SDL_Log("Saved %s", fn);
                    }
                    if (sh) SDL_FreeSurface(sh);
                    break;
                }
                case SDLK_v:
                    if (!recording) {
                        char fn[64], cmd[320]; time_t t0 = time(NULL);
                        strftime(fn, sizeof fn, "ising_%Y%m%d_%H%M%S.mp4", localtime(&t0));
                        snprintf(cmd, sizeof cmd,
                            "ffmpeg -loglevel error -y -f rawvideo -pixel_format bgra "
                            "-video_size %dx%d -framerate 60 -i pipe:0 "
                            "-c:v libx264 -preset fast -crf 18 -pix_fmt yuv420p %s", WIN_W, WIN_H, fn);
                        ffmpeg = popen(cmd, "w");
                        if (ffmpeg) {
                            cap = malloc((size_t)WIN_W * WIN_H * 4);
                            recording = 1; rec_frames = 0;
                            SDL_Log("Recording: %s", fn);
                        }
                    } else {
                        recording = 0;
                        free(cap); cap = NULL;
                        pclose(ffmpeg); ffmpeg = NULL;
                        SDL_Log("Recording stopped (%d frames)", rec_frames);
                    }
                    break;
                default: break;
                }
            }
        }
        if (!running) break;

        if (!sim.paused) {
            if (sim.use_wolff) wolff_calibrate();   /* fixed for this frame */
            for (int k = 0; k < sim.sweeps_per_frame; k++) {
                if (sim.use_wolff) wolff_sweep(&sim.avg_cluster);
                else               metropolis_sweep();
                sim.total_sweeps++;
            }
            sim.cur = measure();
            sim_record(&sim);
        }

        dispatch(CMD_PAINT, 0);
        SDL_UpdateTexture(tex, NULL, pixels, L * 4);

        set_color(rend, 5, 5, 12, 255);
        SDL_RenderClear(rend);
        SDL_Rect dst = {0, 0, SIM_W, WIN_H};
        SDL_RenderCopy(rend, tex, NULL, &dst);

        char tbuf[96];
        snprintf(tbuf, sizeof tbuf, "%s  L=%d  T=%.3f  T/Tc=%.3f", model == MODEL_ISING ? "Ising" : "Potts",
                 L, sim.T, sim.T / model_tc(model, q));
        draw_text(rend, fs, tbuf, 6, 6, 180, 180, 100);
        Uint32 now = SDL_GetTicks();
        snprintf(tbuf, sizeof tbuf, "%.0f fps", 1000.0 / (double)(now - last + 1));
        last = now;
        draw_text(rend, fs, tbuf, SIM_W - 60, WIN_H - 18, 80, 80, 100);

        draw_panel(rend, fl, fs, &sim, recording, rec_frames);
        SDL_RenderPresent(rend);

        if (recording && ffmpeg && cap) {
            SDL_RenderReadPixels(rend, NULL, SDL_PIXELFORMAT_BGRA32, cap, WIN_W * 4);
            fwrite(cap, 4, (size_t)WIN_W * WIN_H, ffmpeg);
            rec_frames++;
        }
        if (o->frames > 0 && ++frames >= o->frames) running = 0;
    }

    if (ffmpeg) pclose(ffmpeg);
    free(cap);
    if (fl) TTF_CloseFont(fl);
    if (fs) TTF_CloseFont(fs);
    TTF_Quit();
    SDL_DestroyTexture(tex);
    SDL_DestroyRenderer(rend);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}

/* ── command line ────────────────────────────────────────────────────────── */

static void usage(FILE *f) {
    fprintf(f,
"Usage: %s [options]\n"
"  -L, --size N         lattice side, even (default 512)\n"
"  -j, --threads N      worker threads (default: online CPUs)\n"
"  -m, --model NAME     ising | potts (default ising)\n"
"  -q, --states Q       Potts states, 2..%d (default 3)\n"
"  -T, --temp T         temperature in J/k_B (default 1.5 Tc)\n"
"  -w, --wolff          use the Wolff cluster algorithm\n"
"  -s, --seed S         RNG seed (default: time)\n"
"  -b, --batch          headless: thermalise, measure, print, exit\n"
"  -n, --sweeps N       measurement sweeps in batch mode (default 10000)\n"
"  -t, --therm N        thermalisation sweeps in batch mode (default 1000)\n"
"      --scan A:B:K     batch: K temperatures from A to B inclusive\n"
"      --csv            batch: CSV output\n"
"      --frames N       GUI: quit after N frames\n"
"  -h, --help           this text\n"
"  -V, --version        version\n", PROG, MAXQ);
}

static bool parse_opts(int argc, char **argv, Opts *o) {
    memset(o, 0, sizeof *o);
    o->L = 512; o->threads = 0; o->model = MODEL_ISING; o->q = 3;
    o->sweeps = 10000; o->therm = 1000; o->seed = (uint64_t)time(NULL) ^ (uint64_t)getpid();

    static const struct option lo[] = {
        {"size", 1, 0, 'L'}, {"threads", 1, 0, 'j'}, {"model", 1, 0, 'm'},
        {"states", 1, 0, 'q'}, {"temp", 1, 0, 'T'}, {"wolff", 0, 0, 'w'},
        {"seed", 1, 0, 's'}, {"batch", 0, 0, 'b'}, {"sweeps", 1, 0, 'n'},
        {"therm", 1, 0, 't'}, {"scan", 1, 0, 1000}, {"csv", 0, 0, 1001},
        {"frames", 1, 0, 1002}, {"help", 0, 0, 'h'}, {"version", 0, 0, 'V'}, {0, 0, 0, 0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "L:j:m:q:T:ws:bn:t:hV", lo, NULL)) != -1) {
        switch (c) {
        case 'L': o->L = atoi(optarg); break;
        case 'j': o->threads = atoi(optarg); break;
        case 'm':
            if (!strcmp(optarg, "ising")) o->model = MODEL_ISING;
            else if (!strcmp(optarg, "potts")) o->model = MODEL_POTTS;
            else { fprintf(stderr, "unknown model '%s'\n", optarg); return false; }
            break;
        case 'q': o->q = atoi(optarg); break;
        case 'T': o->T = atof(optarg); o->T_set = true; break;
        case 'w': o->use_wolff = 1; break;
        case 's': o->seed = strtoull(optarg, NULL, 0); break;
        case 'b': o->batch = 1; break;
        case 'n': o->sweeps = atol(optarg); break;
        case 't': o->therm = atol(optarg); break;
        case 1000:
            if (sscanf(optarg, "%lf:%lf:%d", &o->scan_lo, &o->scan_hi, &o->scan_n) != 3 || o->scan_n < 1) {
                fprintf(stderr, "--scan expects A:B:K\n"); return false;
            }
            break;
        case 1001: o->csv = 1; break;
        case 1002: o->frames = atoi(optarg); break;
        case 'h': usage(stdout); exit(0);
        case 'V': printf("%s %s\n", PROG, VERSION); exit(0);
        default: usage(stderr); return false;
        }
    }
    if (o->L < 4 || o->L & 1) { fprintf(stderr, "L must be even and >= 4\n"); return false; }
    if (o->q < 2 || o->q > MAXQ) { fprintf(stderr, "q must be in 2..%d\n", MAXQ); return false; }
    if (o->threads <= 0) {
        long n = sysconf(_SC_NPROCESSORS_ONLN);
        o->threads = n > 0 ? (int)n : 4;
    }
    if (o->threads > MAXTHREADS) o->threads = MAXTHREADS;
    if (o->threads > o->L / 2)   o->threads = o->L / 2;
    if (o->sweeps < 1 || o->therm < 0) { fprintf(stderr, "bad sweep counts\n"); return false; }
    if (o->batch && !o->T_set && o->scan_n == 0) o->T = model_tc(o->model, o->q);
    return true;
}

int main(int argc, char **argv) {
    Opts o;
    if (!parse_opts(argc, argv, &o)) return 2;

    model = o.model;
    q     = o.q;
    rng_seed(&main_rng, o.seed);
    lattice_alloc(o.L);
    build_potts_palette();
    lattice_hot();
    pool_start(o.threads, o.seed * 0x2545F4914F6CDD1DULL + 1);

    int rc = o.batch ? run_batch(&o) : run_gui(&o);

    pool_stop();
    free(spin); free(pixels); free(wstack); free(wmark);
    return rc;
}

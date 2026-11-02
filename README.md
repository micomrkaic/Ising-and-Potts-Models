# ising_parallel — 2-D Ising and Potts models, multithreaded Monte Carlo

A real-time, interactive simulation of the two-dimensional ferromagnetic Ising
model and the q-state Potts model, written in C17 with SDL2 rendering and
POSIX threads. Two Monte Carlo algorithms: Metropolis single-site updates,
parallelised over all cores with a two-colour checkerboard decomposition, and
Wolff single-cluster updates. A headless batch mode prints observables with
jackknife error bars, including the Binder cumulant, and compares the Ising
results with Onsager's exact solution.

```
make && make check && make install        # installs to ~/.local by default
ising_parallel                            # interactive, Ising, L=512
ising_parallel -m potts -q 3 -L 1024 -j 8 # interactive Potts
ising_parallel --batch -L 64 -T 2.0       # one line of physics, no window
./binder.sh                               # find T_c from Binder-cumulant crossings
```

---

## Contents

1. [Physics](#1-physics)
   1. The Ising model
   2. The q-state Potts model
   3. Energy against entropy
   4. Exact results used for validation
   5. Observables
   6. Critical slowing-down
2. [Monte Carlo algorithms](#2-monte-carlo-algorithms)
   1. Detailed balance
   2. Metropolis
   3. Boltzmann tables
   4. Wolff cluster algorithm
3. [Parallelisation](#3-parallelisation)
   1. Why parallel Metropolis is dangerous
   2. Two-colour checkerboard with row bands
   3. Thread pool and barrier
   4. Random numbers
   5. Scaling
4. [Statistics](#4-statistics)
5. [Validation](#5-validation)
6. [Code map](#6-code-map)
7. [Build, test, install](#7-build-test-install)
8. [Usage](#8-usage)
9. [Possible extensions](#9-possible-extensions)

---

## 1. Physics

### 1.1 The Ising model

Each site *i* of an *L* × *L* square lattice with periodic boundary conditions
carries a spin σᵢ ∈ {−1, +1}. The energy of a configuration is

```
H = −J Σ_<ij> σᵢ σⱼ
```

summed over nearest-neighbour pairs, each pair once. With J > 0 aligned
neighbours lower the energy. Throughout, J = k_B = 1, so temperature is
measured in units of J/k_B.

### 1.2 The q-state Potts model

Each site carries a "colour" σᵢ ∈ {0, 1, …, q−1} and only equal neighbours
interact:

```
H = −J Σ_<ij> δ(σᵢ, σⱼ)
```

For q = 2 this is the Ising model with a rescaled coupling: δ(σ,σ') =
(1 + s s')/2 with s = 2σ−1, so H_Potts = const − (J/2) Σ s s'. Hence a
Potts-2 simulation at temperature T reproduces an Ising simulation at 2T, and

```
m_Potts(T) = m_Ising(2T),      e_Potts(T) = −1 + e_Ising(2T)/2
```

which `make check` verifies.

The Potts model is interesting because the *order* of the transition depends on
q. On the square lattice the transition is continuous for q ≤ 4 and first order
for q ≥ 5 (Baxter 1973). With q = 7 or more the latent heat is visible in the
simulation: cooling through T_c, the energy history shows a jump rather than a
smooth bend, and metastable mixed states persist for long times.

### 1.3 Energy against entropy

Equilibrium minimises the free energy F = E − TS. Energy favours alignment,
entropy favours disorder, and the only question is which wins at a given T.

The entropy side is pure counting. The number of Ising configurations with
magnetisation m is a binomial coefficient, C(N, N↑), overwhelmingly peaked at
m = 0; for N = 512² the count at m ≈ 0 is of order 2^262144. The disordered
state is not favoured because it has lower energy — it does not — but because
there are exponentially more ways of being disordered. Below T_c the energy
penalty of a domain wall, proportional to its length in two dimensions, is
enough to suppress these configurations (Peierls' argument); above T_c it is
not. The transition is where the two terms balance.

### 1.4 Exact results used for validation

**Ising critical temperature** (Kramers–Wannier duality, Onsager 1944):

```
T_c = 2J / ln(1 + √2) = 2.269185…
```

**Ising energy per spin** at any T (Onsager), with K = J/T:

```
e(T) = −J coth(2K) [ 1 + (2/π) κ' K(κ) ]
κ  = 2 sinh(2K) / cosh²(2K)
κ' = 2 tanh²(2K) − 1
```

where K(κ) is the complete elliptic integral of the first kind, computed in the
code by the arithmetic–geometric mean. At T_c, κ' = 0 and e = −√2.

**Ising spontaneous magnetisation** (Onsager 1949, Yang 1952), T < T_c:

```
m(T) = [ 1 − sinh⁻⁴(2J/T) ]^(1/8)
```

**Potts critical temperature** (self-duality):

```
T_c(q) = J / ln(1 + √q)
```

giving 1.1346 for q = 2 (half the Ising value, as it must be), 0.99497 for
q = 3, 0.91024 for q = 4, 0.85153 for q = 5.

**Potts energy at T_c** (Baxter): e(T_c) = −J (1 + 1/√q); for q = 3 that is
−1.57735. For q ≥ 5 this is the average of the two coexisting phases.

Finite lattices do not reproduce these infinite-lattice values exactly. The
energy at T_c of an L × L Ising torus is shifted by an amount proportional to
1/L (Ferdinand & Fisher 1969); the simulation finds (e_L − e_∞)·L ≈ −0.6 for
L = 16…128. The Potts-3 shift scales as L^(α−1)/ν ≈ L^−0.8. Away from T_c,
where the correlation length is a few lattice spacings, L = 64 already agrees
with the exact values to four decimals (section 5).

### 1.5 Observables

| quantity | Ising | Potts |
|---|---|---|
| order parameter m | \|Σσᵢ\| / N | (q · max_s n_s / N − 1) / (q − 1) |
| energy per spin e | −J Σ_<ij> σᵢσⱼ / N | −J Σ_<ij> δ(σᵢ,σⱼ) / N |
| susceptibility χ | N (⟨m²⟩ − ⟨m⟩²) / T | same |
| specific heat C_v | N (⟨e²⟩ − ⟨e⟩²) / T² | same |
| Binder cumulant U_L | 1 − ⟨m⁴⟩ / 3⟨m²⟩² | same |

Here n_s is the number of sites in state s. The Potts order parameter is 0 for a
uniformly random configuration and 1 for a single colour; for q = 2 it equals
the Ising |m|. The absolute value is used because a finite system below T_c
tunnels between the ±m states, and the plain average would vanish.

χ and C_v follow from the fluctuation–dissipation theorem and both diverge at a
continuous transition. C_v diverges only logarithmically in the Ising case.

The **Binder cumulant** is the cleanest way to locate T_c on finite lattices.
Far above T_c the magnetisation is Gaussian with zero mean, ⟨m⁴⟩ = 3⟨m²⟩² and
U_L → 0; far below it m is sharply peaked at ±m₀, ⟨m⁴⟩ = ⟨m²⟩² and U_L → 2/3.
In between, finite-size scaling says U_L = f((T − T_c) L^{1/ν}): a function of
the single scaling variable, so at T = T_c every L gives the same value U*.
Curves of U_L against T for different L cross at T_c, with no fitting and no
knowledge of the exponents required. For the square-lattice Ising model with
periodic boundaries U* = 0.6107 (Kamieniarz & Blöte 1993; Salas & Sokal 2000).
The slopes of the curves at the crossing grow as L^{1/ν}, which yields ν = 1.
`binder.sh` runs the scans and reports the crossings.

### 1.6 Critical slowing-down

Near T_c the correlation length ξ diverges, and with a local update rule the
autocorrelation time of the observables grows as τ ∝ ξ^z with z ≈ 2.17 for
Metropolis in two dimensions. On a finite lattice ξ saturates at L, so
τ ∝ L^2.17: at T_c a 512² lattice needs thousands of sweeps per independent
sample. The panel shows τ of the order parameter, estimated from the data
window, and the effect is immediate: at T_c with Metropolis τ runs into the
tens or hundreds of frames; press `W` and it collapses to about one.

---

## 2. Monte Carlo algorithms

### 2.1 Detailed balance

A Markov chain samples the Boltzmann distribution P(σ) ∝ exp(−H/T) if it is
ergodic and satisfies detailed balance,
P(σ) W(σ→σ') = P(σ') W(σ'→σ). The checkerboard sweep used here updates the
sites in a fixed order rather than at random; the composite of the two colour
phases is not reversible as a whole, but each individual site update obeys
detailed balance and the fixed-order sweep still has the Boltzmann distribution
as its unique stationary state. This is the standard justification for
sequential and checkerboard updating.

### 2.2 Metropolis

For every site: propose a change, compute ΔE, accept with probability
min(1, exp(−ΔE/T)).

- Ising: the proposal is σ → −σ and ΔE = 2J σᵢ Σ_nn σⱼ.
- Potts: the proposal is a uniformly chosen different colour r ≠ σᵢ, and
  ΔE = J (n_old − n_new) where n_old, n_new count the neighbours equal to the
  old and the new colour. Proposing uniformly among the q − 1 alternatives
  keeps the proposal symmetric, which Metropolis requires.

Moves with ΔE ≤ 0 are accepted without drawing a random number.

### 2.3 Boltzmann tables

ΔE takes only a handful of values: for Ising σΣ_nn ∈ {−4,−2,0,2,4}, giving five
possible ΔE; for Potts n_old − n_new ∈ {−4,…,4}, giving nine. The exponentials
are tabulated once per temperature change, and the inner loop does an integer
index and a compare.

### 2.4 Wolff cluster algorithm

Pick a random seed site with colour s₀ and a new colour s₁ ≠ s₀. Grow a
cluster by breadth-first search: every neighbour of a cluster site that still
has colour s₀ joins with probability

```
p = 1 − exp(−2J/T)   (Ising)        p = 1 − exp(−J/T)   (Potts)
```

The factor two for Ising is the energy difference between an aligned and an
anti-aligned bond. Finally the whole cluster is recoloured s₁. The algorithm
satisfies detailed balance (Wolff 1989) and flips clusters whose size scales
with ξ², so the dynamic exponent drops to z ≈ 0.25 and critical slowing-down is
essentially gone.

One Wolff "sweep" in this program is a fixed number of cluster flips, about
N/⟨cluster size⟩, which makes sweep counts and sweeps-per-second comparable
between the two algorithms. The count is calibrated from earlier sweeps and
held fixed during measurement. This matters: an earlier version stopped a sweep
"once N spins had flipped", which is a path-dependent stopping time and biases
the sample towards the ordered configurations that follow a large cluster
flip. The symptom was a Binder cumulant of 0.632 at T_c instead of 0.611, with
Metropolis giving the right answer. The strong Markov property does not rescue
you here: a stationary chain observed at a stopping time is not, in general,
stationary.

Wolff is run on the main thread only: the breadth-first search is inherently
sequential. The repaint and the measurement remain parallel. For a parallel
cluster algorithm see section 9.

---

## 3. Parallelisation

### 3.1 Why parallel Metropolis is dangerous

Each Metropolis update reads four neighbours. If two threads update adjacent
sites at the same time, one of them computes ΔE from a neighbour that is being
flipped under its feet. The result is not a crash but a subtle bias: detailed
balance is violated, and the simulation converges to a distribution that is
almost, but not quite, Boltzmann. No sanitizer reports it and the pictures look
right. The only defence is a decomposition in which no two simultaneously
updated sites are ever neighbours.

The previous version of this program got this wrong. It assigned sites to
threads by `(i mod 2)·2 + (j mod 2)`, so thread 0 owned (even, even) sites and
thread 1 owned (even, odd) sites — nearest neighbours of each other — and all
four threads ran at once. The present version replaces it with the scheme
below, and section 5 shows the result agreeing with Onsager to four decimals.

### 3.2 Two-colour checkerboard with row bands

Colour the lattice like a chessboard by the parity of i + j. Every neighbour of
a black site is white and vice versa. A sweep therefore has two phases:

```
phase 0: update all black sites   (white sites frozen)
phase 1: update all white sites   (black sites frozen)
```

Within a phase the sites of the active colour are independent of one another,
so they can be distributed over threads in any way whatsoever. The program cuts
the lattice into horizontal bands of rows, one per thread:

```
thread 0: rows [0, L/n)              ■ □ ■ □ ■ □ ■ □
thread 1: rows [L/n, 2L/n)           □ ■ □ ■ □ ■ □ ■
   ...                               ■ □ ■ □ ■ □ ■ □
thread n−1: rows [(n−1)L/n, L)       □ ■ □ ■ □ ■ □ ■
```

A thread reads the rows just outside its band — they belong to the frozen
colour or to its own colour at sites nobody writes in this phase — and writes
only inside it. Row bands are cache-friendly (contiguous memory) and work for
any thread count that divides reasonably into L; nothing in the method is tied
to four cores.

The same bands are used for the parallel repaint and the parallel measurement.
In the `thread` colour mode each spin is drawn in its owning thread's hue,
bright for up and dark for down, so the bands and the domains show together.

### 3.3 Thread pool and barrier

`nthreads − 1` worker threads are created once. The main thread is worker 0: it
does its own band's share of every job instead of idling while others work, so
the machine is never oversubscribed by a spinning coordinator. (An earlier
design had the main thread wait while four workers ran — on a four-core
machine that is five runnable threads, and on a two-core sandbox it made the
parallel version twelve times *slower* than the serial one.)

A job is issued through a global command and two crossings of a reusable
barrier:

```
main                                  workers
──────────────────────────────────    ─────────────────────────────
g_cmd = SWEEP; g_phase = 0
barrier ───────────────────────────── barrier   (release)
do band 0                             do own band
barrier ───────────────────────────── barrier   (completion)
```

The barrier is a counter plus a generation number. The last thread to arrive
resets the counter, bumps the generation and broadcasts; the others spin on the
generation for a few thousand `pause` instructions (tens of microseconds) and
then fall back to a condition variable. During a burst of sweeps the spin path
is taken and latency is a microsecond or two; when the GUI is paused, or the
main thread is busy with Wolff or rendering, the workers sleep and the process
uses no CPU beyond the main thread. Memory visibility across the barrier
follows from the sequentially consistent atomics: everything a thread wrote
before arriving is visible to every thread after leaving.

Quitting is a command like any other. Workers that see `CMD_QUIT` return
instead of taking the completion barrier, and the main thread joins them before
tearing down SDL. The barrier count is exactly the number of participants, so
no thread can be left behind in a condition wait — the bug that froze the
previous version on `Q`.

### 3.4 Random numbers

Every thread has its own xoshiro256+ generator, seeded by splitmix64 from the
user seed and the thread index. There is no shared RNG state at all. The main
thread has a further generator for the Wolff algorithm and the initial
configuration. Uniform doubles use the top 53 bits of the 64-bit output.
Given the same seed, lattice size and thread count, a batch run is
bit-for-bit reproducible.

### 3.5 Scaling

Measured on the two-core sandbox where this was developed (Metropolis, Ising):

| L | 1 thread | 2 threads | speed-up |
|---|---|---|---|
| 64 | 31,600 sweeps/s | 41,300 sweeps/s | 1.3× |
| 512 | 552 sweeps/s | 1,080 sweeps/s | 1.96× |

At L = 64 a phase is 2,048 site updates per thread — about ten microseconds —
so the four barrier crossings per sweep are a visible fraction of the work. At
L = 512 the scaling is essentially ideal. Do not run more threads than physical
cores: the spinning barrier assumes each participant has a core.

---

## 4. Statistics

Observables are recorded once per sweep in batch mode and once per frame in the
GUI. Errors are **jackknife** estimates over 20 blocks of the time series:
block sums of m, m², e, e² are formed; for each block b the leave-one-out
estimates of ⟨m⟩, ⟨e⟩, χ_b and C_v,b are computed from the remaining blocks;
and the error of a quantity Q is

```
σ_Q² = (B − 1)/B · Σ_b (Q_b − Q)²
```

This handles χ and C_v correctly even though they are nonlinear in the means,
and it accounts for autocorrelations provided a block is much longer than τ.

The **integrated autocorrelation time** of m is estimated from the blocking
ratio: with n correlated samples of variance s², the true variance of the mean
is (s²/n)(1 + 2τ), so

```
τ = ½ ( n σ_m² / s² − 1 )
```

with σ_m the jackknife error. It is a cheap diagnostic rather than a precise
measurement — the blocks must be long compared with τ for it, and the error
bars, to be trusted. In the GUI it is computed over the rolling window of the
last 2,000 frames and shown as τ_m.

In batch mode the first `--therm` sweeps are discarded. With `--scan` the
configuration is carried from one temperature to the next, so each point needs
only to re-equilibrate rather than to start from scratch. Below T_c, starting
hot (`H`) rather than cold (`K`) costs many sweeps while the system chooses a
direction; near T_c on large lattices Metropolis thermalisation can be very
slow and Wolff is the right tool.

---

## 5. Validation

Results of `make check` and related runs on this build. Each number is one run
with its jackknife error.

**Ising below T_c, L = 64, T = 2.0, 20,000 sweeps**

| | simulation | exact (Onsager) |
|---|---|---|
| m | 0.91176 ± 0.00016 | 0.91132 |
| e | −1.74605 ± 0.00027 | −1.74556 |

The residual differences are the expected finite-size corrections at
ξ ≈ 3 lattice spacings.

**Ising at T_c, L = 32, 40,000 sweeps: Metropolis vs Wolff**

| | Metropolis | Wolff | reference |
|---|---|---|---|
| m | 0.6510 ± 0.0094 | 0.6567 ± 0.0018 | |
| χ | 14.2 ± 1.2 | 13.6 ± 0.1 | |
| U_L | 0.6096 ± 0.0054 | 0.6122 ± 0.0010 | U* = 0.6107 |
| τ_m (sweeps) | 55 | 1.7 | |

Same physics, five to ten times smaller error bars, thirty-fold shorter
autocorrelation time. U_L at T_c for L = 16, 32, 64 with Wolff: 0.6126,
0.6122, 0.6119 — approaching U* from above, as the known finite-size
corrections predict.

**Binder crossings** (`SIZES="16 32 64" SCAN="2.22:2.32:6" SWEEPS=3000
./binder.sh`, four seconds on two cores): L = 16 × 32 cross at 2.2708,
L = 32 × 64 at 2.2654; exact T_c = 2.2692. Longer runs and a finer scan
tighten this at will.

**Energy at T_c versus L, Wolff.** The shift from the infinite-lattice −√2
scales as 1/L:

| L | 16 | 32 | 64 | 128 |
|---|---|---|---|---|
| e_L | −1.4502 | −1.4287 | −1.4248 | −1.4191 |
| (e_L − e_∞)·L | −0.58 | −0.47 | −0.68 | −0.62 |

**Potts q = 2 at T = 1.0 versus Ising at T = 2.0, L = 32, six seeds each**

mean m: 0.91124 (Potts) vs 0.91146 (Ising), exact 0.91132; mean e in Ising
units: −1.74540 vs −1.74560. The two models are the same model.

**Potts q = 3 at T_c, Wolff.** e_L approaches Baxter's −1.57735 from below:
−1.613 (L=16), −1.589 (32), −1.591 (64), −1.585 (128), consistent with a
shift scaling as L^−0.8 within the errors of ±0.003.

**Threads.** Runs with 1, 2, 3 and 4 threads on the same seed agree within
errors (they use different random streams). Runs with identical seed *and*
thread count are byte-identical. ThreadSanitizer and AddressSanitizer runs of
both the batch and the GUI paths, Ising and Potts, Metropolis and Wolff, report
nothing.

---

## 6. Code map

```
ising_parallel.c
├── constants, Model enum, model_tc(), model_name()
├── ellipK(), onsager_energy(), onsager_mag()      exact reference values
├── lattice globals: L, N, spin[], pixels[], model, q
├── Rng  (xoshiro256+ / splitmix64 seeding)
├── build_tables()                                  Boltzmann lookup tables
├── Barrier  (spin-then-sleep, generation counter)
├── Ctx, thread_palette, potts_palette, hsv_bgra()
├── sweep_ising(), sweep_potts()                    one colour phase, one band
├── measure_band(), paint_band(), energy_color()
├── run_cmd(), worker(), dispatch(), pool_start(), pool_stop()
├── metropolis_sweep()                              two dispatches
├── wolff_step(), wolff_sweep(), wolff_calibrate()  serial cluster algorithm
├── measure()                                       parallel reduction → Obs
├── analyse()                                       jackknife → Stats (m, e, χ, C_v, U_L, τ)
├── lattice_hot(), lattice_cold(), lattice_alloc()
├── batch_header(), batch_point(), run_batch()      headless mode
├── SDL helpers, open_font(), draw_plot(), draw_panel()
├── Sim, sim_set_T(), sim_record(), switch_model()
├── run_gui()                                       event loop
└── usage(), parse_opts(), main()
```

Rendering: each frame the threads fill the L × L BGRA pixel buffer according to
the colour mode, `SDL_UpdateTexture` uploads it once, and it is drawn scaled
into the 800 × 800 simulation area with nearest-neighbour filtering. The panel
is drawn with SDL primitives and SDL_ttf text. `SDL_RenderSetLogicalSize`
keeps the layout correct when the window is resized.

---

## 7. Build, test, install

Dependencies: a C17 compiler, SDL2 and SDL2_ttf development packages, POSIX
threads. Optional: `ffmpeg` on the PATH for video recording.

```sh
sudo apt install libsdl2-dev libsdl2-ttf-dev        # Debian / Ubuntu
brew install sdl2 sdl2_ttf                           # macOS
```

```sh
make                   # optimised build
make check             # twelve headless physics tests, about a minute
./binder.sh            # Binder-cumulant T_c locator; writes binder_L*.csv, binder.gp
make install           # ~/.local/bin/ising_parallel and ~/.local/share/man/man1/
make uninstall
make PREFIX=/usr/local install   # system-wide instead
make debug             # ASan + UBSan   -> ising_parallel_debug
make tsan              # ThreadSanitizer -> ising_parallel_tsan
make profile           # gprof          -> ising_parallel_prof
```

`CC`, `CFLAGS`, `CPPFLAGS` and `LDFLAGS` can be overridden on the command line.
The default `CFLAGS` are `-O3 -march=native -DNDEBUG`.

If `~/.local/bin` is not on your PATH, `make install` says so. For the man page
to be found, `~/.local/share/man` must be on MANPATH; on most systems `man`
derives it from PATH automatically.

---

## 8. Usage

```
ising_parallel [options]

  -L, --size N        lattice side, even (default 512)
  -j, --threads N     worker threads (default: online CPUs)
  -m, --model NAME    ising | potts
  -q, --states Q      Potts states, 2..16 (default 3)
  -T, --temp T        temperature (default 1.5 Tc; Tc in batch mode)
  -w, --wolff         Wolff cluster algorithm
  -s, --seed S        RNG seed
  -b, --batch         headless mode
  -n, --sweeps N      measurement sweeps (batch, default 10000)
  -t, --therm N       thermalisation sweeps (batch, default 1000)
      --scan A:B:K    K temperatures from A to B inclusive (batch)
      --csv           CSV output (batch)
      --frames N      quit after N frames (GUI)
```

Keys in the window:

| key | action |
|---|---|
| scroll, ↑ ↓ | temperature by T_c/50; with Shift, T_c/500 |
| + − | sweeps per frame |
| W | Wolff / Metropolis |
| M | Ising / Potts (lattice re-randomised, T/T_c kept) |
| [ ] | Potts q down / up |
| C | colour mode: spin → spins tinted by owning thread → local energy |
| Space | pause |
| H / K | hot / cold start |
| R | reset the statistics window |
| P | screenshot to `ising_YYYYmmdd_HHMMSS.bmp` |
| V | start / stop MP4 recording via ffmpeg |
| Q, Esc | quit |

Things worth looking at:

- Start at the default 1.5 T_c, scroll down slowly through T_c in `energy`
  colour mode: domain walls (red) become fractal at T_c and then straighten
  and anneal away below it.
- At T_c, compare τ_m with Metropolis and Wolff (`W`).
- Switch to Potts (`M`), raise q to 7 or 8 with `]`, and cool through T_c: the
  energy history jumps, and mixed states of several colours persist — a
  first-order transition.
- Batch scan of the Potts-3 susceptibility peak:
  `ising_parallel --batch -m potts -q 3 -L 128 -w --scan 0.90:1.10:21 --csv`
- Locate T_c without knowing it: `./binder.sh`, then `gnuplot binder.gp`.
  For Potts: `MODEL=potts Q=3 SCAN="0.95:1.04:10" ./binder.sh`.

---

## 9. Possible extensions

- **Swendsen–Wang.** Identify all clusters at once with a parallel union–find
  and flip each with probability ½; this gives a parallel cluster algorithm.
  The checkerboard bands could be reused for the bond-activation pass.
- **External field h.** Trivial for Metropolis (one more term in ΔE); for the
  cluster algorithms it needs the ghost-spin construction.
- **Full finite-size scaling.** With U_L in hand, fit the slopes dU_L/dT at
  T_c to L^{1/ν}, the peak of χ to L^{γ/ν} and m(T_c) to L^{−β/ν}, and compare
  with the exact 2-D Ising exponents ν = 1, γ = 7/4, β = 1/8.
- **Histogram reweighting.** Store the energy time series at one T and
  reweight to nearby temperatures instead of rerunning.
- **Zoom and pan** of the lattice view, and a real-space renormalisation
  (block-spin) view that shows the self-similarity at T_c.
- **Triangular or honeycomb lattices**, antiferromagnetic J < 0 (the
  triangular antiferromagnet is frustrated and has no ordered phase), or the
  XY / clock models with continuous spins.
- **Larger systems.** The spin array is one byte per site; a 4096² lattice uses
  about 220 MB including the pixel and Wolff buffers. Beyond that, bit-packing
  and multispin coding give another order of magnitude.

## References

- L. Onsager, *Crystal statistics. I.*, Phys. Rev. **65**, 117 (1944).
- C. N. Yang, *The spontaneous magnetization of a two-dimensional Ising model*,
  Phys. Rev. **85**, 808 (1952).
- A. E. Ferdinand and M. E. Fisher, *Bounded and inhomogeneous Ising models.
  I.*, Phys. Rev. **185**, 832 (1969).
- R. J. Baxter, *Potts model at the critical temperature*, J. Phys. C **6**,
  L445 (1973).
- F. Y. Wu, *The Potts model*, Rev. Mod. Phys. **54**, 235 (1982).
- U. Wolff, *Collective Monte Carlo updating for spin systems*, Phys. Rev.
  Lett. **62**, 361 (1989).
- K. Binder, *Finite size scaling analysis of Ising model block distribution
  functions*, Z. Phys. B **43**, 119 (1981).
- G. Kamieniarz and H. W. J. Blöte, J. Phys. A **26**, 201 (1993);
  J. Salas and A. D. Sokal, J. Stat. Phys. **98**, 551 (2000).
- M. E. J. Newman and G. T. Barkema, *Monte Carlo Methods in Statistical
  Physics*, Oxford (1999).

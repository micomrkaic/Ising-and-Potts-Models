#!/bin/sh
# check.sh — headless regression tests for ising_parallel
# Usage: ./check.sh ./ising_parallel
# Each test prints PASS/FAIL; exit status is the number of failures.
set -u
BIN=${1:-./ising_parallel}
fail=0

run() { "$BIN" --batch --csv "$@" | tail -1; }
col() { echo "$1" | cut -d, -f"$2"; }
near() {  # near VALUE TARGET TOL LABEL
    if awk -v v="$1" -v t="$2" -v tol="$3" 'BEGIN{exit !((v-t)<=tol && (t-v)<=tol)}'; then
        printf 'PASS  %-52s %s (target %s ± %s)\n' "$4" "$1" "$2" "$3"
    else
        printf 'FAIL  %-52s %s (target %s ± %s)\n' "$4" "$1" "$2" "$3"; fail=$((fail+1))
    fi
}

# 1. Ising below Tc against Onsager–Yang / Onsager exact values
r=$(run -L 64 -T 2.0 -n 6000 -t 1500 -s 1)
near "$(col "$r" 9)"  0.91132 0.004 "Ising L=64 T=2.0 magnetisation vs exact"
near "$(col "$r" 11)" -1.74556 0.006 "Ising L=64 T=2.0 energy vs exact"

# 2. Ising above Tc: m -> 0, energy vs exact
r=$(run -L 64 -T 3.0 -n 4000 -t 500 -s 2)
near "$(col "$r" 9)"  0.0 0.06 "Ising L=64 T=3.0 magnetisation ~ 0"
near "$(col "$r" 11)" "$(col "$r" 21)" 0.01 "Ising L=64 T=3.0 energy vs exact"

# 3. Wolff reproduces Metropolis (energy at T=2.1)
a=$(run -L 48 -T 2.1 -n 6000 -t 1000 -s 3)
b=$(run -L 48 -T 2.1 -n 6000 -t 500  -s 3 -w)
near "$(col "$b" 11)" "$(col "$a" 11)" 0.006 "Wolff vs Metropolis energy, Ising T=2.1"

# 4. Wolff kills critical slowing-down at Tc
a=$(run -L 32 -T 2.269185 -n 3000 -t 1000 -s 4)
b=$(run -L 32 -T 2.269185 -n 3000 -t 300  -s 4 -w)
if awk -v m="$(col "$a" 19)" -v w="$(col "$b" 19)" 'BEGIN{exit !(w < m/3)}'; then
    printf 'PASS  %-52s metro %s  wolff %s\n' "tau_m(Wolff) < tau_m(Metropolis)/3 at Tc" "$(col "$a" 19)" "$(col "$b" 19)"
else
    printf 'FAIL  %-52s metro %s  wolff %s\n' "tau_m(Wolff) < tau_m(Metropolis)/3 at Tc" "$(col "$a" 19)" "$(col "$b" 19)"; fail=$((fail+1))
fi

# 5. Potts q=2 at T is Ising at 2T:  m equal,  e_P = -1 + e_I/2
a=$(run -L 48 -T 2.0 -n 6000 -t 1000 -s 5)
b=$(run -m potts -q 2 -T 1.0 -L 48 -n 6000 -t 1000 -s 5)
near "$(col "$b" 9)" "$(col "$a" 9)" 0.004 "Potts q=2 magnetisation = Ising"
ep=$(awk -v e="$(col "$a" 11)" 'BEGIN{printf "%.5f", -1+e/2}')
near "$(col "$b" 11)" "$ep" 0.004 "Potts q=2 energy = -1 + e_Ising/2"

# 6. Potts q=3 at Tc: energy approaches -(1+1/sqrt 3) = -1.57735 (finite-size shift ~ -0.03 at L=64)
r=$(run -m potts -q 3 -L 64 -n 4000 -t 500 -s 6 -w)
near "$(col "$r" 11)" -1.5774 0.05 "Potts q=3 L=64 energy at Tc (finite-size)"

# 7. Thread-count independence (different RNG streams, same physics)
a=$(run -L 64 -T 2.0 -n 4000 -t 1000 -s 7 -j 1)
b=$(run -L 64 -T 2.0 -n 4000 -t 1000 -s 7 -j 4)
near "$(col "$b" 9)" "$(col "$a" 9)" 0.004 "1 thread vs 4 threads, magnetisation"

# 8. Bit reproducibility for a fixed seed and thread count
a=$(run -L 32 -T 2.2 -n 400 -t 100 -s 8 -j 2 | cut -d, -f1-22)
b=$(run -L 32 -T 2.2 -n 400 -t 100 -s 8 -j 2 | cut -d, -f1-22)
if [ "$a" = "$b" ]; then printf 'PASS  %-52s\n' "same seed, same threads -> identical output"
else printf 'FAIL  %-52s\n' "same seed, same threads -> identical output"; fail=$((fail+1)); fi

# 9. Binder cumulant at Tc: U_L -> U* = 0.6107 (2-D Ising, periodic), weak L dependence
r=$(run -L 32 -T 2.269185 -n 6000 -t 500 -s 9 -w)
near "$(col "$r" 17)" 0.611 0.02 "Binder cumulant U_L at Tc, L=32, Wolff"

echo "$fail failure(s)"
exit $fail

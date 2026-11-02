#!/bin/sh
# binder.sh — locate T_c from the crossing of Binder cumulants U_L(T)
#
# Runs a temperature scan for several lattice sizes with the Wolff algorithm,
# writes binder_L<L>.csv for each, prints a table of U_L(T), estimates the
# crossing temperature for each pair of consecutive sizes by linear
# interpolation, and writes binder.gp for gnuplot.
#
# Usage:  ./binder.sh [ising_parallel binary]
# Environment overrides:
#   SIZES="16 32 64 128"   lattice sides
#   SCAN="2.20:2.34:8"     temperature range A:B:K
#   SWEEPS=4000  THERM=500 per temperature
#   MODEL=ising  Q=3       e.g. MODEL=potts Q=3 SCAN="0.95:1.04:10"
#   THREADS=                passed as -j if set
set -eu
BIN=${1:-./ising_parallel}
SIZES=${SIZES:-"16 32 64 128"}
SCAN=${SCAN:-"2.20:2.34:8"}
SWEEPS=${SWEEPS:-4000}
THERM=${THERM:-500}
MODEL=${MODEL:-ising}
Q=${Q:-3}
JOPT=${THREADS:+-j $THREADS}

for L in $SIZES; do
    echo "L=$L ..." >&2
    # shellcheck disable=SC2086
    "$BIN" --batch --csv -w -m "$MODEL" -q "$Q" -L "$L" -n "$SWEEPS" -t "$THERM" \
        --scan "$SCAN" $JOPT > "binder_L$L.csv"
done

echo
printf '%-10s' "T"
for L in $SIZES; do printf '%-22s' "U_$L"; done
echo
paste -d' ' $(for L in $SIZES; do echo "binder_L$L.csv"; done) | tail -n +2 |
awk -v n="$(echo $SIZES | wc -w)" '{
    split($1, a, ","); printf "%-10.5f", a[6];
    for (i = 1; i <= n; i++) { split($i, c, ","); printf "%-10.5f ±%-9.5f ", c[17], c[18] }
    print ""
}'

echo
echo "Crossings (linear interpolation between neighbouring temperatures):"
set -- $SIZES
prev=""
for L in $SIZES; do
    if [ -n "$prev" ]; then
        paste -d' ' "binder_L$prev.csv" "binder_L$L.csv" | tail -n +2 |
        awk -v l1="$prev" -v l2="$L" '{
            split($1, a, ","); split($2, b, ",");
            T[NR] = a[6]; d[NR] = a[17] - b[17];
        } END {
            for (i = 1; i < NR; i++)
                if (d[i] * d[i+1] <= 0 && d[i] != d[i+1]) {
                    tc = T[i] + (T[i+1] - T[i]) * d[i] / (d[i] - d[i+1]);
                    printf "  L=%-4s x L=%-4s  T_c ~ %.5f\n", l1, l2, tc
                }
        }'
    fi
    prev=$L
done
TCEX=$(awk -v m="$MODEL" -v q="$Q" 'BEGIN{ if (m=="ising") printf "%.5f", 2/log(1+sqrt(2)); else printf "%.5f", 1/log(1+sqrt(q)) }')
echo "  exact            T_c = $TCEX"

cat > binder.gp <<EOF
# gnuplot binder.gp
set datafile separator ","
set title "Binder cumulant U_L(T), $MODEL model"
set xlabel "T"; set ylabel "U_L = 1 - <m^4>/(3<m^2>^2)"
set key left bottom
set arrow from $TCEX, graph 0 to $TCEX, graph 1 nohead dt 2 lc rgb "gray"
plot \\
$(for L in $SIZES; do printf '  "binder_L%s.csv" using 6:17:18 with yerrorlines title "L = %s", \\\n' "$L" "$L"; done)
  1/0 notitle
pause -1 "Enter to close"
EOF
echo
echo "gnuplot binder.gp   # to plot"

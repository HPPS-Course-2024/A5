#!/bin/sh
#
# Quick benchmarking script.  Not elegant, but good enough to answer
# the exam questions.  We don't expect that the students automate this
# - it's fine to run things by hand.

set -e
THREADS="1 2 4 8"
MAX_THREADS=8

SMALL=10000
NS="500 1000 5000 10000"
RUNS=10

echo Compiling
make

calc() {
    echo "scale=10; $*" | bc
}

bench_nbody() {
    OMP_NUM_THREADS=${1} ./nbody ${2}.particles /dev/null ${RUNS}
}

bench_nbody_bh() {
    OMP_NUM_THREADS=${1} ./nbody-bh ${2}.particles /dev/null ${RUNS}
}

for t in ${THREADS}; do
    N=$(calc "sqrt($t) * $SMALL")
    ./genparticles ${N} ${N}.particles
done

for N in ${NS}; do
    ./genparticles ${N} ${N}.particles
done

echo
echo Measuring strong scaling
echo
echo "Speedups for nbody:"
base=$(bench_nbody 1 ${SMALL})
for t in ${THREADS}; do
    printf "%2d: " ${t}
    d=$(bench_nbody ${t} ${SMALL})
    calc "$base/$d"
done

echo
echo "Speedups for nbody-bh:"
base=$(bench_nbody_bh 1 ${SMALL})
for t in ${THREADS}; do
    printf "%2d: " ${t}
    d=$(bench_nbody_bh ${t} ${SMALL})
    calc "$base/$d"
done

echo
echo Measuring weak scaling
echo
echo "Speedups for nbody:"
d=$(bench_nbody 1 ${SMALL})
base=$(calc "${SMALL}*${SMALL} / $d")
for t in ${THREADS}; do
    printf "%2d: " ${t}
    N=$(calc "sqrt($t) * $SMALL")
    d=$(bench_nbody ${t} ${N})
    calc "(${N}*${N} / $d) / $base"
done

echo
echo Speedups of nbody-bh vs nbody
echo
for N in ${NS}; do
    printf "%5d: " ${N}
    x=$(bench_nbody ${MAX_THREADS} ${N})
    y=$(bench_nbody_bh ${MAX_THREADS} ${N})
    calc "$x/$y"
done

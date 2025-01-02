#!/bin/bash

# Usage: run_nbody.sh <omp threads> <input> <particle output> <warnings output>

export OMP_NUM_THREADS=${1}
echo "Running on ${2} with ${OMP_NUM_THREADS} threads to ${3} and ${4}"
time ./nbody-bash ${2} ${3} ${4}

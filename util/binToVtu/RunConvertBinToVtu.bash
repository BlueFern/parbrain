#!/usr/bin/env bash
# Bash wrapper to run ConvertBinToVtu with OpenMP
# The first argument specifies number of OpenMP threads to run the proram.
# The second and third arguments are optional and specify the CPU affinitiy to run the program.
# The rest of arguments are arguments required by the program itself.


if [ $# -eq 0 ]; then
	echo "$0 <num of threads> [-c <CPU affinity mask> (optional)] <Data directory> <Final time>"
	exit 0
fi;

export OMP_NUM_THREADS=$1
shift

if [ $1 = "-c" ]; then
	MASK=$2
	shift
	shift
	taskset -c $MASK ./ConvertBinToVtu $@
else
	./ConvertBinToVtu $@
fi;



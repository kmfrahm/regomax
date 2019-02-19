#!/bin/bash
# simple compile script for *.cc in this folder
# requires a working C++ compiler, typically g++

CPP="g++"
OMP="-fopenmp"
# use this for the intel compiler, (only tested for an old version of icc)
#CPP="icc"
#OMP="-openmp"
OPT="-DITER_MODE -O3"

# get name of newest *.h file if it exists
LASTH=`ls -rt *.h 2>/dev/null | tail -1`

if [[ -d TEST ]] ; then
    TARGET=TEST
else
    TARGET=.
fi

# compile a *.cc file if it (or any *.h file) is newer 
# than corresponding *.out file 
CPPFILES="*.cc"

for i in $CPPFILES ; do
    OUT=$TARGET/${i%.*}.out
    if [[ $i -nt $OUT || $LASTH -nt $OUT ]] ; then
	echo "$CPP $OPT $i -o $OUT "
	$CPP $OPT $i -o $OUT
	strip -p $OUT
    fi
done


# compile openmp versions for certain progams
CPPFILES="*reduced.cc"

for i in $CPPFILES ; do
    OUT=$TARGET/${i%.*}_omp.out
    if [[ $i -nt $OUT || $LASTH -nt $OUT ]] ; then
	echo "$CPP $OPT $OMP $i -o $OUT "
	$CPP $OPT $OMP $i -o $OUT 
	strip -p $OUT
    fi
done

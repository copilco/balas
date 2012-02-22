#!/bin/sh

# run.sh
# 
# Created by Alejandro de la Calle on 20/02/12.
# 
#

MKLROOT='/opt/intel/mkl'
COMPOSERLIB='/opt/intel/composer_xe_2011_sp1.9.289/'
DYLD_LIBRARY_PATH='/opt/intel/composer_xe_2011_sp1.9.289/compiler/lib'
MKL_NUM_THREADS=8

for((i=1; i<=20; i++ ))
do	
	export DYLD_LIBRARY_PATH
	icpc mainTestImaginaryHankelProp2D.cpp -L${MKLROOT}/lib -I../../include  -I${MKLROOT}/include ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_intel_thread.a ${MKLROOT}/lib/libmkl_core.a ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_intel_thread.a ${MKLROOT}/lib/libmkl_core.a -liomp5 -lpthread -lm -openmp -o mainTestImaginaryHankelProp2D
	nombre=./ImaginaryROMs/$i
	mkdir -p $nombre
	mv mainTestImaginaryHankelProp2D $nombre/
	cp besselZeros.bin $nombre/
	cd $nombre
	./mainTestImaginaryHankelProp2D $i
	cd ../../
done	
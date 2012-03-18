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

	export DYLD_LIBRARY_PATH
	icpc mainTestXUVplusATTO.cpp -L${MKLROOT}/lib -I../../../include  -I${MKLROOT}/include ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_intel_thread.a ${MKLROOT}/lib/libmkl_core.a ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_intel_thread.a ${MKLROOT}/lib/libmkl_core.a -liomp5 -lpthread -lm -openmp -o mainTestXUVplusATTO


for((i=2; i<=20; i++ ))
do	
	nombre=./$i
	mkdir -p $nombre
	cp mainTestXUVplusATTO $nombre/
	cp besselZeros.bin $nombre/
	cp binwave.bin $nombre/
	cd $nombre
	./mainTestXUVplusATTO $i
	cd ../
done	
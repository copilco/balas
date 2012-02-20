MKLROOT=/opt/intel/mkl
COMPOSERLIB=/opt/intel/composer_xe_2011_sp1.9.289/
DYLD_LIBRARY_PATH=/opt/intel/composer_xe_2011_sp1.9.289/compiler/lib
MKL_NUM_THREADS=8


UniHankel: mainTestUnitaryHankelProp2D.cpp clean
	 export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	 icpc mainTestUnitaryHankelProp2D.cpp -L$(MKLROOT)/lib -I../../include  -I$(MKLROOT)/include $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a -liomp5 -lpthread -lm -openmp -o mainTestUnitaryHankelProp2D
	export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	 ./mainTestUnitaryHankelProp2D

.PHONY: clean
clean: 
	-@rm *.txt
	#-@rm *.out

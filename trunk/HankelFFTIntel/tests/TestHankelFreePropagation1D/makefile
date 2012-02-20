MKLROOT=/opt/intel/mkl
COMPOSERLIB=/opt/intel/composer_xe_2011_sp1.9.289/
DYLD_LIBRARY_PATH=/opt/intel/composer_xe_2011_sp1.9.289/compiler/lib/
MKL_NUM_THREADS=8

HankelProp1D: Hankel1D.cpp clean
	export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	icpc Hankel1D.cpp -L$(MKLROOT)/lib -I../../include  -I$(MKLROOT)/include $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a -liomp5 -lpthread -lm -openmp -o Hankel1D
	export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	./Hankel1D

.PHONY: clean
clean: 
	-@rm *.txt
	#-@rm *.out
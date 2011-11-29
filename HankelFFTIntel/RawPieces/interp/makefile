MKLROOT=/opt/intel/mkl
COMPOSERLIB=/opt/intel/composer_xe_2011_sp1.7.246
DYLD_LIBRARY_PATH=/opt/intel/composer_xe_2011_sp1.7.246/compiler/lib
MKL_NUM_THREADS=8
interp: interp.cpp
	export MKLROOT; \
	export DYLD_LIBRARY_PATH; \
	icpc interp.cpp -L$(MKLROOT)/lib -I$(MKLROOT)/include $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a -liomp5 -lpthread -lm
	export DYLD_LIBRARY_PATH; \
	./a.out
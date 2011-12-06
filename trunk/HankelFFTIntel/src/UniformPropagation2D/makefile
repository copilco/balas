MKLROOT=/opt/intel/mkl
COMPOSERLIB=/opt/intel/composerxe-2011.5.209/
DYLD_LIBRARY_PATH=/opt/intel/composerxe-2011.5.209/compiler/lib/
MKL_NUM_THREADS=8

uni: mainTestPropUniform.cpp
	export MKLROOT; \
	export DYLD_LIBRARY_PATH; \
	icpc mainTestPropUniform.cpp -L$(MKLROOT)/lib -I../../include  -I$(MKLROOT)/include $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a -liomp5 -lpthread -lm
	export DYLD_LIBRARY_PATH; \
	./a.out

uni2D: mainTestPropUniform2D.cpp
	export MKLROOT; \
	export DYLD_LIBRARY_PATH; \
	icpc mainTestPropUniform2D.cpp -L$(MKLROOT)/lib -I../../include  -I$(MKLROOT)/include $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a -liomp5 -lpthread -lm
	export DYLD_LIBRARY_PATH; \
	./a.out


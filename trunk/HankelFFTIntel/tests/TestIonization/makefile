MKLROOT=/opt/intel/mkl
COMPOSERLIB=/opt/intel/composer_xe_2011_sp1.9.289/
DYLD_LIBRARY_PATH=/opt/intel/composer_xe_2011_sp1.9.289/compiler/lib
MKL_NUM_THREADS=4

XUVatto: mainTestXUVplusATTO.cpp clean
	 export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	 icpc mainTestXUVplusATTO.cpp -L$(MKLROOT)/lib -I../../include  -I$(MKLROOT)/include $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a -liomp5 -lpthread -lm -openmp -o mainTestXUVplusATTO
	export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	time -p ./a.out mainTestXUVplusATTO

attoCrank: AttoCrank2D.cpp clean
	 export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	 icpc AttoCrank2D.cpp -L$(MKLROOT)/lib -I../../include  -I$(MKLROOT)/include $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a -liomp5 -lpthread -lm -openmp -o AttoCrank2D
	export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	time -p ./AttoCrank2D

.PHONY: clean
clean:
 
	-@rm *.txt
	#-@rm *.out

MKLROOT=/opt/intel/mkl
COMPOSERLIB=/opt/intel/composer_xe_2011_sp1.8.269/
DYLD_LIBRARY_PATH=/opt/intel/composer_xe_2011_sp1.8.269/compiler/lib/
MKL_NUM_THREADS=8

PrecisionHankel: CreateHankelMatrix.cpp clean
	export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	icpc CreateHankelMatrix.cpp -L$(MKLROOT)/lib -I../../include  -I$(MKLROOT)/include $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a -liomp5 -lpthread -lm -o CreateHankelMatrix
	export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	./CreateHankelMatrix

.PHONY: clean
clean: 
	-@rm *.txt
	#-@rm *.out
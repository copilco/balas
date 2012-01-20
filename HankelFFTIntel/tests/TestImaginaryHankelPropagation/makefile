MKLROOT=/opt/intel/mkl
COMPOSERLIB=/opt/intel/composer_xe_2011_sp1.8.269/
DYLD_LIBRARY_PATH=/opt/intel/composer_xe_2011_sp1.8.269/compiler/lib
MKL_NUM_THREADS=8


ImaginaryHankelProp: mainTestImaginaryHankelProp2D.cpp clean
	 export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	 icpc mainTestImaginaryHankelProp2D.cpp -L$(MKLROOT)/lib -I../../include  -I$(MKLROOT)/include $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a -liomp5 -lpthread -lm -o mainTestImaginaryHankelProp2D
	export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	./mainTestImaginaryHankelProp2D

TestFFT: TestFFT.cpp clean
	 export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	 icpc TestFFT.cpp -L$(MKLROOT)/lib -I../../include  -I$(MKLROOT)/include $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a -liomp5 -lpthread -lm -o TestFFT
	export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	./TestFFT

TESTFFT1D: TESTFFT1D.cpp clean
	 export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	 icpc TESTFFT1D.cpp -L$(MKLROOT)/lib -I../../include  -I$(MKLROOT)/include $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a -liomp5 -lpthread -lm -o TESTFFT1D
	export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	./TESTFFT1D

.PHONY: clean
clean: 
	-@rm *.txt
	#-@rm *.out

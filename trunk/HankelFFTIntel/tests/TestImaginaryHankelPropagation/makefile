MKLROOT=/opt/intel/mkl
COMPOSERLIB=/opt/intel/composerxe-2011.5.209/
DYLD_LIBRARY_PATH=/opt/intel/composerxe-2011.5.209/compiler/lib
MKL_NUM_THREADS=8


ImaginaryHankelProp: mainTestImaginaryHankelProp2D.cpp clean
	 export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	 icpc mainTestImaginaryHankelProp2D.cpp -L$(MKLROOT)/lib -I../../include  -I$(MKLROOT)/include $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a -liomp5 -lpthread -lm -o mainTestImaginaryHankelProp2D
	export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	./mainTestImaginaryHankelProp2D

.PHONY: clean
clean: 
	-@rm *.txt
	#-@rm *.out

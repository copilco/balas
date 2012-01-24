MKLROOT=/opt/intel/mkl
COMPOSERLIB=/opt/intel/composer_xe_2011_sp1.8.269/
DYLD_LIBRARY_PATH=/opt/intel/composer_xe_2011_sp1.8.269/compiler/lib
MKL_NUM_THREADS=8


LoopCEP: StrongFieldATI_loopCEP.cpp clean
	 export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	 icpc StrongFieldATI_loopCEP.cpp -L$(MKLROOT)/lib -I../../include  -I$(MKLROOT)/include $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a -liomp5 -lpthread -lm -o StrongFieldATI_loopCEP
	export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	./StrongFieldATI_loopCEP

.PHONY: clean
clean: 
	-@rm *.txt
	#-@rm *.out

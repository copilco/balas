MKLROOT=/opt/intel/mkl
COMPOSERLIB=/opt/intel/composer_xe_2011_sp1.9.289/
DYLD_LIBRARY_PATH=/opt/intel/composer_xe_2011_sp1.9.289/compiler/lib
MKL_NUM_THREADS=12
n=0

multiphoton: ScanMultiPhotonOpemmpH1D.cpp #clean
	 export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	 icpc ScanMultiPhotonOpemmpH1D.cpp -L$(MKLROOT)/lib -I /Users/alexis/Dropbox/BALAS/include  -I$(MKLROOT)/include $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a -liomp5 -lpthread -lm -fopenmp -o multi
	export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	#time -p ./multi $n &


.PHONY: clean
clean: 
	-@rm *.txt
	#-@rm *.out

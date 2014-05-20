MKLROOT=/opt/intel/mkl
COMPOSERLIB=/opt/intel/composer_xe_2011_sp1.9.289/
DYLD_LIBRARY_PATH=/opt/intel/composer_xe_2011_sp1.9.289/compiler/lib
OWNER_ROOT=/Users/alexis/Dropbox/CODE/BALAS
MKL_NUM_THREADS=2


XUVatto:  Hyd_HHG1s.cpp  #clean
	 export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	 icpc  Hyd_HHG1s.cpp -L$(MKLROOT)/lib -I$(OWNER_ROOT)/include -I$(MKLROOT)/include $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a -liomp5 -lpthread -lm -openmp 

clean:
	rm multiphoton

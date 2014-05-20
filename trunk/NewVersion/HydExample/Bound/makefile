MKLROOT=/opt/intel/mkl
COMPOSERLIB=/opt/intel/composer_xe_2011_sp1.9.289/
DYLD_LIBRARY_PATH=/opt/intel/composer_xe_2011_sp1.9.289/compiler/lib
#COMPOSERLIB=/opt/intel/composer_xe_2013.3.171
#DYLD_LIBRARY_PATH=/opt/intel/composer_xe_2013.3.171/compiler/lib
#OWNER=/Users/alexis/Desktop/BALAS/Source
OWNER=/Users/alexis/Dropbox/CODE/BALAS
MKL_NUM_THREADS=8

XUVatto: BoundGround3D.cpp #clean
	 export DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH); \
	 icpc BoundGround3D.cpp -L$(MKLROOT)/lib -I$(OWNER)/include -I$(MKLROOT)/include $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_intel_lp64.a $(MKLROOT)/lib/libmkl_intel_thread.a $(MKLROOT)/lib/libmkl_core.a -liomp5 -lpthread -lm -openmp 

clean:
	rm multiphoton

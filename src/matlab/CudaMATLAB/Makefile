# Define installation location for CUDA and compilation flags compatible
# with the CUDA include files.
INCLUDEDIR  = -I$(CUDAHOME)/include -I$(CUDAHOME)/NVIDIA_GPU_Computing_SDK/C/common/inc -I$(CUDAHOME)/NVIDIA_GPU_Computing_SDK/shared/inc -I/$(TDA)


INCLUDELIB  = -L$(CUDAHOME)/lib64 -lcufft -lcurand -lcublas  -lcudart -Wl,-rpath,$(CUDAHOME)/lib64
CFLAGS      = -fPIC -D_GNU_SOURCE -pthread -fexceptions -arch glnxa64 
COPTIMFLAGS = -O3 -funroll-loops -msse2



# Define installation location for MATLAB.
MEX           = $(MATLAB)/bin/mex
MEXEXT        = $(MATLAB)/bin/mexext

# nvmex is a modified mex script that knows how to handle CUDA .cu files.
NVMEX = $(TDA)/CudaMATLAB/nvmex
NVOPTS = $(TDA)/CudaMATLAB/nvopts.sh

# List the mex files to be built.  The .mex extension will be replaced with the
# appropriate extension for this installation of MATLAB, e.g. .mexglx or
# .mexa64.
MEXFILES = fft2_cuda.mex       \
           fft2_cuda_sp_dp.mex \
           ifft2_cuda.mex      \
           Szeta.mex

all: $(MEXFILES:.mex=$(MEXEXT))

clean:
	rm -f $(MEXFILES:.mex=$(MEXEXT))

.SUFFIXES: .cu .cu_o .mexglx .mexa64 .mexmaci

.c.mexglx:
	$(MEX) CFLAGS='$(CFLAGS)' COPTIMFLAGS='$(COPTIMFLAGS)' $< \
        $(INCLUDEDIR) $(INCLUDELIB)

.cu.mexglx:
	$(NVMEX) -f $(NVOPTS) $< $(INCLUDEDIR) $(INCLUDELIB)

.c.mexa64:
	export ARCH=glnxa64
	$(MEX) ARCH=glnxa64 -v CFLAGS='$(CFLAGS)' COPTIMFLAGS='$(COPTIMFLAGS)' $< \
        $(INCLUDEDIR) $(INCLUDELIB) 

.cu.mexa64:
	export ARCH=glnxa64
	$(NVMEX) -f $(NVOPTS) $< $(INCLUDEDIR) $(INCLUDELIB) 

.c.mexmaci:
	$(MEX) CFLAGS='$(CFLAGS)' COPTIMFLAGS='$(COPTIMFLAGS)' $< \
        $(INCLUDEDIR) $(INCLUDELIB)

.cu.mexmaci:
	$(NVMEX) -f nvopts.sh $< $(INCLUDEDIR) $(INCLUDELIB)

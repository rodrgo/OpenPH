#!/bin/sh
# Dirty hack to automatically get matlab path
MATLAB_PATH=$(matlab -nosplash -nodesktop -r "try; f=matlabroot; fprintf('%s:',f); catch; end; exit" | tail -n1 | awk -F: '{print $1}')
echo "export MATLAB=${MATLAB_PATH}" >> config

# CUDAHOME
CUDA_VERSION=$(nvcc --version | grep -oP '(?<=release )[0-9.]+')
CUDAHOME=$(whereis "cuda-${CUDA_VERSION}" | awk '{print $2}')
if [ -z "$CUDAHOME" ]; then
    CUDAHOME=$(whereis "cuda" | awk '{print $2}')
fi
echo "export CUDAHOME=${CUDAHOME}" >> config

# SM
CUDA_SM=$(python sm.py)
echo "export GPUARCH=${CUDA_SM}" >> config

# CUDAMATLAB
BASEDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
echo "CUDAMATLAB=${BASEDIR}CudaMatlab" >> config

# Get nvmex and nvopts.sh
wget http://developer.download.nvidia.com/compute/cuda/1_1/Matlab_Cuda_1.1.tgz
tar -xvzf Matlab_Cuda_1.1.tgz
mv Matlab_Cuda_1.1/nvmex CudaMatlab
mv Matlab_Cuda_1.1/nvopts.sh CudaMatlab
rm -rf Matlab_Cuda_1.1
rm Matlab_Cuda_1.1.tgz

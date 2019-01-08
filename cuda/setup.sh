
# Dirty hack to automatically get matlab path
MATLAB_PATH=$(matlab -nosplash -nodesktop -r "try; f=matlabroot; fprintf(2,'%s',f); catch; end; exit" > /dev/null 2 | sed -e 's/^[ \t\n\r]*//')

# Dirty hack to get cuda version
CUDA_VERSION=$(nvcc --version | grep -oP '(?<=release )[0-9.]+')

# Get nvmex and nvopts.sh
wget http://developer.download.nvidia.com/compute/cuda/1_1/Matlab_Cuda_1.1.tgz
tar -xvzf Matlab_Cuda_1.1.tgz
mv Matlab_Cuda_1.1/nvmex CudaMatlab
mv Matlab_Cuda_1.1/nvopts.sh CudaMatlab
rm -rf Matlab_Cuda_1.1
rm Matlab_Cuda_1.1.tgz

# Create config file
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "export TDA=${DIR}/cuda" >> config
echo "export MATLAB=${MATLAB_PATH}" >> config
echo "export CUDAHOME=/usr/local/cuda-${CUDA_VERSION}" >> config
echo "export GPUARCH=sm_61" >> config
echo "CUDAMATLAB=${DIR}/cuda/CudaMatlab" >> config


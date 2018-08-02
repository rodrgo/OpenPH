
void __global__ indexShiftDown(int *d_rows, const int m){
  unsigned int xIndex = blockDim.x * blockIdx.x + threadIdx.x;
  if (xIndex < m) d_rows[xIndex] = d_rows[xIndex]-1;
}

void __global__ indexShiftUp(int *d_rows, const int m){
  unsigned int xIndex = blockDim.x * blockIdx.x + threadIdx.x;
  if (xIndex < m) d_rows[xIndex] = d_rows[xIndex]+1;
}

__global__ void zero_vector_float(float *vec, const int n){
  unsigned int xIndex = blockDim.x * blockIdx.x + threadIdx.x;
  if (xIndex < n) vec[xIndex]=0.0f;
}

__global__ void fill(float *vec, int value, int n){
  unsigned int xIndex = blockDim.x * blockIdx.x + threadIdx.x;
  if (xIndex < n) vec[xIndex]=value;
}

__global__ void minusone_vector_int(int *vec, const int n){
  unsigned int xIndex = blockDim.x * blockIdx.x + threadIdx.x;
  if ( xIndex < n ){
    vec[xIndex]=-1;
  }
}

__global__ void zero_vector_int(int *vec, const int n){
  unsigned int xIndex = blockDim.x * blockIdx.x + threadIdx.x;
  if ( xIndex < n ){
    int z=0;
    vec[xIndex]=z;
  }
}

__global__ void one_vector_float(float *vec, const int n){
  unsigned int xIndex = blockDim.x * blockIdx.x + threadIdx.x;
  if (xIndex < n) vec[xIndex]=1.0f;
}

__global__ void one_vector_int(int *vec, const int n){
  unsigned int xIndex = blockDim.x * blockIdx.x + threadIdx.x;
  if (xIndex < n) vec[xIndex]=1;
}


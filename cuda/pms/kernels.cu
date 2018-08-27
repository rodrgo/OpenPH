
// -------------------
// Generic
// -------------------

__global__ void indexShiftDown(int *d_rows, const int m){
  unsigned int xIndex = blockDim.x * blockIdx.x + threadIdx.x;
  if (xIndex < m) d_rows[xIndex] = d_rows[xIndex]-1;
}

__global__ void indexShiftUp(int *d_rows, const int m){
  unsigned int xIndex = blockDim.x * blockIdx.x + threadIdx.x;
  if (xIndex < m) d_rows[xIndex] = d_rows[xIndex]+1;
}

__global__ void fill(int *vec, int value, int n){
  unsigned int xIndex = blockDim.x * blockIdx.x + threadIdx.x;
  if (xIndex < n) vec[xIndex]=value;
}

// -------------------
// Boundary matrix and PH set-up & preprocessing
// -------------------

__global__ void create_rows_mp(int *d_rows, int *d_cols, int *d_rows_mp, int m, int p, int nnz){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    // We assume col major order
    if(tid < nnz){ // tid == col
        if (tid == 0 || d_cols[tid] != d_cols[tid-1]){
            int i = 0;
            int col = d_cols[tid];
            while ((tid+i < nnz) && col == d_cols[tid+i]){
                d_rows_mp[col*p+i] = d_rows[tid+i];
                i++;
            }
        }
    }
}

 __global__ void compute_low_mp(int *d_rows_mp, int *d_low, int m, int p){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    // We assume col major order
    if(tid < m){
        int low = -1;
        int idx = tid*p;
        int idx_MAX = (tid+1)*p;
        while (idx < idx_MAX && d_rows_mp[idx] != -1){
            low = d_rows_mp[idx++];
        }
        d_low[tid] = low;
    }
}


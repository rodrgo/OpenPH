
__device__ int d_lock = 0;
__global__ void compute_dims_order(int *d_dims, int *d_dims_order, int *d_last_pos, int m, int *d_sentinel){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        int j = tid;
        // set lock
        //printf("{tid=%d, lock=%d}, ", tid, lock);
        do {} while(d_lock != j);
        // do stuff
        int dim_j = d_dims[tid];
        d_dims_order[tid] = d_last_pos[dim_j+1]+1;
        d_last_pos[dim_j+1] += 1;
        //printf("[j=%d, dim_j=%d, lock=%d, d_last_pos=%d], ", j, dim_j, lock, d_last_pos[dim_j+1]);
        // free lock
        d_lock = j+1;
        __syncthreads();
    }
}

__global__ void get_dims_order_start(int *d_dims, int *d_dims_order, int *d_dims_order_start, int m){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        if (d_dims_order[tid] == 0){
            int cdim_pos = d_dims[tid] + 1;
            d_dims_order_start[cdim_pos] = tid;
        }
    }
}

__global__ void count_simplices_dim(int *d_dim_count, int *d_dims, int m){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        if (d_dims[tid] > -1){
            atomicAdd(d_dim_count+d_dims[tid], 1);
        }
    }
}

__global__ void get_simplex_dimensions(int *d_simplex_dimensions, int *d_rows_mp, int m, int p){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if(tid < m){
        int idx = tid*p;
        int idx_MAX = (tid+1)*p;
        int dim = -1;
        while (idx < idx_MAX && d_rows_mp[idx] != -1){
            dim++;
            idx++;
        }
        d_simplex_dimensions[tid] = dim;
    }
}

inline void compute_simplex_dimensions_h(int *h_cols, int m, int p, int nnz, int *d_dims, int *d_dims_order, int *d_dims_order_next, int *p_complex_dimension){
    // h_rows, h_cols are MATLAB vectors, so are supported on [m] 
    // This one we compute on the host for the moment
    int *h_dims = (int*)malloc( sizeof(int) * m );
    int *h_dims_order = (int*)malloc( sizeof(int) * m );
    // Get simplex dimensions
    for (int i = 0; i < m; i++)
        h_dims[i] = -1;
    for (int i = 0; i < nnz; i++)
        h_dims[h_cols[i]-1] += 1;

    int complex_dim = -1;
    for (int i = 0; i < m; i++)
        complex_dim = h_dims[i] > complex_dim ? h_dims[i] : complex_dim;
    *p_complex_dimension = complex_dim;

    // Dimensions are {-1, 0, 1, ..., complex_dim}
    int cdim = complex_dim + 2;
    int *h_dims_order_next = (int*)malloc( sizeof(int) * m );
    int *h_past_cdim = (int*)malloc( sizeof(int) * cdim );
    int *h_dims_order_aux = (int*)malloc( sizeof(int) * cdim );

    for (int i = 0; i < cdim; i++)
        h_dims_order_aux[i] = 0;
    for (int i = 0; i < m; i++)
        h_dims_order_next[i] = -1;
    for (int i = 0; i < cdim; i++)
        h_past_cdim[i] = -1;
    int cdim_pos;
    for (int i = 0; i < m; i++){
        cdim_pos = h_dims[i]+1;
        h_dims_order[i] = h_dims_order_aux[cdim_pos];
        h_dims_order_aux[cdim_pos] += 1;
        if (h_past_cdim[cdim_pos] > -1){
            h_dims_order_next[h_past_cdim[cdim_pos]] = i;
        }
        h_past_cdim[cdim_pos] = i;
    }
    // Copy to device
    cudaMemcpy(d_dims, h_dims, m*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dims_order, h_dims_order, m*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_dims_order_next, h_dims_order_next, m*sizeof(int), cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();

    // free
    free(h_past_cdim);
    free(h_dims_order_aux);
    free(h_dims_order_next);
    free(h_dims);
    free(h_dims_order);
}


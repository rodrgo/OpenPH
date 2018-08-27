
__device__ int d_pos  = 0;
__device__ int d_lock_twist  = 0;
__global__ void find_simplices_with_dim(int *d_simplex_dimensions, int *d_simplices_dim, int dim, int m){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        // set lock
        do {} while(atomicCAS(&d_lock_twist, tid, -1) == tid);
        // do stuff
        if(d_simplex_dimensions[tid] == dim){
            int x = atomicAdd(&d_pos, 1); 
            d_simplices_dim[x] = tid ;
        }
        // free lock
        d_lock_twist = tid+1;
    }
}

__global__ void twist_step(int j, int *d_rows_mp, int *d_low, int m, int p){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if(j == tid && j < m){
        int low_j = d_low[j];
        if (low_j > -1){
            clear_column(low_j, d_rows_mp, p);
            d_low[low_j] = -1;
        }
    }
}

inline void twist(int *d_rows_mp, int *d_aux_mp, int *d_low, int *d_arglow, int *d_dims, int complex_dim, const int m, const int p, float *resRecord, float *timeRecord, int *p_iter, dim3 NBm, dim3 TPBm){

    // Copy to host
    int *h_simplices_with_dim;
    h_simplices_with_dim = (int*)malloc( m * sizeof(int) );
    cudaMemcpy(h_simplices_with_dim, d_dims, m*sizeof(int), cudaMemcpyDeviceToHost);

    int *d_simplices_with_dim;
    cudaMalloc((void**)&d_simplices_with_dim, m * sizeof(int));

    for(int dim = complex_dim; dim > 0; dim--){

        // Bring simplices of dimension d to host 
        int zero = 0;
        cudaMemcpyToSymbol(d_lock_twist, &zero, sizeof(int));
        cudaMemcpyToSymbol(d_pos,  &zero, sizeof(int));
        fill<<< NBm, TPBm >>>(d_simplices_with_dim, -1, m);
        find_simplices_with_dim<<<NBm, TPBm>>>(d_dims, d_simplices_with_dim, dim, m);

        // Bring d_simplices_dim back to host
        cudaMemcpy(h_simplices_with_dim, d_simplices_with_dim, m*sizeof(int), cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();

        int idx = 0;
        while (idx < m && h_simplices_with_dim[idx] > -1){
            int j = h_simplices_with_dim[idx]; 
            reduce_col<<<NBm, TPBm>>>(j, d_rows_mp, d_aux_mp, d_low, d_arglow, m, p);
            twist_step<<<NBm, TPBm>>>(j, d_rows_mp, d_low, m, p);
            cudaDeviceSynchronize();
            idx++;
        }
    } 

    free(h_simplices_with_dim);
    cudaFree(d_simplices_with_dim);

}


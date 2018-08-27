
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

inline void twist(int *d_rows_mp, int *d_aux_mp, int *d_low, int *d_arglow, int *d_dims, int *d_dims_order, int *d_dims_order_next, int *d_dims_order_start, int complex_dim, const int m, const int p, float *resRecord, float *timeRecord, int *p_iter, dim3 NBm, dim3 TPBm){

    int cdim = complex_dim + 2;

    int *h_dims_order_start;
    h_dims_order_start = (int*)malloc( cdim * sizeof(int) );
    cudaMemcpy(h_dims_order_start, d_dims_order_start, cdim*sizeof(int), cudaMemcpyDeviceToHost);

    int *h_dims_order_next;
    h_dims_order_next = (int*)malloc( m * sizeof(int) );
    cudaMemcpy(h_dims_order_next, d_dims_order_next, m*sizeof(int), cudaMemcpyDeviceToHost);

    for(int dim = complex_dim; dim > 0; dim--){

        int dim_pos = dim+1;
        int j = h_dims_order_start[dim_pos];
        int iterate = 1;

        while (iterate && j > -1){
            reduce_col<<<NBm, TPBm>>>(j, d_rows_mp, d_aux_mp, d_low, d_arglow, m, p);
            cudaDeviceSynchronize();

            twist_step<<<NBm, TPBm>>>(j, d_rows_mp, d_low, m, p);
            cudaDeviceSynchronize();

            // Iterator
            if (h_dims_order_next[j] == -1){
                iterate = 0;
            }else{
                j = h_dims_order_next[j];
            }
        }
    } 

    free(h_dims_order_start);
    free(h_dims_order_next);

}


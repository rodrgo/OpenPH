
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

inline void twist_parallel(int *d_low, int *d_arglow, int *d_classes,
        int *d_ess, int *d_rows_mp, const int m, const int p,
        int *d_dims, int *d_dims_order, int *d_dims_order_next,
        int *d_dims_order_start, int complex_dim,
        int *d_aux_mp, int *d_low_true, int *d_ess_true, float *d_float_m,
        float *error_lone, float *error_linf, float *error_redu,
        float *error_ess, float *time_track, int *p_iter, dim3 NBm,
        dim3 TPBm){

    // time
    float time = 0.0;

    // iter and trackers
    track(0, m, d_low, d_ess, d_classes,
            d_low_true, d_ess_true, d_float_m,
            error_lone, error_linf, error_redu,
            error_ess, time_track, time, NBm, TPBm);

    // start algo

    int cdim = complex_dim + 2;

    int *h_dims_order_start;
    h_dims_order_start = (int*)malloc( cdim * sizeof(int) );
    cudaMemcpy(h_dims_order_start, d_dims_order_start, cdim*sizeof(int), cudaMemcpyDeviceToHost);

    int *h_dims_order_next;
    h_dims_order_next = (int*)malloc( m * sizeof(int) );
    cudaMemcpy(h_dims_order_next, d_dims_order_next, m*sizeof(int), cudaMemcpyDeviceToHost);

    int iter = 1;
    for(int dim = complex_dim; dim > 0; dim--){

        int dim_pos = dim+1;
        int j = h_dims_order_start[dim_pos];
        int iterate = 1;

        while (iterate && j > -1){

            // TIC
            cudaEvent_t start, stop;
            tic(&start, &stop);

            // Work on column "j"
            reduce_col<<<NBm, TPBm>>>(j, d_rows_mp, d_aux_mp, d_low, d_arglow, m, p);
            cudaDeviceSynchronize();

            twist_step<<<NBm, TPBm>>>(j, d_rows_mp, d_low, m, p);
            cudaDeviceSynchronize();

            // update classes (Not necessary for algo to work)
            update_classes<<<NBm, TPBm>>>(d_classes, d_low, d_arglow, m);
            cudaDeviceSynchronize();

            // Essential estimation
            ess_hat<<<NBm, TPBm>>>(d_ess, d_low, d_arglow, m);
            cudaDeviceSynchronize();

            // TOC
            toc(start, stop, &time);

            // meausre progress
            track(iter, m, d_low, d_ess, d_classes,
                    d_low_true, d_ess_true, d_float_m,
                    error_lone, error_linf, error_redu,
                    error_ess, time_track, time, NBm, TPBm);

            // iter
            iter++;

            // Iterator
            if (h_dims_order_next[j] == -1){
                iterate = 0;
            }else{
                j = h_dims_order_next[j];
            }
        }
    } 
    p_iter[0] = iter;

    free(h_dims_order_start);
    free(h_dims_order_next);

}


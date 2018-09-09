
__device__ int d_pivot;
__global__ void mark_neighbours(int i, int *d_is_neighbour, int *d_low, int m){
    int j = threadIdx.x + blockDim.x*blockIdx.x;
    if (j < m){
        if (d_low[j] == i){
            d_is_neighbour[j] = 1;
            atomicMin(&d_pivot, j); 
        }
    }
}

__global__ void left_to_right_neighbours(int *d_is_neighbour, int *d_rows_mp, int *d_aux_mp, int *d_low, int *d_arglow, int m, int p){
    int j = threadIdx.x + blockDim.x*blockIdx.x;
    if (j < m && j > d_pivot && d_is_neighbour[j] == 1 ){
        left_to_right(d_pivot, j, d_rows_mp, d_aux_mp, d_low, m, p);
    }
}

inline void ph_row(int *d_low, int *d_arglow, int *d_classes,
        int *d_ess, int *d_rows_mp, const int m, const int p,
        int *d_aux_mp, int *d_low_true, int *d_ess_true, float * d_float_m,
        float *error_lone, float *error_linf, float *error_redu,
        float *error_ess, float *time_track, int *p_iter, dim3 NBm, dim3 TPBm){

    // time
    float time = 0.0;

    // iter and trackers
    track(0, m, d_low, d_ess, d_classes,
            d_low_true, d_ess_true, d_float_m,
            error_lone, error_linf, error_redu,
            error_ess, time_track, time, NBm, TPBm);

    // d_is_neighbour
    int *d_is_neighbour;
    cudaMalloc((void**)&d_is_neighbour, m * sizeof(int));
    fill<<<NBm, TPBm>>>(d_is_neighbour, 0, m);
    cudaDeviceSynchronize();

    int m_value = m;
    cudaMemcpyToSymbol(d_pivot, &m_value, sizeof(int));

    int iter = 1;
    for(int i = m-1; i > -1; i--){

        // TIC
        cudaEvent_t start, stop;
        tic(&start, &stop);

        // Mark neighbours
        mark_neighbours<<<NBm, TPBm>>>(i, d_is_neighbour, d_low, m);
        cudaDeviceSynchronize();

        // Reduce neighbours
        left_to_right_neighbours<<<NBm, TPBm>>>(d_is_neighbour, d_rows_mp, d_aux_mp, d_low, d_arglow, m, p);
        cudaDeviceSynchronize();

        // Reset
        cudaMemcpyToSymbol(d_pivot, &m_value, sizeof(int));
        fill<<<NBm, TPBm>>>(d_is_neighbour, 0, m);
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

    } 
    p_iter[0] = iter;

    cudaFree(d_is_neighbour);

}


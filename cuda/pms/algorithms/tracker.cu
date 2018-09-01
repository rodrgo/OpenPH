
__device__ int d_int_tracker;

__global__ void compute_norm_linf(int *d_vec_1, int *d_vec_2, int m){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        int diff = abs(d_vec_1[tid] - d_vec_2[tid]);
        atomicMax(&d_int_tracker, diff);
    }
} 

inline int norm_linf(int *d_lows, int *d_lows_baseline, int m, dim3 NBm, dim3 TPBm){
    int zero = 0;
    int ninf;
    cudaMemcpyToSymbol(d_int_tracker, &zero, sizeof(int));
    compute_norm_linf<<<NBm, TPBm>>>(d_lows, d_lows_baseline, m);
    cudaDeviceSynchronize();
    cudaMemcpyFromSymbol(&ninf, d_int_tracker, sizeof(int));
    return ninf;
}

/*
   In our design d_classes[j] != 0 iff "we know j is reduced"
*/
__global__ void count_values_equal_to(int value, int *d_vec, int m){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        if (d_vec[tid] == value){
            atomicAdd(&d_int_tracker, 1);
        }
    }
} 

inline int get_percentage_unreduced(int *d_classes, int m, dim3 NBm, dim3 TPBm){
    int zero = 0;
    int progress;
    cudaMemcpyToSymbol(d_int_tracker, &zero, sizeof(int));
    count_values_equal_to<<<NBm, TPBm>>>(0, d_classes, m);
    cudaDeviceSynchronize();
    cudaMemcpyFromSymbol(&progress, d_int_tracker, sizeof(int));
    return progress;
}

inline void track(int iter, int m, 
        int *d_low, int *d_classes, int *d_low_true, 
        int *d_ess_true, float *error_lone, 
        float *error_linf, float *error_redu, 
        float *error_ess, float *time_track, float time,
        dim3 NBm, dim3 TPBm){

    // Track time
    time_track[iter] = time;

    // Track inf-norm of low estimate against true vector.
    int ninf = norm_linf(d_low, d_low_true, m, NBm, TPBm);
    error_linf[iter] = ((float) ninf);

    // Track percentage of unreduced columns
    int progress = get_percentage_unreduced(d_classes, m, NBm, TPBm);
    error_redu[iter] = ((float) progress);
}



inline float asum_host(int m, float *h_float_m){
    float res = 0.0f;
    for (int i=0; i<m; i++)
        res+=abs(h_float_m[i]);
    return res;
}

__global__ void diff(float *d_res, int *d_vec_1, int *d_vec_2, int m){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        d_res[tid] = ((float) d_vec_1[tid] - d_vec_2[tid]);
    }
} 

inline void diff_host(float *h_res, int *h_vec_1, int *h_vec_2, int m){
    for (int i=0; i<m; i++){
        h_res[i] = ((float) h_vec_1[i] - h_vec_2[i]);
    }
}

__global__ void eq_value(float *d_res, int *d_v, int val, int m){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        if (d_v[tid] == val){
            d_res[tid] = 1.0f;
        }else{
            d_res[tid] = 0.0f;
        }
    }
} 

inline void eq_value_host(float *h_res, int *h_v, int val, int m){
    for (int i=0; i<m; i++){
        if (h_v[i] == val){
            h_res[i] = 1.0f;
        }else{
            h_res[i] = 0.0f;
        }
    }
}

__global__ void eq_vectors(float *d_res, int *d_v1, int *d_v2, int m){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        if (d_v1[tid] == d_v2[tid]){
            d_res[tid] = 1;
        }else{
            d_res[tid] = 0;
        }
    }
} 

__global__ void to_float(float *d_float_m, int *d_v, int m){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        float val;
        val = (float) d_v[tid];
        d_float_m[tid] = val;
    }
} 

inline void to_float_host(float *h_float_m, int *h_v, int m){
    for (int i=0; i<m; i++){
        h_float_m[i] = ((float) h_v[i]);
    }
}

inline float norm_1(float *d_float_m, int m){
    return cublasSasum(m, d_float_m, 1);
}

inline float norm_1_host(float *h_float_m, int m){
    return asum_host(m, h_float_m);
}

inline float norm_inf(float *d_float_m, int m){
    float ninf;
    int pos = cublasIsamax(m, d_float_m, 1)-1;
    cudaMemcpy(&ninf, d_float_m+pos, sizeof(float), cudaMemcpyDeviceToHost);
    return abs(ninf); 
}

inline float norm_inf_host(float *h_float_m, int m){
    float ninf;
    for (int i=0; i<m; i++)
        if (abs(h_float_m[i]) > ninf)
            ninf = abs(h_float_m[i]);
    return abs(ninf);
}

inline void track(int iter, int m, 
        int *d_low, int *d_ess, int *d_classes, int *d_low_true, 
        int *d_ess_true, float *d_float_m, float *error_lone, 
        float *error_linf, float *error_redu, 
        float *error_ess, float *time_track, float time,
        dim3 NBm, dim3 TPBm){

    // Time
    time_track[iter] = time;

    // (float) d_low_err = d_low - d_low_true
    diff<<<NBm, TPBm>>>(d_float_m, d_low, d_low_true, m);
    cudaDeviceSynchronize();

    // norm_lone(d_low, d_low_true)
    error_lone[iter] = norm_1(d_float_m, m);

    // norm_linf(d_low, d_low_true)
    error_linf[iter] = norm_inf(d_float_m, m);

    // |j : d_classes[j] = 0| (number of unreduced columns)
    eq_value<<<NBm, TPBm>>>(d_float_m, d_classes, 0, m);
    cudaDeviceSynchronize();
    error_redu[iter] = cublasSasum(m, d_float_m, 1) / ((float) m);
    cudaDeviceSynchronize();

    // |j : d_ess[j] = d_ess_true[j]|/sum(d_ess_true)
    to_float<<<NBm, TPBm>>>(d_float_m, d_ess_true, m);
    cudaDeviceSynchronize();
    float num_ess_true = cublasSasum(m, d_float_m, 1);

    to_float<<<NBm, TPBm>>>(d_float_m, d_ess, m);
    cudaDeviceSynchronize();
    float num_ess_hat = cublasSasum(m, d_float_m, 1);

    if (num_ess_hat > 0){
        error_ess[iter] = num_ess_true/num_ess_hat;
    }else{
        error_ess[iter] = -1.0;
    }

}

inline void track_host(int iter, int m, 
        int *h_low, int *h_ess, int *h_classes, int *h_low_true, 
        int *h_ess_true, float *h_float_m, float *error_lone, 
        float *error_linf, float *error_redu, 
        float *error_ess, float *time_track, float time){

    // Time
    time_track[iter] = time;

    // (float) d_low_err = d_low - d_low_true
    diff_host(h_float_m, h_low, h_low_true, m);

    // norm_lone(d_low, d_low_true)
    error_lone[iter] = norm_1_host(h_float_m, m);

    // norm_linf(d_low, d_low_true)
    error_linf[iter] = norm_inf_host(h_float_m, m);

    // |j : d_classes[j] = 0| (number of unreduced columns)
    eq_value_host(h_float_m, h_classes, 0, m);
    error_redu[iter] = asum_host(m, h_float_m) / ((float) m);

    // |j : d_ess[j] = d_ess_true[j]|/sum(d_ess_true)
    to_float_host(h_float_m, h_ess_true, m);
    float num_ess_true = asum_host(m, h_float_m);

    to_float_host(h_float_m, h_ess, m);
    float num_ess_hat = asum_host(m, h_float_m);

    if (num_ess_hat > 0){
        error_ess[iter] = num_ess_true/num_ess_hat;
    }else{
        error_ess[iter] = -1.0;
    }

}


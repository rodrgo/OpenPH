
inline void h2d(int *d_v, int j, int v){
    cudaMemcpy(d_v+j, &v, sizeof(int), cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
    return;
}

inline int d2h(int *d_v, int j){
    int v_j = 0;
    cudaMemcpy(&v_j, d_v+j, sizeof(int), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    return v_j;
}

inline int arglow(int *d_arglow, int low_j){
    int v = -1;
    if (low_j > -1){
        cudaMemcpy(&v, d_arglow+low_j, sizeof(int), cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
    }
    return v;
}

inline void reduce_col(int j, int *d_rows_mp, int *d_aux_mp, int *d_low, int *d_arglow, int m, int p, dim3 numBlocks_m, dim3 threadsPerBlock_m){
    // j = 0, 1, ..., m-1
    int j0 = -1;
    int low_j = d2h(d_low, j); // low_j = -1, 0, 1, ..., m-1
    int low_j_aux = 0;
    while (low_j > -1 && d2h(d_arglow, low_j) != -1){
        j0 = d2h(d_arglow, low_j);
        // left_to_right also updates low
        left_to_right<<<numBlocks_m, threadsPerBlock_m>>>(j0, j, d_rows_mp, d_aux_mp, d_low, m, p);
        cudaDeviceSynchronize();
        low_j_aux = low_j;
        low_j = d2h(d_low, j);
    }
    low_j = d2h(d_low, j);
    if (low_j > -1){
        h2d(d_arglow, low_j, j);
    } 
}


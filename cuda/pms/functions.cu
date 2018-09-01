
inline int max_nnz(int *h_index, int m, int nnz){
    // Get maximum nnz in rows/columns from row/column index.
    // h_index_dim is either h_rows or h_cols

    int *h_nnz_index = (int*)malloc( sizeof(int) * m );

    for(int l=0; l<m; l++)   h_nnz_index[l] = 0;
    for(int l=0; l<nnz; l++) h_nnz_index[h_index[l]-1] += 1;

    int max_nnz_index = 0;

    for(int l=0; l<m; l++){
        if (h_nnz_index[l] > max_nnz_index)
            max_nnz_index = h_nnz_index[l];
    }

    free(h_nnz_index);

    return max_nnz_index;
} 

inline int num_blocks(int m, int threads_perblock){
    return (int)ceil((float)m/(float)threads_perblock);
}

inline int assert_col_order(int *h_cols, int nnz){
    int is_col_order = 1;
    for (int l=1; l<nnz; l++)
        if (h_cols[l] < h_cols[l-1]){
            is_col_order = 0; 
            break;
        }
    return is_col_order;
}

inline void set_gpu_device(int gpuNumber, int *p_threads_perblock_m, int *p_threads_perblock_nnz, int *p_threads_perblock_mp, int m, int nnz, int p){

    // gpuNumber
    unsigned int max_threads_per_block;

    cudaDeviceProp dp;
    cudaSetDevice(gpuNumber);
    cudaGetDeviceProperties(&dp,gpuNumber);
    max_threads_per_block = dp.maxThreadsPerBlock;

    int devCount;
    cudaGetDeviceCount(&devCount);
    if ((gpuNumber >= devCount) && (gpuNumber != 0)){
        cout << "This computer has " << devCount 
        << " gpus and gpuNumber was" << endl << "selected at "
        << gpuNumber << " which is larger than admissible." 
        << endl << "gpuNumber has been reset to 0." << endl; 
        gpuNumber = 0;
    }
    cudaSetDevice(gpuNumber);
    cudaGetDeviceProperties(&dp,gpuNumber);
    max_threads_per_block = dp.maxThreadsPerBlock;

    int mp = m * p;

    *p_threads_perblock_m   = min(m, max_threads_per_block); 
    *p_threads_perblock_nnz = min(nnz, max_threads_per_block);
    *p_threads_perblock_mp  = min(mp, max_threads_per_block);

}

inline void create_rows(int *h_rows, int *h_cols, int* d_rows_mp, int m, int p, int nnz, dim3 NBnnz, dim3 TPBnnz, dim3 NBmp, dim3 TPBmp){

    int *d_rows_in, *d_cols_in;

    cudaMalloc((void**)&d_rows_in, nnz * sizeof(int));
    cudaMalloc((void**)&d_cols_in, nnz * sizeof(int));

    cudaMemcpy(d_rows_in, h_rows, sizeof(int)*nnz, cudaMemcpyHostToDevice);
    cudaMemcpy(d_cols_in, h_cols, sizeof(int)*nnz, cudaMemcpyHostToDevice);

    indexShiftDown<<<NBnnz, TPBnnz>>>(d_rows_in, nnz);
    indexShiftDown<<<NBnnz, TPBnnz>>>(d_cols_in, nnz); 
    cudaDeviceSynchronize();

    int mp = m * p;

    fill<<<NBmp, TPBmp>>>(d_rows_mp, -1, mp);
    create_rows_mp<<<NBnnz, TPBnnz>>>(d_rows_in, d_cols_in, d_rows_mp, m, p, nnz);
    cudaDeviceSynchronize();

    cudaFree(d_cols_in);
    cudaFree(d_rows_in);

}

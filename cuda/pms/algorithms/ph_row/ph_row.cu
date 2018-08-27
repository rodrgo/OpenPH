
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

inline void ph_row(int *d_rows_mp, int *d_aux_mp, int *d_low, int *d_arglow, const int m, const int p, float *resRecord, float *timeRecord, int *p_iter, dim3 NBm, dim3 TPBm){

    // d_is_neighbour
    int *d_is_neighbour;
    cudaMalloc((void**)&d_is_neighbour, m * sizeof(int));
    fill<<<NBm, TPBm>>>(d_is_neighbour, 0, m);
    cudaDeviceSynchronize();

    int m_value = m;
    cudaMemcpyToSymbol(d_pivot, &m_value, sizeof(int));

    for(int i = m-1; i > -1; i--){

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

    } 

    cudaFree(d_is_neighbour);

}


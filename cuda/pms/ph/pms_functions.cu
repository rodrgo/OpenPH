
__device__ int is_reduced = 1;
void __global__ matrix_is_reduced(int *d_lows, int *d_aux, int m){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        int low_j = d_lows[tid];
        if (low_j > -1){
            atomicAdd(d_aux+low_j, 1);
            if (d_aux[low_j] > 1)
                is_reduced = 0;
        }
    }
} 

void __global__ alpha_beta_reduce(int *d_lows, int *d_beta, int *d_classes, int *d_rows_mp, int *d_arglow, int *d_lowstar, int m){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        int alpha = d_lows[tid];
        int beta  = d_beta[tid];
        if (alpha == beta && beta > -1){
            // tid is "negative"
            d_classes[tid] = -1;
            int pos_pair = d_beta[tid];
            clear_column(pos_pair, d_rows_mp, p);
            d_arglow[pos_pair] = beta;
            d_lowstar[pos_pair] = -1;
            d_lowstar[beta] = pos_pair;
        }
    }
}

void __global__ count_simplices_dim(int *d_dim_count, int *d_simplex_dimensions, int complex_dimension){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        if (d_simplex_dimensions[tid] > -1){
            atomicAdd(d_dim_count[d_simplex_dimensions[tid]], 1);
        }
    }
}

void __global__ get_dim_pos(int *d_simplex_dims, int *d_dim_count, int *d_dim_pos, int complex_dim){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        if (d_simplex_dimensions[tid] > -1){
            atomicAdd(d_dim_count[d_simplex_dimensions[tid]], 1);
        }
    }
}

inline void compute_simplex_dimensions(int *d_simplex_dimensions, int *d_rows_mp, int m, int p, int *p_complex_dimension, dim3 numBlocks_m, dim3 threadsPerBlock_m){
    get_simplex_dimensions<<<numBlocks_m, threadsPerBlock_m>>>(d_simplex_dimensions, d_rows_mp, m, p);
    thrust::device_ptr<int> dev_ptr = thrust::device_pointer_cast(d_simplex_dimensions);
    thrust::device_ptr<int> max_ptr = thrust::max_element(dev_ptr, dev_ptr + m);
    *p_complex_dimension = max_ptr[0];
}

int is_reduced(int *d_aux, int *d_lows, int m, dim3 numBlocks_m, dim3 threadsPerBlock_m){
    int one = 1;
    int is_reduced;
    cudaMemcpyToSymbol(d_is_reduced, &one, sizeof(int));
    zero_vector_int<<<numBlocks_m, threadsPerBlock_m>>>(d_aux, m);
    matrix_is_reduced<<<numBlocks_m, threadsPerBlock_m>>>(d_lows, d_aux, m);
    cudaMemcpyFromSymbol(&is_reduced, d_is_reduced, sizeof(int));
    return is_reduced;
}

inline void create_beta(int *d_beta, int *h_rows, int *h_cols, int m, int nnz){
    int *h_beta;
    h_beta = (int*)malloc( sizeof(int) * m );
    create_beta_h(h_beta, h_rows, h_cols, m, nnz);
    cudaMemcpy(d_beta, h_beta, m*sizeof(int), cudaMemcpyHostToDevice);
    free(h_beta);
}

in


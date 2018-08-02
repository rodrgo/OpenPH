
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

void __global__ compute_dims_order(int *d_dims, int *d_dims_order, int *d_last_pos, int m){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        int dim_j = d_dims[tid];
        if (tid == 0){
            int pos = d_last_pos[dim_j]+1;
            d_dims_order[tid] = pos;
            d_last_pos[dim_j] = pos;
        }else{
            do {} while (d_dims_order[tid-1] == -1);
            int pos = d_last_pos[dim_j]+1;
            d_dims_order[tid] = pos;
            d_last_pos[dim_j] = pos;
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

void __global__ count_simplices_dim(int *d_dim_count, int *d_dims){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        if (d_dims[tid] > -1){
            atomicAdd(d_dim_count[d_dims[tid]], 1);
        }
    }
}

void __global__ phase_i(int *d_ceilings, int *d_dims, int *d_dims_order, int *d_low, int *d_arglow, int *d_classes, int *d_clear, int *d_visited, int *d_ceil_cdim, int *d_locks_cdim, int m){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        int dim_j = d_dims[tid];
        int j_ord = d_dims_order[tid]; // 0, 1, ...
        int dim_ceil = d_ceil_cdim[dim_j];
        // set lock
        do {} while(atomicCAS(d_locks_cdim[dim_j], j_ord, -1) == j_ord);
        // do stuff
        if (low_j > -1){
            if (d_visited[low_j] == 0){
                d_arglow[low_j] = tid;
                d_classes[tid] = -1;
                d_clear[low_j] = 1;
                d_ceilings[tid] = dim_ceil;
            }else{
                d_ceilings[tid] = dim_ceil;
                d_ceil_cdim[dim_j] = low_j > dim_ceil ? low_j : dim_ceil;
            }
        }
        // free lock
        //d_lock = tid+1;
        d_locks_cdim[dim_j] = j_ord + 1;
    }
}

void __global__ phase_ii(int *d_low, int *d_arglow, int *d_rows_mp, int *d_aux_mp, int m, int p){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        int low_j = d_low[tid];
        if (d_arglow[low_j] > -1){
            int pivot = d_arglow[low_j];
            if (pivot < j){
                left_to_right_device(pivot, j, d_rows_mp, d_aux_mp, d_low, m, p);
                d_updated[tid] = 1;
                // alpha_beta_check 
                low_j = d_low[tid];
                if (low_j > -1){
                    if (d_beta[low_j] == tid){
                        // is lowstar, do a twist clearing
                        d_arglow[low_j] = tid;
                        d_classes[tid] = -1;
                        clear_column(low_j, d_rows_mp, p);
                    }
                }else{
                    d_classes[tid] = 1;
                }
            }
        }
    }
}

void __global__ set_unmarked(int *d_classes, int *d_low, int *d_arglow, int *d_rows_mp, int m, int p){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        if (d_classes[tid] == 0)
            if (d_low[tid] > -1){
                d_arglow[tid] = tid;
                d_classes[tid] = -1;
                clear_column(tid, d_rows_mp, p);
            }
        }else{
            d_classes[tid] = 2; 
        }
    }
}

inline void compute_simplex_dimensions(int *d_dims, int *d_dims_order, int *p_complex_dimension, int *d_rows_mp, int m, int p, dim3 numBlocks_m, dim3 threadsPerBlock_m){
    // d_dims
    get_simplex_dimensions<<<numBlocks_m, threadsPerBlock_m>>>(d_dims, d_rows_mp, m, p);
    // d_complex_dim
    thrust::device_ptr<int> dev_ptr = thrust::device_pointer_cast(d_dims);
    thrust::device_ptr<int> max_ptr = thrust::max_element(dev_ptr, dev_ptr + m);
    *p_complex_dimension = max_ptr[0];
}

inline void compute_dimension_order(int *d_dims, int *d_dims_order, int *d_last_pos, int cdim, int m, dim3 numBlocks_m, dim3 threadsPerBlock_m){
    fill<<<numBlocks_m, threadsPerBlock_m>>>(d_last_pos, -1, cdim);
    fill<<<numBlocks_m, threadsPerBlock_m>>>(d_dims_order, -1, m);
    compute_dims_order<<<numBlocks_m, threadsPerBlock_m>>>(d_dims, d_dims_order, d_last_pos, m);
    fill<<<numBlocks_m, threadsPerBlock_m>>>(d_last_pos, -1, cdim);
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


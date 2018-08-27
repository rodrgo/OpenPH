
// -----------------------
// Check if matrix is reduced
// -----------------------

__device__ int d_is_reduced = 1;
void __global__ matrix_is_reduced(int *d_lows, int *d_aux, int m){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        int low_j = d_lows[tid];
        if (low_j > -1){
            atomicAdd(d_aux+low_j, 1);
            if (d_aux[low_j] > 1)
                d_is_reduced = 0;
        }
    }
} 

int is_reduced(int *d_aux, int *d_lows, int m, dim3 numBlocks_m, dim3 threadsPerBlock_m){
    int one = 1;
    int is_reduced;
    cudaMemcpyToSymbol(d_is_reduced, &one, sizeof(int));
    fill<<<numBlocks_m, threadsPerBlock_m>>>(d_aux, 0, m);
    matrix_is_reduced<<<numBlocks_m, threadsPerBlock_m>>>(d_lows, d_aux, m);
    cudaMemcpyFromSymbol(&is_reduced, d_is_reduced, sizeof(int));
    return is_reduced;
}

// -----------------------
// Phase 0
// -----------------------

void __global__ mark_pivots_and_clear(int *d_low, int *d_beta, int *d_classes, int *d_rows_mp, int *d_arglow, int m, int p){
    int j = threadIdx.x + blockDim.x*blockIdx.x;
    if (j < m){
        int low_j = d_low[j]; 
        int beta_j = d_beta[j];
        // Check if is pivot
        if (low_j == beta_j && beta_j > -1){
            // j is "negative"
            d_classes[j] = -1;
            // low_j is positive
            clear_column(low_j, d_rows_mp, p);
            d_low[low_j] = -1;
            d_classes[low_j] = 1;
            // Record j as pivot
            d_arglow[low_j] = j;
        }
    }
}

void __global__ transverse_dimensions(int *d_dims, int *d_dims_order, int *d_dims_order_next, int *d_dims_order_start, int *d_low, int *d_arglow, int *d_classes, int *d_clear, int *d_visited, int *d_ceil_cdim, int cdim){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < cdim){
        int dim_j; // -1, 0, 1, ..., complex_dim
        int cdim_pos;
        int dim_ceil; // initialized at -1 
        int low_j;
        int j = d_dims_order_start[tid];
        int iterate = 1;
        while (iterate && j > -1){
            dim_j = d_dims[j];
            cdim_pos = dim_j + 1; 
            dim_ceil = d_ceil_cdim[cdim_pos];
            low_j = d_low[j];
            if (low_j > -1){
                if (d_visited[low_j] == 0){
                    if (d_classes[j] == 0 && low_j > dim_ceil){
                        d_arglow[low_j] = j;
                        d_classes[j] = -1;
                        d_clear[low_j] = 1;
                    }
                }else{
                    d_ceil_cdim[cdim_pos] = low_j > dim_ceil ? low_j : dim_ceil;
                }
                d_visited[low_j] = 1;
            }
            // Iterator
            if (d_dims_order_next[j] == -1){
                iterate = 0;
            }else{
                j = d_dims_order_next[j];
            }
        }
    }
}

void __global__ phase_ii(int *d_low, int *d_left, int *d_classes, int *d_arglow, int *d_rows_mp, int *d_aux_mp, int m, int p){
    int j = threadIdx.x + blockDim.x*blockIdx.x;
    if (j < m){
        int low_j = d_low[j];
        int pivot = d_arglow[low_j];
        if (-1 < pivot && pivot < j && d_classes[j] == 0){
            left_to_right(pivot, j, d_rows_mp, d_aux_mp, d_low, m, p);
            // alpha_beta_check 
            low_j = d_low[j];
            if (low_j > -1){
                if (d_left[low_j] == j){
                    // is lowstar, do a twist clearing
                    d_arglow[low_j] = j;
                    d_classes[j] = -1;
                    clear_column(low_j, d_rows_mp, p);
                    d_classes[low_j] = 1;
                }
            }else{
                d_classes[j] = 1;
            }
        }
    }
}

void __global__ set_unmarked(int *d_classes, int *d_low, int *d_arglow, int *d_rows_mp, int m, int p){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        if (d_classes[tid] == 0){
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

void __global__ clear_positives(int *d_clear, int *d_low, int *d_classes, int *d_rows_mp, int m, int p){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        if (d_clear[tid] == 1){
            d_low[tid] = -1;
            d_classes[tid] = 1;
            clear_column(tid, d_rows_mp, p);
        }
    }
}


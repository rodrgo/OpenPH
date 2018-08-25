   
__device__ int d_lock = 0;
__device__ int d_pos  = 0;
__global__ void find_simplices_with_dim(int *d_simplex_dimensions, int *d_simplices_dim, int dim, int m){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if (tid < m){
        // set lock
        do {} while(atomicCAS(&d_lock, tid, -1) == tid);
        // do stuff
        if(d_simplex_dimensions[tid] == dim){
            int x = atomicAdd(&d_pos, 1); 
            d_simplices_dim[x] = tid ;
        }
        // free lock
        d_lock = tid+1;
    }
}

 __global__ void get_simplex_dimensions(int *d_simplex_dimensions, int *d_rows_mp, int m, int p){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if(tid < m){
        int idx = tid*p;
        int idx_MAX = (tid+1)*p;
        int dim = -1;
        while (idx < idx_MAX && d_rows_mp[idx] != -1){
            dim++;
            idx++;
        }
        d_simplex_dimensions[tid] = dim;
    }
}

 __global__ void clear_cols_mp(int j, int *d_rows_mp, int m, int p){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if(tid == j){
        int idx = tid*p;
        int idx_MAX = (tid+1)*p;
        while (idx < idx_MAX && d_rows_mp[idx] != -1){
            d_rows_mp[idx++] = -1;
        }
    }
}

 __global__ void compute_low_mp(int *d_rows_mp, int *d_low, int m, int p){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    // We assume col major order
    if(tid < m){
        int low = -1;
        int idx = tid*p;
        int idx_MAX = (tid+1)*p;
        while (idx < idx_MAX && d_rows_mp[idx] != -1){
            low = d_rows_mp[idx++];
        }
        d_low[tid] = low;
    }
}

__global__ void create_rows_mp(int *d_rows, int *d_cols, int *d_rows_mp, int m, int p, int nnz){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    // We assume col major order
    if(tid < nnz){ // tid == col
        if (tid == 0 || d_cols[tid] != d_cols[tid-1]){
            int i = 0;
            int col = d_cols[tid];
            while ((tid+i < nnz) && col == d_cols[tid+i]){
                d_rows_mp[col*p+i] = d_rows[tid+i];
                i++;
            }
        }
    }
}

__device__ void left_to_right_device(int j0, int j1, int *d_rows_mp, int *d_aux_mp, int *d_low, int m, int p){
    // Compute symmetric difference of supp(j0) and supp(j1) and store in d_aux
    // If rows are initially sorted, this returns a sorted list
    int idx0 = j0*p; 
    int idx1 = j1*p; 
    int idx0_MAX = (j0+1)*p; 
    int idx1_MAX = (j1+1)*p; 
    int idx = idx1;
    bool idx0_ok = d_rows_mp[idx0] != -1 && idx0 < idx0_MAX;
    bool idx1_ok = d_rows_mp[idx1] != -1 && idx1 < idx1_MAX;
    while (idx0_ok || idx1_ok){
        if (idx0_ok && idx1_ok){
            if (d_rows_mp[idx0] < d_rows_mp[idx1]){
                d_aux_mp[idx++] = d_rows_mp[idx0++];
            }else if (d_rows_mp[idx1] < d_rows_mp[idx0]){
                d_aux_mp[idx++] = d_rows_mp[idx1++];
            }else{
                idx0++;
                idx1++;
                if (idx0 == idx0_MAX-1 || idx1 == idx1_MAX-1)
                    printf("WARNING: Column reaching memalloc limit\n");
            }
        }else{
            if (idx0_ok){
                d_aux_mp[idx++] = d_rows_mp[idx0++];
            }
            if (idx1_ok){
                d_aux_mp[idx++] = d_rows_mp[idx1++];
            }
        }
        idx0_ok = d_rows_mp[idx0] != -1 && idx0 < idx0_MAX;
        idx1_ok = d_rows_mp[idx1] != -1 && idx1 < idx1_MAX;
    }
    int low_j1 = -1;
    // At least one value was written in d_aux_mp
    for (idx1 = j1*p; idx1 < idx1_MAX; idx1++){
        d_rows_mp[idx1] = d_aux_mp[idx1];
        d_aux_mp[idx1] = -1;
        if (d_rows_mp[idx1] > -1)
            low_j1 = d_rows_mp[idx1];
    }
    d_low[j1] = low_j1;
}

__device__ void clear_column(int j, int *d_rows_mp, int p){
    int idx = j*p; 
    int idx_MAX = (j+1)*p; 
    while (idx < idx_MAX && d_rows_mp[idx] != -1){
        d_rows_mp[idx++] = -1;
    }
}

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

__global__ void reduce_col_gpu(int j, int *d_rows_mp, int *d_aux_mp, int *d_low, int *d_arglow, int m, int p){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if(j == tid && j < m){
        int j0 = -1;
        int low_j = d_low[j]; // low_j = -1, 0, 1, ..., m-1
        while (low_j > -1 && d_arglow[low_j] != -1){
            j0 = d_arglow[low_j];
            left_to_right_device(j0, j, d_rows_mp, d_aux_mp, d_low, m, p);
            low_j = d_low[j];
        }
        low_j = d_low[j];
        if (low_j > -1){
            d_arglow[low_j] = j;
        }
    }
}

 __global__ void left_to_right(int j0, int j1, int *d_rows_mp, int *d_aux_mp, int *d_low, int m, int p){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    // Assume:
    //  j0 < j1
    //  col major order
    //  nonzeros in d_rows are contiguous and in order
    //  d_aux is set to -1
    if(j1 == tid & tid < m){
        // Compute symmetric difference of supp(j0) and supp(j1) and store in d_aux
        // If rows are initially sorted, 
        // this returns a sorted list
        int idx0 = j0*p; 
        int idx1 = j1*p; 
        int idx0_MAX = (j0+1)*p; 
        int idx1_MAX = (j1+1)*p; 
        int idx = idx1;
        bool idx0_ok = d_rows_mp[idx0] != -1 && idx0 < idx0_MAX;
        bool idx1_ok = d_rows_mp[idx1] != -1 && idx1 < idx1_MAX;
        bool DEBUG = 1 == 0;
        if (DEBUG){
            for (int i = 0; i < p; i++){
                printf("(j0, %d, %d), ", i, d_rows_mp[j0*p+i]);
            } 
            printf("\n");
            for (int i = 0; i < p; i++){
                printf("(j1, %d, %d), ", i, d_rows_mp[j1*p+i]);
            } 
            printf("\n");
        }
        while (idx0_ok || idx1_ok){
            if (idx0_ok && idx1_ok){
                if (d_rows_mp[idx0] < d_rows_mp[idx1]){
                    d_aux_mp[idx++] = d_rows_mp[idx0++];
                }else if (d_rows_mp[idx1] < d_rows_mp[idx0]){
                    d_aux_mp[idx++] = d_rows_mp[idx1++];
                }else{
                    idx0++;
                    idx1++;
                }
            }else{
                if (idx0_ok){
                    d_aux_mp[idx++] = d_rows_mp[idx0++];
                }
                if (idx1_ok){
                    d_aux_mp[idx++] = d_rows_mp[idx1++];
                }
            }
            idx0_ok = d_rows_mp[idx0] != -1 && idx0 < idx0_MAX;
            idx1_ok = d_rows_mp[idx1] != -1 && idx1 < idx1_MAX;
        }
        int low_j1 = -1;
        // At least one value was written in d_aux_mp
        for (idx1 = j1*p; idx1 < idx1_MAX; idx1++){
            d_rows_mp[idx1] = d_aux_mp[idx1];
            d_aux_mp[idx1] = -1;
            if (d_rows_mp[idx1] > -1)
                low_j1 = d_rows_mp[idx1];
        }
        d_low[j1] = low_j1;
        // DEBUG
        if (DEBUG){
            for (int i = 0; i < p; i++){
                printf("(j0+j1, %d, %d), ", i, d_rows_mp[j1*p+i]);
            } 
            printf("\n");
        }
    }
}


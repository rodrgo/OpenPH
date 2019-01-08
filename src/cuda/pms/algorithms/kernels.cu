
__device__ void left_to_right(int j0, int j1, int *d_rows_mp, int *d_aux_mp, int *d_low, int m, int p){
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
                if (idx0 == idx0_MAX || idx1 == idx1_MAX){
                    printf("WARNING: Reached memalloc limit\n");
                    asm("trap;");
                }
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
    idx1 = j1*p;
    while(d_aux_mp[idx1] != -1 && idx1 < idx1_MAX){
        d_rows_mp[idx1] = d_aux_mp[idx1];
        d_aux_mp[idx1] = -1;
        if (d_rows_mp[idx1] > -1)
            low_j1 = d_rows_mp[idx1];
        idx1++;
    }
    while(d_rows_mp[idx1] != -1 && idx1 < idx1_MAX){
        d_rows_mp[idx1] = -1;
        idx1++;
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

__global__ void reduce_col(int j, int *d_rows_mp, int *d_aux_mp, int *d_low, int *d_arglow, int m, int p){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if(j == tid && j < m){
        int j0 = -1;
        int low_j = d_low[j]; // low_j = -1, 0, 1, ..., m-1
        while (low_j > -1 && d_arglow[low_j] != -1){
            j0 = d_arglow[low_j];
            left_to_right(j0, j, d_rows_mp, d_aux_mp, d_low, m, p);
            low_j = d_low[j];
        }
        low_j = d_low[j];
        if (low_j > -1){
            d_arglow[low_j] = j;
        }
    }
}

__global__ void update_classes_final(int *d_classes, int *d_low, int *d_ess, int m){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if(tid < m){
        if (d_low[tid] > -1){
            d_classes[tid] = -1;
            d_classes[d_low[tid]] = 1;
        }else{
            if (d_ess[tid] == 1){
                d_classes[tid] = 2;
            }
        }
    }
}

__global__ void update_classes(int *d_classes, int *d_low, int *d_arglow, int m){
    int tid = threadIdx.x + blockDim.x*blockIdx.x;
    if(tid < m){
        if (d_arglow[tid] > -1){
            d_classes[d_arglow[tid]] = -1;
            d_classes[tid] = 1;
        }
    }
}

__global__ void ess_hat_final(int *d_essential_hat, int *d_low, int *d_arglow, int m){
    int j = threadIdx.x + blockDim.x*blockIdx.x;
    if(j < m){
        if (d_low[j] > -1){
            d_essential_hat[j] = 0;
            d_essential_hat[d_low[j]] = 0;
        }
    }
}

__global__ void ess_hat(int *d_essential_hat, int *d_low, int *d_arglow, int m){
    int j = threadIdx.x + blockDim.x*blockIdx.x;
    if(j < m){
        if (d_low[j] > -1){
            d_essential_hat[d_low[j]] = 0;
        }
        if (d_arglow[j] > -1){
            d_essential_hat[d_arglow[j]] = 0;
            d_essential_hat[j] = 0;
        }
    }
}

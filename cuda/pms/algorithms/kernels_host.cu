
inline void left_to_right_host(int j0, int j1, int *h_rows_mp, int *h_aux_mp, int *h_low, int m, int p){
    // Compute symmetric difference of supp(j0) and supp(j1) and store in d_aux
    // If rows are initially sorted, this returns a sorted list
    int idx0 = j0*p; 
    int idx1 = j1*p; 
    int idx0_MAX = (j0+1)*p; 
    int idx1_MAX = (j1+1)*p; 
    int idx = idx1;
    bool idx0_ok = h_rows_mp[idx0] != -1 && idx0 < idx0_MAX;
    bool idx1_ok = h_rows_mp[idx1] != -1 && idx1 < idx1_MAX;
    while (idx0_ok || idx1_ok){
        if (idx0_ok && idx1_ok){
            if (h_rows_mp[idx0] < h_rows_mp[idx1]){
                h_aux_mp[idx++] = h_rows_mp[idx0++];
            }else if (h_rows_mp[idx1] < h_rows_mp[idx0]){
                h_aux_mp[idx++] = h_rows_mp[idx1++];
            }else{
                idx0++;
                idx1++;
                if (idx0 == idx0_MAX-1 || idx1 == idx1_MAX-1)
                    printf("WARNING: Column reaching memalloc limit\n");
            }
        }else{
            if (idx0_ok){
                h_aux_mp[idx++] = h_rows_mp[idx0++];
            }
            if (idx1_ok){
                h_aux_mp[idx++] = h_rows_mp[idx1++];
            }
        }
        idx0_ok = h_rows_mp[idx0] != -1 && idx0 < idx0_MAX;
        idx1_ok = h_rows_mp[idx1] != -1 && idx1 < idx1_MAX;
    }
    int low_j1 = -1;
    // At least one value was written in d_aux_mp
    for (idx1 = j1*p; idx1 < idx1_MAX; idx1++){
        h_rows_mp[idx1] = h_aux_mp[idx1];
        h_aux_mp[idx1] = -1;
        if (h_rows_mp[idx1] > -1)
            low_j1 = h_rows_mp[idx1];
    }
    h_low[j1] = low_j1;
}

inline void clear_column_host(int j, int *h_rows_mp, int p){
    int idx = j*p; 
    int idx_MAX = (j+1)*p; 
    while (idx < idx_MAX && h_rows_mp[idx] != -1){
        h_rows_mp[idx++] = -1;
    }
}

inline void reduce_col_host(int j, int *h_rows_mp, int *h_aux_mp, int *h_low, int *h_arglow, int m, int p){
    for (int tid=0; tid<m; tid++){
        int j0 = -1;
        int low_j = h_low[j]; // low_j = -1, 0, 1, ..., m-1
        while (low_j > -1 && h_arglow[low_j] != -1){
            j0 = h_arglow[low_j];
            left_to_right_host(j0, j, h_rows_mp, h_aux_mp, h_low, m, p);
            low_j = h_low[j];
        }
        low_j = h_low[j];
        if (low_j > -1){
            h_arglow[low_j] = j;
        }
    }
}

inline void update_classes_host(int *h_classes, int *h_low, int *h_arglow, int m){
    for (int tid=0; tid<m; tid++){
        if (h_arglow[tid] > -1){
            h_classes[h_arglow[tid]] = -1;
            h_classes[tid] = 1;
        }
    }
}

inline void ess_hat_host(int *h_essential_hat, int *h_low, int *h_arglow, int m){
    for (int j=0; j<m; j++){
        if (h_low[j] > -1){
            h_essential_hat[h_low[j]] = 0;
        }
        if (h_arglow[j] > -1){
            h_essential_hat[h_arglow[j]] = 0;
            h_essential_hat[j] = 0;
        }
    }
}

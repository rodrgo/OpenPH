
inline void left_to_right_neighbours_host(int h_pivot, int *h_pivots, int *h_neighbours, int *h_rows_mp, int *h_aux_mp, int *h_low, int *h_arglow, int *h_classes, int *h_ess, int m, int p){
    int l = 0;
    while(h_neighbours[l] != -1){
        int j = h_neighbours[l];
        // essential estimation
        if (h_low[j] > -1)
            h_ess[h_low[j]] = 0;
        h_ess[h_pivot] = 0;
        h_ess[h_low[h_pivot]] = 0;
        // left to right
        left_to_right_host(h_pivot, j, h_rows_mp, h_aux_mp, h_low, m, p);
        if ((h_low[j] != -1) && (h_pivots[h_low[j]] == -1 || j < h_pivots[h_low[j]]))
            h_pivots[h_low[j]] = j;
        l++;
        // Update classes and ess
        if (h_pivots[h_low[j]] == j){
            h_classes[j] = -1;
            h_classes[h_low[j]] = 1;
            h_ess[j] = 0;
            h_ess[h_low[j]] = 0;
        }
    }
}

inline void get_neighbours_position_host(int pivot, int *h_neighbours, int *h_low, int *h_dim_next, int m){
    int pos = 0;
    int col = h_dim_next[pivot];
    while((col != -1)){
        if (h_low[col] == h_low[pivot]){
            h_neighbours[pos] = col;
            pos++;
        }
        col = h_dim_next[col];
    }
}

inline void ph_row(int *h_low, int *h_arglow, int *h_classes,
        int *h_ess, int *h_rows_mp, const int m, const int p,
        int *h_dim, int *h_dim_order, int *h_dim_next, int *h_dim_start,
        int *h_aux_mp, int *h_low_true, int *h_ess_true, float * h_float_m,
        float *error_lone, float *error_linf, float *error_redu,
        float *error_ess, float *time_track, int *p_iter){

    // time
    float time = 0.0;

    // iter and trackers
    track_host(0, m, h_low, h_ess, h_classes,
            h_low_true, h_ess_true, h_float_m,
            error_lone, error_linf, error_redu,
            error_ess, time_track, time);

    // d_is_neighbour
    int *h_neighbours = (int*)malloc( m * sizeof(int) );
    for (int i = 0; i < m; i++) h_neighbours[i] = -1;

    int *h_pivots = (int*)malloc( m * sizeof(int) );
    for (int i = 0; i < m; i++) h_pivots[i] = -1;
    for (int j = 0; j < m; j++)
        if ((h_low[j] != -1) && (h_pivots[h_low[j]] == -1 || j < h_pivots[h_low[j]])){
                h_pivots[h_low[j]] = j;
        }

    int iter = 1;
    for(int i = m-1; i > -1; i--){

        // TIC
        //clock_t tic = clock();
        cudaEvent_t start, stop;
        tic(&start, &stop);

        // Mark neighbours
        // Contrary to parallel case, we store indices of neighbours
        int pivot = h_pivots[i];
        if (pivot != -1){
            get_neighbours_position_host(pivot, h_neighbours, h_low, h_dim_next, m);
            
            // Reduce neighbours
            left_to_right_neighbours_host(pivot, h_pivots, h_neighbours, h_rows_mp, h_aux_mp, h_low, h_arglow, h_classes, h_ess, m, p);

            int l = 0;
            while (l<m && h_neighbours[l] != -1)
                h_neighbours[l++] = -1;

        }

        // TOC
        toc(start, stop, &time);

        // If last iteration, we found all essentials
        if(i == 0){
            for (int jj = 0; jj < m; jj++)
                if (h_low[jj] > -1){
                    h_ess[jj] = 0;
                    h_ess[h_low[jj]] = 0;
                    h_classes[jj] = -1;
                    h_classes[h_low[jj]] = 1;
                }
            for (int jj = 0; jj < m; jj++)
                if (h_ess[jj] ==  1 && h_classes[jj] == 0){
                    h_classes[jj] = 2;
                }
        }

        // meausre progress
        track_host(iter, m, h_low, h_ess, h_classes,
                h_low_true, h_ess_true, h_float_m,
                error_lone, error_linf, error_redu,
                error_ess, time_track, time);

        // iter
        iter++;

    } 
    p_iter[0] = iter;

    free(h_neighbours);
    free(h_pivots);

}


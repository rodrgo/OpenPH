inline void standard(int *h_low, int *h_arglow, int *h_classes, 
        int *h_ess, int *h_rows_mp, const int m, const int p, 
        int *h_aux_mp, int *h_low_true, int *h_ess_true, 
        float *h_float_m, float *error_lone,
        float *error_linf, float *error_redu, float *error_ess, 
        float *time_track, int *p_iter){

    // time
    float time = 0.0;

    // iter and trackers
    track_host(0, m, h_low, h_ess, h_classes, 
            h_low_true, h_ess_true, h_float_m, 
            error_lone, error_linf, error_redu, 
            error_ess, time_track, time);

    int iter = 1;
    for(int j = 0; j < m; j++){

        // TIC
        //clock_t tic = clock();
        cudaEvent_t start, stop;
        tic(&start, &stop);

        // Work on column "j"
        reduce_col_host(j, h_rows_mp, h_aux_mp, h_low, h_arglow, m, p, h_ess);

        // Update classes host
        if (h_low[j] > -1){
            h_classes[j] = -1;
            h_classes[h_low[j]] = 1;
        }else{
            h_classes[j] = 1;
        }

        // Essential estimation
        if (h_low[j] > -1){
            h_ess[j] = 0;
            h_ess[h_low[j]] = 0;
        }

        // TOC
        toc(start, stop, &time);

        // meausre progress
        track_host(iter, m, h_low, h_ess, h_classes, 
                h_low_true, h_ess_true, h_float_m, 
                error_lone, error_linf, error_redu, 
                error_ess, time_track, time);

        // iter
        iter++;

    }
    p_iter[0] = iter;

}

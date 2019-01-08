inline void standard_parallel(int *d_low, int *d_arglow, int *d_classes, 
        int *d_ess, int *d_rows_mp, const int m, const int p, 
        int *d_aux_mp, int *d_low_true, int *d_ess_true, 
        float *d_float_m, float *error_lone,
        float *error_linf, float *error_redu, float *error_ess, 
        float *time_track, int *p_iter, dim3 NBm, dim3 TPBm){

    // time
    float time = 0.0;

    // iter and trackers
    track(0, m, d_low, d_ess, d_classes, 
            d_low_true, d_ess_true, d_float_m, 
            error_lone, error_linf, error_redu, 
            error_ess, time_track, time, NBm, TPBm);

    int iter = 1;
    for(int j = 0; j < m; j++){

        // TIC
        cudaEvent_t start, stop;
        tic(&start, &stop);

        // Work on column "j"
        reduce_col<<<NBm, TPBm>>>(j, d_rows_mp, d_aux_mp, 
                d_low, d_arglow, m, p);
        cudaDeviceSynchronize();

        update_classes<<<NBm, TPBm>>>(d_classes, d_low, d_arglow, m);
        cudaDeviceSynchronize();

        // Essential estimation
        ess_hat<<<NBm, TPBm>>>(d_ess, d_low, d_arglow, m);
        cudaDeviceSynchronize();

        // TOC
        toc(start, stop, &time);

        // meausre progress
        track(iter, m, d_low, d_ess, d_classes, 
                d_low_true, d_ess_true, d_float_m, 
                error_lone, error_linf, error_redu, 
                error_ess, time_track, time, NBm, TPBm);

        // iter
        iter++;

    }
    p_iter[0] = iter;

}

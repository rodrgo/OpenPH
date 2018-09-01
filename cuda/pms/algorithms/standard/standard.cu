inline void standard(int *d_low, int *d_arglow, int *d_classes, 
        int *d_ess, int *d_rows_mp, const int m, const int p, 
        int *d_aux_mp, int *d_low_true, int *d_ess_true, float *error_lone, 
        float *error_linf, float *error_redu, float *error_ess, float *time_track,
        int *p_iter, dim3 NBm, dim3 TPBm){

    int iter = 1;
    for(int j = 0; j < m; j++){

        // Create timing
        float time;
        cudaEvent_t start, stop;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        // Reduce column
        reduce_col<<<NBm, TPBm>>>(j, d_rows_mp, d_aux_mp, d_low, d_arglow, m, p);
        cudaDeviceSynchronize();

        update_classes<<<NBm, TPBm>>>(d_classes, d_low, d_arglow, m);
        cudaDeviceSynchronize();

        // end timing
        cudaThreadSynchronize();
        cudaEventRecord(stop, 0);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&time, start, stop);
        cudaEventDestroy(start);
        cudaEventDestroy(stop);

        // Essential estimation
        ess_hat<<<NBm, TPBm>>>(d_ess, d_low, d_arglow, m);
        cudaDeviceSynchronize();

        // iter and trackers
        track(iter, m, d_low, d_classes, d_low_true, d_ess_true, error_lone, error_linf, error_redu, error_ess, time_track, time, NBm, TPBm);

        iter++;

    }
    p_iter[0] = iter;

}

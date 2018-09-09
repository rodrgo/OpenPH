inline void pms(int *d_low, int *d_arglow, int *d_classes, int *d_ess,
        int *d_rows_mp, const int m, const int p,
        int *d_dims, int *d_dims_order, int *d_dims_order_next, int *d_dims_order_start,
        int complex_dimension, int *d_left, int *d_beta, int *d_aux_mp,
        int *d_low_true, int *d_ess_true, float *d_float_m,
        float *error_lone, float *error_linf, float *error_redu, float *error_ess,
        float *time_track, int *p_iter,
        dim3 NBm, dim3 TPBm, dim3 NBcdim, dim3 TPBcdim){
    //  d_pivots[j] = 1  <=> d_arglow[d_low[j]] = j 

    // time
    float time = 0.0;

    // iter and trackers
    track(0, m, d_low, d_ess, d_classes,
            d_low_true, d_ess_true, d_float_m,
            error_lone, error_linf, error_redu,
            error_ess, time_track, time, NBm, TPBm);

    // Auxiliary variables
    int *d_aux;
    int *d_is_positive;
    cudaMalloc((void**)&d_aux, m * sizeof(int));
    cudaMalloc((void**)&d_is_positive, m * sizeof(int));

    // -----------------------
    // Do some pre-processing work
    // -----------------------

    // Compute simplex dimensions (on device)
    // Get maximum dimension 
    int cdim = complex_dimension + 2; // -1, 0, 1, 2, ..., complex_dim
    int *d_aux_cdim;    // Auxiliary vector of size cdim 
    cudaMalloc((void**)&d_aux_cdim, cdim * sizeof(int));

    // -----------------------
    // Phase 0
    // -----------------------

    // Mark pivots and clear corresponding positives
    mark_pivots_and_clear<<<NBm, TPBm>>>(d_low, d_beta, 
            d_classes, d_rows_mp, d_arglow, m, p);
    cudaDeviceSynchronize();

    int converged = is_reduced(d_aux, d_low, m, NBm, TPBm);
    int iter = 1;

    while (! converged ){

        // TIC
        cudaEvent_t start, stop;
        tic(&start, &stop);

        // -----------------------
        // Main iteration : Phase I 
        // -----------------------

        fill<<<NBm, TPBm>>>(d_aux, 0, m);
        fill<<<NBm, TPBm>>>(d_aux_cdim, -1, cdim); // d_ceil
        fill<<<NBm, TPBm>>>(d_is_positive, 0, m);
        cudaDeviceSynchronize();

        transverse_dimensions<<<NBcdim, TPBcdim>>>(d_dims, 
                d_dims_order, d_dims_order_next, d_dims_order_start, 
                d_low, d_arglow, d_classes, d_is_positive,   
                d_aux, d_aux_cdim, cdim);
        cudaDeviceSynchronize();

        clear_positives<<<NBm, TPBm>>>(d_is_positive, 
                d_low, d_classes, d_rows_mp, m, p);
        cudaDeviceSynchronize();

        // -----------------------
        // Main iteration : Phase II 
        // -----------------------

        fill<<<NBm, TPBm>>>(d_is_positive, 0, m);
        cudaDeviceSynchronize();

        phase_ii<<<NBm, TPBm>>>(d_low, d_left, d_classes, 
                d_is_positive, d_arglow, d_rows_mp, d_aux_mp, m, p);
        cudaDeviceSynchronize();

        clear_positives<<<NBm, TPBm>>>(d_is_positive, 
                d_low, d_classes, d_rows_mp, m, p);
        cudaDeviceSynchronize();

        // Check again if its reduced
        converged = is_reduced(d_aux, d_low, m, NBm, TPBm);

        // update classes (Not necessary for algo to work)
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

    set_unmarked<<<NBm, TPBm>>>(d_classes, d_low, d_arglow, d_rows_mp, m, p);
    cudaDeviceSynchronize();

    p_iter[0] = iter;

    cudaFree(d_aux_cdim);
    cudaFree(d_aux);
    cudaFree(d_is_positive);

}


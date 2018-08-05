inline void pms(int *d_rows_mp, int *d_aux_mp, int *d_low, int *d_arglow, int *d_dims, int *d_dims_order, int *d_dims_order_next, int *d_dims_order_start, const int m, const int p, int complex_dimension, int *d_beta, float *resRecord, float *timeRecord, int *p_iter, dim3 NBm, dim3 TPBm, dim3 NBcdim, dim3 TPBcdim){

    // Auxiliary variables
    int *d_aux;
    int *d_clear;
    cudaMalloc((void**)&d_aux, m * sizeof(int));
    cudaMalloc((void**)&d_clear, m * sizeof(int));

    // -----------------------
    // Do some pre-processing work
    // -----------------------

    // d_classes
    int *d_classes;
    cudaMalloc((void**)&d_classes, m * sizeof(int));
    fill<<<NBm, TPBm>>>(d_classes, 0, m);

    // Compute simplex dimensions (on device)
    // Get maximum dimension (TODO: Change Twist code too)
    int cdim = complex_dimension + 2; // -1, 0, 1, 2, ..., complex_dim
    int *d_aux_cdim;    // Auxiliary vector of size cdim 
    cudaMalloc((void**)&d_aux_cdim, cdim * sizeof(int));
    printf("passed compute_dimension_order!!\n");
    printf("cdim = %d\n", cdim);

    // locks
    int *d_locks_cdim;
    cudaMalloc((void**)&d_locks_cdim, cdim * sizeof(int));

    // updated

    // -----------------------
    // Phase 0
    // -----------------------

    printvec(d_low, m, "d_low");
    printvec(d_beta, m, "d_beta");

    alpha_beta_reduce<<<NBm, TPBm>>>(d_low, d_beta, d_classes, 
            d_rows_mp, d_arglow, m, p);

    printvec(d_classes, m, "d_classes");
    printvec(d_arglow, m, "d_arglow");

    int converged = is_reduced(d_aux, d_low, m, NBm, TPBm);
    printf("converged %d\n", converged);

    while (! converged ){

        // -----------------------
        // Main iteration : Phase I 
        // -----------------------

        fill<<<NBm, TPBm>>>(d_aux, 0, m);
        fill<<<NBm, TPBm>>>(d_locks_cdim, 0, cdim); // d_next_cdim
        fill<<<NBm, TPBm>>>(d_aux_cdim, -1, cdim); // d_ceil

        printvec(d_dims, m, "d_dims");
        printvec(d_dims_order, m, "d_dims_order");
        printvec(d_low, m, "d_low");
        printvec(d_arglow, m, "d_arglow");
        printvec(d_classes, m, "d_classes");
        printvec(d_locks_cdim, cdim, "d_next_cdim");
        printvec(d_aux_cdim, cdim, "d_ceil");

        phase_i_cdim<<<NBcdim, TPBcdim>>>(d_dims, d_dims_order,
                d_dims_order_next, d_dims_order_start, 
                d_low, d_arglow, d_classes, d_clear,   
                d_aux, d_aux_cdim, d_locks_cdim, m, cdim);
        /*
        phase_i<<<NBm, TPBm>>>(d_dims, d_dims_order, 
                d_low, d_arglow, d_classes, d_clear,   
                d_aux, d_aux_cdim, d_locks_cdim, m);
        */
        cudaDeviceSynchronize();
        exit(0);

        // -----------------------
        // Main iteration : Phase II 
        // -----------------------

        phase_ii<<<NBm, TPBm>>>(d_low, d_beta, d_classes, 
                d_arglow, d_rows_mp, d_aux_mp, m, p);

        // record iteration
        //record_iteration();

        // Check again if its reduced
        converged = is_reduced(d_aux, d_low, m, NBm, TPBm);
    }

    set_unmarked<<<NBm, TPBm>>>(d_classes, d_low, d_arglow, 
            d_rows_mp, m, p);

    cudaFree(d_locks_cdim);
    cudaFree(d_aux_cdim);
    cudaFree(d_aux);
    cudaFree(d_clear);
    cudaFree(d_classes);

}


inline void pms(int *d_rows_mp, int *d_aux_mp, int *d_low, int *d_arglow, const int m, const int p, int *d_beta, float *resRecord, float *timeRecord, int *p_iter, dim3 NBm, dim3 TPBm){

    // Auxiliary variables
    int *d_aux;
    cudaMalloc((void**)&d_aux, m * sizeof(int));

    // -----------------------
    // Do some pre-processing work
    // -----------------------

    // d_classes
    int *d_classes;
    cudaMalloc((void**)&d_classes, m * sizeof(int));
    fill<<<NBm, TPBm>>>(d_classes, 0, m);

    // Compute simplex dimensions (on device)
    // Get maximum dimension (TODO: Change Twist code too)
    int *d_dims;        // Dimension of simplex j
    int *d_dims_order;  // j is "d_dims_order[j]"-th simplex in dim_j
    int *d_aux_cdim;    // Auxiliary vector of size cdim 
    int complex_dim;    // max(d_dims), dimension of complex

    cudaMalloc((void**)&d_dims, m * sizeof(int));
    cudaMalloc((void**)&d_dim_pos, m * sizeof(int));
    compute_simplex_dimensions(d_dims, d_dims_order, &complex_dim, 
            d_rows_mp, m, p, NBm, TPBm);

    int cdim = complex_dim + 1;

    cudaMalloc((void**)&d_aux_cdim, cdim * sizeof(int));
    compute_dimension_order(d_dims, d_dims_order, d_aux_cdim, 
            cdim, m, NBm, TPBm);

    // ceilings
    int *d_ceilings;
    cudaMalloc((void**)&d_ceilings, m * sizeof(int));
    fill<<<NBm, TPBm>>>(d_ceilings, 0, m);

    // locks
    int *d_locks_cdim;
    cudaMalloc((void**)&d_locks_cdim, cdim * sizeof(int));

    // updated

    // -----------------------
    // Phase 0
    // -----------------------

    alpha_beta_reduce<<<NBm, TPBm>>>(d_lows, d_beta, d_classes, 
            d_rows_mp, d_arglow, d_lowstar, m);

    int converged = is_reduced(d_aux, d_lows, m, NBm, TPBm);
    while (! converged ){

        // -----------------------
        // Main iteration : Phase I 
        // -----------------------

        // TODO: In get_ceilings, atomicCAS needs pointer in first
        // argument. Check if this is being given correctly
        // TODO: atomicMAX?
        fill<<<NBm, TPBm>>>(d_aux, 0, m);
        fill<<<NBm, TPBm>>>(d_locks_cdim, 0, cdim);
        fill<<<NBm, TPBm>>>(d_aux_cdim, 0, cdim);
        phase_i<<<NBm, TPBm>>>(d_ceilings, d_dims, d_dims_order, 
                d_low, d_arglow, d_classes, d_clear,   
                d_aux, d_aux_cdim, d_locks_cdim, m);

        // -----------------------
        // Main iteration : Phase II 
        // -----------------------

        phase_ii<<<NBm, TPBm>>>(d_low, d_arglow, d_rows_mp, 
                d_aux_mp, m, p);

        // record iteration
        record_iteration();

        // Check again if its reduced
        converged = is_reduced(d_aux, d_lows, m, NBm, TPBm);
    }

    set_unmarked<<<NBm, TPBm>>>(d_classes, d_low, d_arglow, 
            d_rows_mp, m, p);

    /*
    int *d_dim_count;
    cudaMalloc((void**)&d_dim_count, cdim * sizeof(int)); // [0, 1, ..., complex_dim]
    fill<<<NBm, TPBm>>>(d_dim_count, 0, cdim);
    count_simplices_dim<<<NBm, TPBm>>>(d_dim_count, d_dims);
    */

    cudaFree(d_ceilings);
    cudaFree(d_dims);
    cudaFree(d_dims_order);
    cudaFree(d_locks_cdim);
    cudaFree(d_aux_cdim);
    cudaFree(d_aux);
    cudaFree(d_classes);

}


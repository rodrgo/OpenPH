inline void pms(int *d_rows_mp, int *d_aux_mp, int *d_low, int *d_arglow, int *d_dims, int *d_dims_order, int *d_dims_order_next, int *d_dims_order_start, const int m, const int p, int complex_dimension, int *d_left, int *d_beta, float *resRecord, float *timeRecord, int *p_iter, dim3 NBm, dim3 TPBm, dim3 NBcdim, dim3 TPBcdim){
    //  d_pivots[j] = 1  <=> d_arglow[d_low[j]] = j 

    // Auxiliary variables
    int *d_aux;
    int *d_clear;
    cudaMalloc((void**)&d_aux, m * sizeof(int));
    cudaMalloc((void**)&d_clear, m * sizeof(int));

    // -----------------------
    // Do some pre-processing work
    // -----------------------

    // d_pivots (Binary)
    int *d_pivots;
    cudaMalloc((void**)&d_pivots, m * sizeof(int));
    fill<<<NBm, TPBm>>>(d_pivots, 0, m);

    // d_classes
    int *d_classes;
    cudaMalloc((void**)&d_classes, m * sizeof(int));
    fill<<<NBm, TPBm>>>(d_classes, 0, m);

    // Compute simplex dimensions (on device)
    // Get maximum dimension (TODO: Change Twist code too)
    int cdim = complex_dimension + 2; // -1, 0, 1, 2, ..., complex_dim
    int *d_aux_cdim;    // Auxiliary vector of size cdim 
    cudaMalloc((void**)&d_aux_cdim, cdim * sizeof(int));
    printf("cdim = %d\n", cdim);

    // locks
    int *d_locks_cdim;
    cudaMalloc((void**)&d_locks_cdim, cdim * sizeof(int));

    // -----------------------
    // Phase 0
    // -----------------------

    // only mark pivots. Do clearing in next routine
    mark_pivots<<<NBm, TPBm>>>(d_pivots, d_low, d_beta, d_classes, 
            d_rows_mp, d_arglow, m, p);

    printvec(d_classes, m, "d_classes");
    printvec(d_arglow, m, "d_arglow");

    int converged = is_reduced(d_aux, d_low, m, NBm, TPBm);
    printf("converged %d\n", converged);

    int iter = 0;
    thrust::device_ptr<int> d_classes_ptr = thrust::device_pointer_cast(d_classes);
    int num_zeros = thrust::count(d_classes_ptr, d_classes_ptr + m, 0);
    printf("num_zeros=%d\n", num_zeros);

    while (! converged ){

        // -----------------------
        // Main iteration : Phase I 
        // -----------------------

        fill<<<NBm, TPBm>>>(d_aux, 0, m);
        fill<<<NBm, TPBm>>>(d_locks_cdim, 0, cdim); // d_next_cdim
        fill<<<NBm, TPBm>>>(d_aux_cdim, -1, cdim); // d_ceil
        fill<<<NBm, TPBm>>>(d_clear, 0, m);
        cudaDeviceSynchronize();

        phase_i_cdim<<<NBcdim, TPBcdim>>>(d_dims, d_dims_order,
                d_dims_order_next, d_dims_order_start, 
                d_low, d_arglow, d_classes, d_clear,   
                d_aux, d_aux_cdim, d_locks_cdim, m, cdim);
        cudaDeviceSynchronize();

        clear_phase_i<<<NBm, TPBm>>>(d_low, d_classes, d_rows_mp, d_clear, m, p);
        cudaDeviceSynchronize();

        // -----------------------
        // Main iteration : Phase II 
        // -----------------------

        phase_ii<<<NBm, TPBm>>>(d_low, d_left, d_classes, 
                d_arglow, d_rows_mp, d_aux_mp, m, p);
        cudaDeviceSynchronize();

        // Check again if its reduced
        converged = is_reduced(d_aux, d_low, m, NBm, TPBm);
        printf("converged %d\n", converged);

        count_nonzeros<<<NBm, TPBm>>>(d_classes, m);
        cudaDeviceSynchronize();

        num_zeros = thrust::count(d_classes_ptr, d_classes_ptr + m, 0);
        //cudaMemcpyFromSymbol(&num_zeros, "d_atomicsum", sizeof(int), 0, cudaMemcpyDeviceToHost);
        
        printf("num_zeros = %d\n", num_zeros);
        // record iteration
        //record_iteration();
        iter++;
        if (iter > 4)
            converged = 1;

    }

    set_unmarked<<<NBm, TPBm>>>(d_classes, d_low, d_arglow, 
            d_rows_mp, m, p);

    cudaFree(d_locks_cdim);
    cudaFree(d_aux_cdim);
    cudaFree(d_aux);
    cudaFree(d_clear);
    cudaFree(d_classes);
    cudaFree(d_pivots);

}


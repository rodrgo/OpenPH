inline void algorithm_factory(char *algstr,
        int *d_low, int *d_arglow, int *d_classes, int *d_ess,
        int *d_rows, int *d_dim, int *d_dim_order, int *d_dim_next, int *d_dim_start,
        int *d_beta, int *d_left,
        int m, int p, int complex_dim, 
        int *d_aux_mp, int *d_low_true, int *d_ess_true, float *d_float_m,
        float *error_lone, float *error_linf, float *error_redu, float *error_ess, float *time_track, int *p_iter,
        dim3 NBm, dim3 TPBm, dim3 NBcdim, dim3 TPBcdim){

    int is_serial = ((strcmp(algstr, "standard") == 0) || \
        (strcmp(algstr, "twist") == 0) || \
		(strcmp(algstr, "ph_row") == 0)) ? 1 : 0;

    if (is_serial){

        int mp = m * p;
        int cdim = complex_dim + 2;

        // Copy data from device to host
        int *h_low      = (int*)malloc( sizeof(int)*m );
        int *h_arglow   = (int*)malloc( sizeof(int)*m );
        int *h_classes  = (int*)malloc( sizeof(int)*m );
        int *h_ess      = (int*)malloc( sizeof(int)*m );
        int *h_rows_mp  = (int*)malloc( sizeof(int)*mp );
        int *h_aux_mp   = (int*)malloc( sizeof(int)*mp );
        int *h_low_true = (int*)malloc( sizeof(int)*m );
        int *h_ess_true = (int*)malloc( sizeof(int)*m );

        int *h_dim       = (int*)malloc( m * sizeof(int) );
        int *h_dim_order = (int*)malloc( m * sizeof(int) );
        int *h_dim_next  = (int*)malloc( m * sizeof(int) );
        int *h_dim_start = (int*)malloc( cdim * sizeof(int) );

        float *h_float_m = (float*)malloc( sizeof(float)*m );

        fill_host(h_arglow, -1, m);
        fill_host(h_classes, 0, m);
        fill_host(h_ess, 1, m);
        fill_host(h_aux_mp, -1, mp);

        cudaMemcpy(h_dim, d_dim_start, m*sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_dim_order, d_dim_order, m*sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_dim_next, d_dim_next, m*sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_dim_start, d_dim_start, cdim*sizeof(int), cudaMemcpyDeviceToHost);

        cudaMemcpy(h_low, d_low, sizeof(int)*m, cudaMemcpyDeviceToHost);
        cudaMemcpy(h_rows_mp, d_rows, sizeof(int)*mp, cudaMemcpyDeviceToHost);
        cudaMemcpy(h_low_true, d_low_true, sizeof(int)*m, cudaMemcpyDeviceToHost);
        cudaMemcpy(h_ess_true, d_ess_true, sizeof(int)*m, cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();

        // Start serial algorithms

        if (strcmp(algstr, "standard")==0){
            standard(h_low, h_arglow, h_classes, h_ess, 
                    h_rows_mp, m, p, h_aux_mp, h_low_true,
                    h_ess_true, h_float_m,
                    error_lone, error_linf, error_redu, 
                    error_ess, time_track, p_iter);
        }else if(strcmp(algstr, "twist")==0){
            twist(h_low, h_arglow, h_classes, h_ess,
                    h_rows_mp, m, p, h_dim, h_dim_order,
                    h_dim_next, h_dim_start,
                    complex_dim, h_aux_mp, h_low_true,
                    h_ess_true, h_float_m,
                    error_lone, error_linf, error_redu,
                    error_ess, time_track, p_iter);
        }else if(strcmp(algstr, "ph_row")==0){
            ph_row(h_low, h_arglow, h_classes, h_ess,
                    h_rows_mp, m, p,
                    h_dim, h_dim_order, h_dim_next, h_dim_start,
                    h_aux_mp, h_low_true, h_ess_true, h_float_m,
                    error_lone, error_linf, error_redu, error_ess,
                    time_track, p_iter);

        }else{
            printf("Not recognised");
        }

        // Export h_low and h_ess
        cudaMemcpy(d_low, h_low, sizeof(int)*m, cudaMemcpyHostToDevice);
        cudaMemcpy(d_ess, h_ess, sizeof(int)*m, cudaMemcpyHostToDevice);
        cudaDeviceSynchronize();

        // Free pointers
        free(h_low);
        free(h_arglow);
        free(h_classes); 
        free(h_ess);
        free(h_rows_mp);
        free(h_aux_mp);
        free(h_low_true);
        free(h_ess_true);

        free(h_dim);
        free(h_dim_order);
        free(h_dim_next);
        free(h_dim_start);

        free(h_float_m);

    }else{
        if (strcmp(algstr, "standard_parallel")==0){
            standard_parallel(d_low, d_arglow, d_classes, d_ess, 
                    d_rows, m, p, d_aux_mp, d_low_true,
                    d_ess_true, d_float_m,
                    error_lone, error_linf, error_redu, 
                    error_ess, time_track, p_iter, NBm, 
                    TPBm); 
        } else if (strcmp(algstr, "twist_parallel")==0){
            twist_parallel(d_low, d_arglow, d_classes, d_ess,
                    d_rows, m, p,
                    d_dim, d_dim_order, d_dim_next, d_dim_start,
                    complex_dim,
                    d_aux_mp, d_low_true, d_ess_true, d_float_m,
                    error_lone, error_linf, error_redu, error_ess,
                    time_track, p_iter, NBm, TPBm);
        } else if (strcmp(algstr, "ph_row_parallel")==0){
            ph_row_parallel(d_low, d_arglow, d_classes, d_ess,
                    d_rows, m, p,
                    d_aux_mp, d_low_true, d_ess_true, d_float_m,
                    error_lone, error_linf, error_redu, error_ess,
                    time_track, p_iter, NBm, TPBm);
        } else if (strcmp(algstr, "pms")==0){
            pms(d_low, d_arglow, d_classes, d_ess,
                    d_rows, m, p,
                    d_dim, d_dim_order, d_dim_next, d_dim_start,
                    complex_dim, d_left, d_beta, d_aux_mp,
                    d_low_true, d_ess_true, d_float_m,
                    error_lone, error_linf, error_redu, error_ess,
                    time_track, p_iter, NBm, TPBm, NBcdim, TPBcdim);
        }else
            printf("Not recognised");
    }

}

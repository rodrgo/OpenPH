inline void twist(int *d_rows_mp, int *d_aux_mp, int *d_low, int *d_arglow, const int m, const int p, float *resRecord, float *timeRecord, int *p_iter, dim3 numBlocks_m, dim3 threadsPerBlock_m){

    // Compute initial dimensions (on device)
    int *d_simplex_dimensions;
    cudaMalloc((void**)&d_simplex_dimensions, m * sizeof(int));
    get_simplex_dimensions<<<numBlocks_m, threadsPerBlock_m>>>(d_simplex_dimensions, d_rows_mp, m, p);

    // Copy to host
    int *h_simplices_dim;
    h_simplices_dim = (int*)malloc( m * sizeof(int) );
    cudaMemcpy(h_simplices_dim, d_simplex_dimensions, m*sizeof(int), cudaMemcpyDeviceToHost);

    // Get dimension of complex (maximum initial dimension)
    int complex_dim = 0;
    for (int i = 0; i < m; i++){
        complex_dim = h_simplices_dim[i] > complex_dim ? h_simplices_dim[i] : complex_dim;
    }

    int *d_simplices_dim;
    cudaMalloc((void**)&d_simplices_dim, m * sizeof(int));

    for(int dim = complex_dim; dim > 0; dim--){

        // Bring simplices of dimension d to host 
        int zero = 0;
        cudaMemcpyToSymbol(d_lock, &zero, sizeof(int));
        cudaMemcpyToSymbol(d_pos,  &zero, sizeof(int));
        minusone_vector_int<<< numBlocks_m, threadsPerBlock_m >>>((int*)d_simplices_dim, m);
        get_simplices_dim<<<numBlocks_m, threadsPerBlock_m>>>(d_simplex_dimensions, d_simplices_dim, dim, m);

        // Bring d_simplices_dim back to host
        cudaMemcpy(h_simplices_dim, d_simplices_dim, m*sizeof(int), cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();

        int idx = 0;
        while (idx < m && h_simplices_dim[idx] > -1){
            int j = h_simplices_dim[idx]; 
            reduce_col_gpu<<<numBlocks_m, threadsPerBlock_m>>>(j, d_rows_mp, d_aux_mp, d_low, d_arglow, m, p);
            twist_step<<<numBlocks_m, threadsPerBlock_m>>>(j, d_rows_mp, d_low, m, p);
            cudaDeviceSynchronize();
            idx++;
        }
    } 

    free(h_simplices_dim);
    cudaFree(d_simplex_dimensions);
    cudaFree(d_simplices_dim);

}


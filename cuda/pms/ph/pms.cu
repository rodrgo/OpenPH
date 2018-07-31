inline void pms(int *d_rows_mp, int *d_aux_mp, int *d_low, int *d_arglow, const int m, const int p, int *d_beta, float *resRecord, float *timeRecord, int *p_iter, dim3 numBlocks_m, dim3 threadsPerBlock_m){

    // Check if matrix is reduced
    int *d_aux;
    cudaMalloc((void**)&d_aux, m * sizeof(int));

    // Returns 1 if matrix is reduced, or 0 otherwise
    is_reduced(d_aux, d_lows, m, numBlocks_m, threadsPerBlock_m);

    // Create classes vector
    int *d_classes;
    cudaMalloc((void**)&d_classes, m * sizeof(int));
    zero_vector_int<<<numBlocks_m, threadsPerBlock_m>>>(d_classes, m);

    // Compute initial dimensions (on device)
    // Get maximum dimension (TODO: Change Twist code too)
    int *d_simplex_dimensions;
    int complex_dim;
    cudaMalloc((void**)&d_simplex_dimensions, m * sizeof(int));
    compute_simplex_dimensions(d_simplex_dimensions, 
            d_rows_mp, m, p, &complex_dim, 
            numBlocks_m, threadsPerBlock_m);

    // Vector of ceilings?
    /*
       // need vector with order of simplex in each dimension
        if (tid < m){
            int dim_j = d_simplex_dims[tid];
            int order = d_


        }


    */


    // Alpha-Beta reduce
    alpha_beta_reduce<<<numBlocks_m, threadsPerBlock_m>>>(d_lows, d_beta, d_classes, d_rows_mp, d_arglow, d_lowstar, m);

    // mark_unique_unreduced_per_dimension
    /*
       1. For each dimension count how many entries the dimension has
       2. Create one vector per dimension of that size
       3. On kernel, use an atomicAdd to accumulate all the entries into the vector
       4. Sort each vector

       1. For each dimension get indices of columns with that dimension and store in vector
       2. Sort

       2. Parallelise over number of dimensions and run a sentinel across each list, keeping a "ceiling" variable for each dimension  
            2.1 If column meets condition, put an indicator variable
       3. In another kernel do the clearing and updating in parallel 
    */
    int *d_dim_count;
    cudaMalloc((void**)&d_dim_count, (complex_dim+1) * sizeof(int)); // [0, 1, ..., complex_dim]
    init_vector_int<<<numBlocks_m, threadsPerBlock_m>>>(d_dim_count, 0, complex_dim);
    count_simplices_dim<<<numBlocks_m, threadsPerBlock_m>>>(d_dim_count, d_simplex_dimensions, complex_dim);

    // Create array of arrays containing positions of simplices at each dimension
    // https://stackoverflow.com/questions/12406948/array-of-vectors-using-thrust

    /*
        The elements for dimension "i" are
        idx_start = sum(d_dim_count[j] : j < i), if i > 0 
        idx_start = 0, i == 0
        idx_end   = sum(d_dim_count[j] : j < i+1), i < complex_dim
        idx_end   = m, i == complex_dim
    */
    int *d_dim_pos; // ordered vector of positions for each dimension
    cudaMalloc((void**)&d_dim_pos, m * sizeof(int));
    zero_vector_int<<<numBlocks_m, threadsPerBlock_m>>>(d_dim_pos, m);
    get_dim_pos<<<numBlocks_m, threadsPerBlock_m>>>(d_simplex_dimensions, d_dim_count, d_dim_pos, complex_dim);

    // Continue
    int *d_simplices_dim;
    cudaMalloc((void**)&d_simplices_dim, m * sizeof(int));

    cudaFree(d_simplex_dimensions);
    cudaFree(d_simplices_dim);
    cudaFree(d_aux);
    cudaFree(d_classes);

}


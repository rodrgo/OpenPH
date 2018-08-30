
#include "init.cu"

// PMS
// Pass boundary matrix in row-column format 
// Convert to GPU
// Compute low and lowstar in GPU

// Host function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

    if ( nlhs!=3 & nrhs != 5)
        printf("[Error nlhs not 3]");
    else {

        // -------------------------------
        // Define I/O
        // -------------------------------

        // inputs
        int *h_rows, *h_cols, *h_vals;
        int m;

        // outputs
        int *h_low;
        float *resRecord, *timeRecord;

        // -------------------------------
        // Read Inputs
        // -------------------------------

        // reading in the string to determine the algorithm
        int strlen = mxGetN(prhs[0])+1;
        char algstr[strlen+100];
        int algerr = mxGetString(prhs[0], algstr, strlen);

        // Read data
        h_rows = (int*)mxGetData(prhs[1]);
        h_cols = (int*)mxGetData(prhs[2]);
        h_vals = (int*)mxGetData(prhs[3]);   // Ignore h_vals
        m      = (int)mxGetScalar(prhs[4]);

        // number of non-zeros
        int nnz;
        nnz = mxGetM(prhs[3]);

        // -------------------------------
        // Set value of p
        // -------------------------------

        // Assert col order
        int is_col_order = assert_col_order(h_cols, m, nnz);
        // printf("is_col_order = %d\n", is_col_order);

        // Count nonzeros in columns and rows
        int max_nnz_rows = get_max_nnz(h_rows, m, nnz);
        int max_nnz_cols = get_max_nnz(h_cols, m, nnz);

        //DEBUG
        if (1 == 0){
            printf("max_nnz_rows : %d\n", max_nnz_rows);
            printf("max_nnz_cols : %d\n", max_nnz_cols);
        }
          
        int p = 5*max_nnz_cols;
        int mp = m * p;

        //DEBUG
        if (1 == 0){
            printf("p = %d\n", p);
        }

        // -------------------------------
        // GPU
        // -------------------------------

        int gpu_number = 3;
        int threads_perblock_m  = 0;
        int threads_perblock_nnz = 0;
        int threads_perblock_mp = 0;
        set_gpu_device(gpu_number, &threads_perblock_m, &threads_perblock_nnz, &threads_perblock_mp, m, nnz, p);

        //DEBUG
        if (1 == 0){
            printf("threads_perblock_m:   %d\n", threads_perblock_m);
            printf("threads_perblock_nnz: %d\n", threads_perblock_nnz);
            printf("threads_perblock_mp:  %d\n", threads_perblock_mp);
        }

        dim3 threadsPerBlock_nnz(threads_perblock_nnz);
        int num_blocks_nnz = (int)ceil((float)(nnz)/(float)threads_perblock_nnz);
        dim3 numBlocks_nnz(num_blocks_nnz);

        dim3 threadsPerBlock_m(threads_perblock_m);
        int num_blocks_m = (int)ceil((float)m/(float)threads_perblock_m);
        dim3 numBlocks_m(num_blocks_m);

        dim3 threadsPerBlock_mp(threads_perblock_mp);
        int num_blocks_mp = (int)ceil((float)mp/(float)threads_perblock_mp);
        dim3 numBlocks_mp(num_blocks_mp);

        // -------------------------------
        // Outputs
        // -------------------------------
        
        // Create output
        plhs[0] = mxCreateNumericMatrix(1, m, mxINT32_CLASS, mxREAL);
        h_low  = (int*) mxGetData(plhs[0]);  

        plhs[1] = mxCreateNumericMatrix((m + 1), 1, mxSINGLE_CLASS, mxREAL);
        resRecord = (float*) mxGetData(plhs[1]);

        plhs[2] = mxCreateNumericMatrix((m + 1), 1, mxSINGLE_CLASS, mxREAL);
        timeRecord = (float*) mxGetData(plhs[2]);

        // initialise options at default values
        unsigned int seed = clock();

        // -------------------------------
        // Create data on device
        // -------------------------------

        // d_rows, d_cols
        int *d_rows, *d_cols;

        cudaMalloc((void**)&d_rows, nnz * sizeof(int));
        cudaMalloc((void**)&d_cols, nnz * sizeof(int));

        cudaMemcpy(d_rows, h_rows, sizeof(int)*nnz, cudaMemcpyHostToDevice);
        cudaDeviceSynchronize();
        cudaMemcpy(d_cols, h_cols, sizeof(int)*nnz, cudaMemcpyHostToDevice);
        cudaDeviceSynchronize();

        indexShiftDown<<<numBlocks_nnz,threadsPerBlock_nnz>>>(d_rows, nnz);
        indexShiftDown<<<numBlocks_nnz,threadsPerBlock_nnz>>>(d_cols, nnz); 
        cudaDeviceSynchronize();

        // d_low, d_arglow
        int *d_low, *d_arglow;

        cudaMalloc((void**)&d_low,    m * sizeof(int));
        cudaMalloc((void**)&d_arglow, m * sizeof(int));

        fill<<< numBlocks_m, threadsPerBlock_m >>>(d_low, -1, m);
        fill<<< numBlocks_m, threadsPerBlock_m >>>(d_arglow, -1, m);
        cudaDeviceSynchronize();

        // d_aux_mp, d_rows_mp
        int *d_rows_mp, *d_aux_mp;

        cudaMalloc((void**)&d_aux_mp,   mp * sizeof(int));
        cudaMalloc((void**)&d_rows_mp,  mp * sizeof(int));

        fill<<< numBlocks_mp, threadsPerBlock_mp >>>(d_aux_mp, -1, mp);
        fill<<< numBlocks_mp, threadsPerBlock_mp >>>(d_rows_mp, -1, mp);
        create_rows_mp<<< numBlocks_nnz, threadsPerBlock_nnz>>>(d_rows, d_cols, d_rows_mp, m, p, nnz);
        cudaDeviceSynchronize();

        compute_low_mp<<<numBlocks_m,threadsPerBlock_m>>>(d_rows_mp, d_low, m, p);
        cudaDeviceSynchronize();

        //DEBUG
        if (1 == 0){
            printvec(d_arglow, m, "ARGLOW");
            printvec(d_low, m, "LOW");
        }

        // DEBUG
        if (1 == 0){
            print_matrix_cols(d_rows, d_cols, nnz);
            print_matrix_cols_mp(d_rows_mp, m, p);
        }

        // -------------------------------
        // Get PH vectors
        // -------------------------------

        cudaEvent_t startIHT, stopIHT;
        float timeIHT;
        cudaEventCreate(&startIHT);
        cudaEventCreate(&stopIHT);
        cudaEventRecord(startIHT,0);

        // -------------------------------
        // simplex dimensions, complex dimension, dimension order
        // -------------------------------
        
        // d_dims: Dim of simplex j (-1, 0, 1, ..., complex_dim)
        // j is "d_dims_order[j]"-th simplex in dim_j
        // max(d_dims), dimension of complex

        int *d_dims, *d_dims_order; 
        int *d_dims_order_next;
        int complex_dim = -1;
        cudaMalloc((void**)&d_dims, m * sizeof(int));
        cudaMalloc((void**)&d_dims_order, m * sizeof(int));
        cudaMalloc((void**)&d_dims_order_next, m * sizeof(int));

        compute_simplex_dimensions_h(h_cols, m, p, nnz,
                d_dims, d_dims_order, d_dims_order_next, &complex_dim);
        int cdim = complex_dim + 2;
        int threads_perblock_cdim = min(cdim, threads_perblock_m);
        dim3 threadsPerBlock_cdim(threads_perblock_cdim);
        int num_blocks_cdim = (int)ceil((float)cdim/(float)threads_perblock_cdim);
        dim3 numBlocks_cdim(num_blocks_cdim);

        int *d_dims_order_start;
        cudaMalloc((void**)&d_dims_order_start, cdim * sizeof(int));
        fill<<<numBlocks_cdim, threadsPerBlock_cdim>>>(d_dims_order_start, -1, cdim);
        get_dims_order_start<<<numBlocks_m, threadsPerBlock_m>>>(d_dims, d_dims_order, d_dims_order_start, m);
        cudaDeviceSynchronize();

        //DEBUG
        if (1 == 0){
            printvec(d_dims, m, "d_dims");
            printvec(d_dims_order, m, "d_dims_order");
            printvec(d_dims_order_next, m, "d_dims_order_next");
            printvec(d_dims_order_start, cdim, "d_dims_order_start");
            printf("complex_dim = %d\n", complex_dim);
        }

        // beta
        int *d_beta, *d_left;
        cudaMalloc((void**)&d_beta, m * sizeof(int));
        cudaMalloc((void**)&d_left, m * sizeof(int));

        create_beta(d_beta, d_left, h_rows, h_cols, m, nnz);

        //DEBUG
        if (1 == 0){
            printvec(d_low, m, "d_low");
            printvec(d_beta, m, "d_beta");
            printvec(d_left, m, "d_left");
            for (int i = 0; i < nnz; i++){
                printf("(%d, %d) ", h_rows[i], h_cols[i]);
            }
        }

        // -------------------------------
        // Algorithms
        // -------------------------------

        int iter  = 0;
        if (strcmp(algstr, "std")==0){
            standard(d_rows_mp, d_aux_mp, d_low, d_arglow, m, p, resRecord, timeRecord, &iter, numBlocks_m, threadsPerBlock_m);
        } else if (strcmp(algstr, "twist")==0){
            twist(d_rows_mp, d_aux_mp, d_low, d_arglow, d_dims, d_dims_order, d_dims_order_next, d_dims_order_start, complex_dim, m, p, resRecord, timeRecord, &iter, numBlocks_m, threadsPerBlock_m);
        } else if (strcmp(algstr, "ph_row")==0){
            ph_row(d_rows_mp, d_aux_mp, d_low, d_arglow, m, p, resRecord, timeRecord, &iter, numBlocks_m, threadsPerBlock_m);
        } else if (strcmp(algstr, "pms")==0){
            pms(d_rows_mp, d_aux_mp, d_low, d_arglow, d_dims, d_dims_order, d_dims_order_next, d_dims_order_start, m, p, complex_dim, d_left, d_beta, resRecord, timeRecord, &iter, numBlocks_m, threadsPerBlock_m, numBlocks_cdim, threadsPerBlock_cdim);
        }else
            printf("Not recognised");

        // matrix: device to host
        indexShiftUp<<<numBlocks_m,threadsPerBlock_m>>>(d_low, m); 
        cudaMemcpy(h_low,  d_low,  sizeof(int)*m, cudaMemcpyDeviceToHost);

        cudaThreadSynchronize();
        cudaEventRecord(stopIHT,0);
        cudaEventSynchronize(stopIHT);
        cudaEventElapsedTime(&timeIHT, startIHT, stopIHT);
        cudaEventDestroy(startIHT);
        cudaEventDestroy(stopIHT);

        // CLEANUP
        // free up the allocated memory on the device

        // clear beta

        cudaFree(d_low);
        cudaFree(d_arglow);

        cudaFree(d_beta);
        cudaFree(d_left);

        cudaFree(d_cols);
        cudaFree(d_rows);
        cudaFree(d_rows_mp);
        cudaFree(d_aux_mp);

        cudaFree(d_dims);
        cudaFree(d_dims_order);
        cudaFree(d_dims_order_next);
        cudaFree(d_dims_order_start);

        cudaDeviceSynchronize();
        cublasShutdown();

    }  //closes the else ensuring a correct number of input and output arguments

    return;
}


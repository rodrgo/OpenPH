
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
        printf("is_col_order = %d\n", is_col_order);

        // Count nonzeros in columns and rows
        int max_nnz_rows = get_max_nnz(h_rows, m, nnz);
        int max_nnz_cols = get_max_nnz(h_cols, m, nnz);

        printf("max_nnz_rows : %d\n", max_nnz_rows);
        printf("max_nnz_cols : %d\n", max_nnz_cols);
          
        int p = 3*max_nnz_cols;
        int mp = m * p;

        printf("p = %d\n", p);

        // -------------------------------
        // GPU
        // -------------------------------

        int gpu_number = 3;
        int threads_perblock_m  = 0;
        int threads_perblock_nnz = 0;
        int threads_perblock_mp = 0;
        set_gpu_device(gpu_number, &threads_perblock_m, &threads_perblock_nnz, &threads_perblock_mp, m, nnz, p);

        printf("threads_perblock_m:   %d\n", threads_perblock_m);
        printf("threads_perblock_nnz: %d\n", threads_perblock_nnz);
        printf("threads_perblock_mp:  %d\n", threads_perblock_mp);

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
        // d_low, d_arglow
        int *d_low, *d_arglow;

        cudaMalloc((void**)&d_low,    m * sizeof(int));
        cudaMalloc((void**)&d_arglow, m * sizeof(int));

        minusone_vector_int<<< numBlocks_m, threadsPerBlock_m >>>((int*)d_low, m);
        minusone_vector_int<<< numBlocks_m, threadsPerBlock_m >>>((int*)d_arglow, m);
        cudaDeviceSynchronize();

        // d_aux_mp, d_rows_mp
        int *d_rows_mp, *d_aux_mp;

        cudaMalloc((void**)&d_aux_mp,   mp * sizeof(int));
        cudaMalloc((void**)&d_rows_mp,  mp * sizeof(int));

        minusone_vector_int<<< numBlocks_mp, threadsPerBlock_mp >>>((int*)d_aux_mp, mp);
        minusone_vector_int<<< numBlocks_mp, threadsPerBlock_mp >>>((int*)d_rows_mp, mp);
        create_rows_mp<<< numBlocks_nnz, threadsPerBlock_nnz>>>(d_rows, d_cols, d_rows_mp, m, p, nnz);
        cudaDeviceSynchronize();

        compute_low_mp<<<numBlocks_m,threadsPerBlock_m>>>(d_rows_mp, d_low, m, p);
        cudaDeviceSynchronize();

        //DEBUG
        if (1 == 0){
            printf("ARGLOW\n");
            printvec(d_arglow, m);
            printf("LOW\n");
            printvec(d_low, m);
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
        // beta
        // -------------------------------


        int iter  = 0;
        if (strcmp(algstr, "std")==0){
            standard(d_rows_mp, d_aux_mp, d_low, d_arglow, m, p, resRecord, timeRecord, &iter, numBlocks_m, threadsPerBlock_m);
        } else if (strcmp(algstr, "twist")==0){
            twist(d_rows_mp, d_aux_mp, d_low, d_arglow, m, p, resRecord, timeRecord, &iter, numBlocks_m, threadsPerBlock_m);
        } else if (strcmp(algstr, "pms")==0){

            // beta
            int *d_beta;
            cudaMalloc((void**)&d_beta, m * sizeof(int));
            create_beta(d_beta, h_rows, h_cols, m, nnz);

            // simplex dimensions

            pms(d_rows_mp, d_aux_mp, d_low, d_arglow, m, p, d_beta, resRecord, timeRecord, &iter, numBlocks_m, threadsPerBlock_m);
        }else
            printf("Not recognised");

        // matrix: device to host
        indexShiftUp<<<numBlocks_nnz,threadsPerBlock_nnz>>>(d_rows, nnz);
        indexShiftUp<<<numBlocks_nnz,threadsPerBlock_nnz>>>(d_cols, nnz); 

        cudaMemcpy(h_rows, d_rows, sizeof(int)*m, cudaMemcpyDeviceToHost);
        cudaMemcpy(h_cols, d_cols, sizeof(int)*m, cudaMemcpyDeviceToHost);
        cudaMemcpy(h_low,  d_low,  sizeof(int)*m, cudaMemcpyDeviceToHost);

        cudaThreadSynchronize();
        cudaEventRecord(stopIHT,0);
        cudaEventSynchronize(stopIHT);
        cudaEventElapsedTime(&timeIHT, startIHT, stopIHT);
        cudaEventDestroy(startIHT);
        cudaEventDestroy(stopIHT);

        // CLEANUP
        // free up the allocated memory on the device

        cudaFree(d_beta);
        cudaFree(d_rows);
        cudaFree(d_cols);
        
        cudaFree(d_rows_mp);
        cudaFree(d_aux_mp);

        cublasShutdown();

    }  //closes the else ensuring a correct number of input and output arguments

    return;
}


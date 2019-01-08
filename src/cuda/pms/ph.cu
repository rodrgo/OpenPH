
#include "init.cu"

// PMS
// Pass boundary matrix in row-column format 
// Convert to GPU
// Compute low and lowstar in GPU

// Host function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

    if ( nlhs!=3 & nrhs!=6 )
        printf("[Error in number of inputs or number of outputs]");
    else {

        // -------------------------------
        // Read Inputs
        // -------------------------------

        // inputs
        int *h_rows_in, *h_cols_in;
        int m;
        int col_width;
        int *h_low_true;

        // reading in the string to determine the algorithm
        int strlen = mxGetN(prhs[0])+1;
        char algstr[strlen+100];
        int algerr = mxGetString(prhs[0], algstr, strlen);

        // Read data
        h_rows_in  = (int*)mxGetData(prhs[1]);
        h_cols_in  = (int*)mxGetData(prhs[2]);
        m          = (int)mxGetScalar(prhs[3]);
        col_width  = (int)mxGetScalar(prhs[4]);  
        h_low_true = (int*)mxGetData(prhs[5]);  

        // -------------------------------
        // Get nnz and p
        // -------------------------------

        int nnz = mxGetM(prhs[1]);
        int p   = col_width * max_nnz(h_cols_in, m, nnz);
        int mp  = m * p;

        if (assert_col_order(h_cols_in, nnz) == 0){
            printf("WARNING: Matrix incorrect!\n");
        }

        // -------------------------------
        // GPU
        // -------------------------------

        int gpu_number = 3;

        int tpb_m   = 0; // threads per block (m)
        int tpb_nnz = 0; // threads per block (nnz)
        int tpb_mp  = 0; // threads per block (mp)

        set_gpu_device(gpu_number, &tpb_m, &tpb_nnz, &tpb_mp, m, nnz, p);

        dim3 TPBnnz(tpb_nnz);
        dim3 NBnnz(num_blocks(nnz, tpb_nnz));

        dim3 TPBm(tpb_m);
        dim3 NBm(num_blocks(m, tpb_m));

        dim3 TPBmp(tpb_mp);
        dim3 NBmp(num_blocks(mp, tpb_mp));

        // -------------------------------
        // Outputs
        // -------------------------------

        // outputs
        int *h_low;
        float *error_linf, *error_lone, *error_redu, *error_ess;
        float *time_track; 
        int *h_ess;
        int *num_iters;
        
        plhs[0] = mxCreateNumericMatrix(1, m, mxINT32_CLASS, mxREAL);
        plhs[1] = mxCreateNumericMatrix(1, m, mxINT32_CLASS, mxREAL);
        plhs[2] = mxCreateNumericMatrix(1, m+1, mxSINGLE_CLASS, mxREAL);
        plhs[3] = mxCreateNumericMatrix(1, m+1, mxSINGLE_CLASS, mxREAL);
        plhs[4] = mxCreateNumericMatrix(1, m+1, mxSINGLE_CLASS, mxREAL);
        plhs[5] = mxCreateNumericMatrix(1, m+1, mxSINGLE_CLASS, mxREAL);
        plhs[6] = mxCreateNumericMatrix(1, m+1, mxSINGLE_CLASS, mxREAL);
        plhs[7] = mxCreateNumericMatrix(1, 1  , mxINT32_CLASS, mxREAL); 

        h_low      = (int*) mxGetData(plhs[0]);  
        h_ess      = (int*) mxGetData(plhs[1]);
        error_linf = (float*) mxGetData(plhs[2]);
        error_lone = (float*) mxGetData(plhs[3]);
        error_redu = (float*) mxGetData(plhs[4]);
        error_ess  = (float*) mxGetData(plhs[5]); 
        time_track = (float*) mxGetData(plhs[6]);
        num_iters  = (int*) mxGetData(plhs[7]);

        // Initialise

        for (int i=0; i<m  ; i++) h_low[i]      = -1;
        for (int i=0; i<m  ; i++) h_ess[i]      = -1;
        for (int i=0; i<m+1; i++) error_linf[i] = -1;
        for (int i=0; i<m+1; i++) error_lone[i] = -1;
        for (int i=0; i<m+1; i++) error_redu[i] = -1;
        for (int i=0; i<m+1; i++) error_ess[i]  = -1;
        for (int i=0; i<m+1; i++) time_track[i] = 0;

        // -------------------------------
        // Create data on device
        // -------------------------------

        // d_rows, d_cols
        int *d_rows; 
        cudaMalloc((void**)&d_rows, mp * sizeof(int));
        create_rows(h_rows_in, h_cols_in, d_rows, m, p, nnz, NBnnz, TPBnnz, NBmp, TPBmp);

        // device vectors
        int *d_low, *d_arglow, *d_classes, *d_ess;
        int *d_aux_mp;

        cudaMalloc((void**)&d_low, m * sizeof(int));
        cudaMalloc((void**)&d_arglow, m * sizeof(int));
        cudaMalloc((void**)&d_classes, m * sizeof(int));
        cudaMalloc((void**)&d_ess, m * sizeof(int));
        cudaMalloc((void**)&d_aux_mp, mp * sizeof(int));

        fill<<<NBm, TPBm>>>(d_low, -1, m);
        fill<<<NBm, TPBm>>>(d_arglow, -1, m);
        fill<<<NBm, TPBm>>>(d_classes, 0, m);
        fill<<<NBm, TPBm>>>(d_ess, 1, m);
        fill<<<NBmp,TPBmp>>>(d_aux_mp, -1, mp);
        cudaDeviceSynchronize();

        compute_low<<<NBm, TPBm>>>(d_rows, d_low, m, p);
        cudaDeviceSynchronize();

        // d_float_m
        float *d_float_m;
        cudaMalloc((void**)&d_float_m, m * sizeof(float));

        // d_low_true
        int *d_low_true;
        cudaMalloc((void**)&d_low_true, m * sizeof(int));
        cudaMemcpy(d_low_true, h_low_true, sizeof(int)*m, cudaMemcpyHostToDevice);
        cudaDeviceSynchronize();

        // ... compute norm of d_low_true
        to_float<<<NBm, TPBm>>>(d_float_m, d_low_true, m);
        cudaDeviceSynchronize();

        float norm1_low_true = norm_1(d_float_m, m); 
        float norminf_low_true = norm_inf(d_float_m, m); 

        // ... Now shift index down
        indexShiftDown<<<NBm, TPBm>>>(d_low_true, m);
        cudaDeviceSynchronize();

        // d_ess_true
        int *d_ess_true;
        cudaMalloc((void**)&d_ess_true, m * sizeof(int));

        fill<<<NBm, TPBm>>>(d_ess_true, 1, m);
        cudaDeviceSynchronize();

        compute_ess_true<<<NBm, TPBm>>>(d_low_true, d_ess_true, m);
        cudaDeviceSynchronize();

        // -------------------------------
        // Get PH vectors
        // -------------------------------

        float time;
        cudaEvent_t start, stop;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        // -------------------------------
        // Dimensions
        // -------------------------------
        
        // d_dim:           Dim of simplex j (-1, 0, 1, ..., complex_dim)
        // d_dim_order[j]:  Order of simplex j in dimension d_dim[j]
        // d_dim_start[d]:  Position of first simplex in dimension "d" in d_dim
        // d_dim_next[j]:   Simplex of dimension d_dim[j] after "j"

        int *d_dim, *d_dim_order, *d_dim_next, *d_dim_start;
        int complex_dim = -1;

        cudaMalloc((void**)&d_dim, m * sizeof(int));
        cudaMalloc((void**)&d_dim_order, m * sizeof(int));
        cudaMalloc((void**)&d_dim_next, m * sizeof(int));

        compute_simplex_dimensions_h(h_cols_in, m, p, nnz, d_dim, d_dim_order, d_dim_next, &complex_dim);

        int cdim = complex_dim + 2;
        cudaMalloc((void**)&d_dim_start, cdim * sizeof(int));

        int threads_perblock_cdim = min(cdim, tpb_m);
        dim3 TPBcdim(threads_perblock_cdim);
        dim3 NBcdim(num_blocks(cdim, threads_perblock_cdim));

        create_dim_start(d_dim, d_dim_order, d_dim_start, cdim, m, NBm, TPBm, NBcdim, TPBcdim);

        // left, beta
        int *d_beta, *d_left;
        cudaMalloc((void**)&d_beta, m * sizeof(int));
        cudaMalloc((void**)&d_left, m * sizeof(int));
        create_beta(d_beta, d_left, h_rows_in, h_cols_in, m, nnz);

        // -------------------------------
        // Algorithms
        // -------------------------------

        int iter  = 0;
        algorithm_factory(algstr, d_low, d_arglow,
                d_classes, d_ess, d_rows, d_dim,
                d_dim_order, d_dim_next, d_dim_start, 
                d_beta, d_left,
                m, p, complex_dim, d_aux_mp, d_low_true, d_ess_true, d_float_m,
                error_lone, error_linf, error_redu, error_ess,
                time_track, &iter, NBm, TPBm, NBcdim,
                TPBcdim);

        // Record iters to output
        num_iters[0] = iter;

        // scale remaining trackers
        for(int i=0; i<m+1; i++)
            error_lone[i] = error_lone[i]/norm1_low_true;

        for(int i=0; i<m+1; i++)
            error_linf[i] = error_linf[i]/norminf_low_true;

        // matrix: device to host
        indexShiftUp<<<NBm, TPBm>>>(d_low, m); 
        cudaDeviceSynchronize();

        cudaMemcpy(h_low, d_low, sizeof(int)*m, cudaMemcpyDeviceToHost);
        cudaMemcpy(h_ess, d_ess, sizeof(int)*m, cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();

        cudaThreadSynchronize();
        cudaEventRecord(stop,0);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&time, start, stop);
        cudaEventDestroy(start);
        cudaEventDestroy(stop);

        // cudaFree

        cudaFree(d_low);
        cudaFree(d_arglow);
        cudaFree(d_classes);
        cudaFree(d_ess);

        cudaFree(d_low_true);
        cudaFree(d_ess_true);

        cudaFree(d_beta);
        cudaFree(d_left);

        cudaFree(d_rows);
        cudaFree(d_aux_mp);

        cudaFree(d_float_m);

        cudaFree(d_dim);
        cudaFree(d_dim_order);
        cudaFree(d_dim_next);
        cudaFree(d_dim_start);

        cudaDeviceSynchronize();
        cublasShutdown();

    }  //closes the else ensuring a correct number of input and output arguments

    return;
}


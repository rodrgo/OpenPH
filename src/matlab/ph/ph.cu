
#include <mex.h>
#include "../../cuda/openph.cu"

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
        // Inputs
        // -------------------------------

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

        int nnz = mxGetM(prhs[1]);

        // -------------------------------
        // Outputs
        // -------------------------------

        int *h_low;
        int *h_ess;
        float *error_linf, *error_lone, *error_redu, *error_ess;
        float *time_track; 
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
        // OpenPH
        // -------------------------------

        openph(algstr, h_rows_in, h_cols_in, m, col_width,
                h_low_true, nnz, h_low, h_ess, error_linf,
                error_lone, error_redu, error_ess, time_track,
                num_iters);

    }  //closes the else ensuring a correct number of input and output arguments

    return;
}


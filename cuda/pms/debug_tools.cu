
inline void h2d(int *d_v, int j, int v){
    cudaMemcpy(d_v+j, &v, sizeof(int), cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
    return;
}

int d2h(int *d_v, int j){
    int v_j = 0;
    cudaMemcpy(&v_j, d_v+j, sizeof(int), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    return v_j;
}

inline void printvec_float(float *d_v, int m, char* s){
    float val = 0;
    printf("%s\n",s);
    for (int i = 0; i < m; i++){
        cudaMemcpy(&val, d_v+i, sizeof(float), cudaMemcpyDeviceToHost);
        printf("[%d: %5.6f], ", i, val);
    }
    printf("\n\n");
}

inline void printvec(int *d_v, int m, char* s){
    int val = 0;
    printf("%s\n",s);
    for (int i = 0; i < m; i++){
        cudaMemcpy(&val, d_v+i, sizeof(int), cudaMemcpyDeviceToHost);
        printf("[%d: %d], ", i, val);
    }
    printf("\n\n");
}

inline void print_matrix_cols_mp(int *d_rows_mp, int m, int p){
    int val = 0;
    for (int j = 0; j < m; j++){
        printf("col (%d): ", j);
        for (int i = 0; i < p; i++){
            cudaMemcpy(&val, d_rows_mp+(j*p+i), sizeof(int), cudaMemcpyDeviceToHost);
            printf("%d, ", val);
        } 
        printf("\n");
    }
}

inline void print_matrix_cols(int *d_rows, int *d_cols, int nnz){
    int col = -1;
    int val = 0;
    for (int idx = 0; idx < nnz; idx++){
        cudaMemcpy(&val, d_cols+idx, sizeof(int), cudaMemcpyDeviceToHost);
        if (val != col){
            printf("\n");
            printf("col (%d): ", val);
            col = val;
        }
        cudaMemcpy(&val, d_rows+idx, sizeof(int), cudaMemcpyDeviceToHost);
        printf("%d, ", val);
        cudaMemcpy(&val, d_cols+idx, sizeof(int), cudaMemcpyDeviceToHost);
    } 
    printf("\n");
}

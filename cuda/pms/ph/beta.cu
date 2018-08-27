
inline void create_beta_h(int *h_beta, int *h_left, int *h_rows, int *h_cols, int m, int nnz){
    // Note: Index in h_rows and h_cols starts at 1
    int *h_visited = (int*)malloc( sizeof(int) * m );
    for(int i=0; i<m; i++) h_left[i] = -1;
    for(int i=0; i<m; i++) h_visited[i] = 0;
    for(int i=0; i<m; i++) h_beta[i] = -1;
    for(int l=0; l<nnz; l++)
        if (h_visited[h_rows[l]-1] == 0){
            h_beta[h_cols[l]-1] = h_beta[h_cols[l]-1] > h_rows[l]-1 ? h_beta[h_cols[l]-1] : h_rows[l]-1;
            h_visited[h_rows[l]-1] = 1;
            h_left[h_rows[l]-1] = h_cols[l]-1;
        }
    free(h_visited);
}

inline void create_beta(int *d_beta, int *d_left, int *h_rows, int *h_cols, int m, int nnz){
    int *h_beta = (int*)malloc( sizeof(int) * m );
    int *h_left = (int*)malloc( sizeof(int) * m );
    create_beta_h(h_beta, h_left, h_rows, h_cols, m, nnz);
    cudaMemcpy(d_beta, h_beta, m*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_left, h_left, m*sizeof(int), cudaMemcpyHostToDevice);
    cudaDeviceSynchronize();
    free(h_beta);
    free(h_left);
}


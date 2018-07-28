inline void standard(int *d_rows_mp, int *d_aux_mp, int *d_low, int *d_arglow, const int m, const int p, float *resRecord, float *timeRecord, int *p_iter, dim3 numBlocks_m, dim3 threadsPerBlock_m){
    for(int j = 0; j < m; j++){
        reduce_col_gpu<<<numBlocks_m, threadsPerBlock_m>>>(j, d_rows_mp, d_aux_mp, d_low, d_arglow, m, p);
    } 
}

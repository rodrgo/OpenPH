inline void standard(int *d_rows_mp, int *d_aux_mp, int *d_low, int *d_arglow, const int m, const int p, float *resRecord, float *timeRecord, int *p_iter, dim3 NBm, dim3 TPBm){
    for(int j = 0; j < m; j++){
        reduce_col<<<NBm, TPBm>>>(j, d_rows_mp, d_aux_mp, d_low, d_arglow, m, p);
    }
}

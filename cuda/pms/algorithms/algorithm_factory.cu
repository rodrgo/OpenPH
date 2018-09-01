inline void algorithm_factory(char *algstr,
        int *d_low, int *d_arglow, int *d_classes, int *d_ess,
        int *d_rows, int *d_dim, int *d_dim_order, int *d_dim_next, int *d_dim_start,
        int *d_beta, int *d_left,
        int m, int p, int complex_dim, 
        int *d_aux_mp, int *d_low_true, int *d_ess_true,
        float *error_lone, float *error_linf, float *error_redu, float *error_ess, float *time_track, int *p_iter,
        dim3 NBm, dim3 TPBm, dim3 NBcdim, dim3 TPBcdim){

    if (strcmp(algstr, "standard")==0){
        standard(d_low, d_arglow, d_classes, d_ess, 
                d_rows, m, p, d_aux_mp, d_low_true,
                d_ess_true,
                error_lone, error_linf, error_redu, 
                error_ess, time_track, p_iter, NBm, 
                TPBm); 
    } else if (strcmp(algstr, "twist")==0){
        twist(d_rows, d_aux_mp, d_low, d_arglow, d_dim, d_dim_order, d_dim_next, d_dim_start, complex_dim, m, p, error_linf, time_track, p_iter, NBm, TPBm);
    } else if (strcmp(algstr, "ph_row")==0){
        ph_row(d_rows, d_aux_mp, d_low, d_arglow, m, p, error_linf, time_track, p_iter, NBm, TPBm);
    } else if (strcmp(algstr, "pms")==0){
        pms(d_rows, d_aux_mp, d_low, d_arglow, d_dim, d_dim_order, d_dim_next, d_dim_start, m, p, complex_dim, d_left, d_beta, error_linf, time_track, p_iter, NBm, TPBm, NBcdim, TPBcdim);
    }else
        printf("Not recognised");

}

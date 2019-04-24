
inline void tic(cudaEvent_t *p_start, cudaEvent_t *p_stop){
    //cudaEvent_t start, stop;
    cudaEventCreate(p_start);
    cudaEventCreate(p_stop);
    cudaEventRecord(p_start[0], 0);
}

inline void toc(cudaEvent_t start, cudaEvent_t stop, float *p_time){
    cudaThreadSynchronize();
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(p_time, start, stop);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
}

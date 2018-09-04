% DO_ALL_CUDA.M
% Perform all analyses and create all plots for paper

tic;
benchmark_pms
fileID = fopen('do_all_log.txt','a');
fprintf(fileID,'benchmark_pms: %6.6f secs\n', toc);
fclose(fileID);


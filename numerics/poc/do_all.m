% DO_ALL.M
% Perform all analyses and create all plots for paper


tic;
test_vr_ensembles
fileID = fopen('do_all_log.txt','w');
fprintf(fileID,'test_vr_ensembles: %6.6f secs\n', toc);
fclose(fileID);

tic;
test_vr_ensembles_by_numpoints
fileID = fopen('do_all_log.txt','a');
fprintf(fileID,'test_vr_ensembles_by_numpoints: %6.6f secs\n', toc);
fclose(fileID);

tic;
test_speed_standard_reduction
fileID = fopen('do_all_log.txt','a');
fprintf(fileID,'test_speed_standard_reduction: %6.6f secs\n', toc);
fclose(fileID);

tic;
compare_std_vs_twist
fileID = fopen('do_all_log.txt','a');
fprintf(fileID,'compare_std_vs_twist: %6.6f secs\n', toc);
fclose(fileID);

tic;
plot_morozov_matrix
fileID = fopen('do_all_log.txt','a');
fprintf(fileID,'plot_morozov_matrix: %6.6f secs\n', toc);
fclose(fileID);

tic;
test_morozov_cubic_time
fileID = fopen('do_all_log.txt','a');
fprintf(fileID,'test_morozov_cubic_time: %6.6f secs\n', toc);
fclose(fileID);

tic;
alpha_beta_curves
fileID = fopen('do_all_log.txt','a');
fprintf(fileID,'alpha_beta_curves: %6.6f secs\n', toc);
fclose(fileID);

tic;
rho_curve_example
fileID = fopen('do_all_log.txt','a');
fprintf(fileID,'rho_curve_example: %6.6f secs\n', toc);
fclose(fileID);

tic;
test_essential_reduction
fileID = fopen('do_all_log.txt','a');
fprintf(fileID,'test_essential_reduction: %6.6f secs\n', toc);
fclose(fileID);

tic;
test_curiosity_8_red
fileID = fopen('do_all_log.txt','a');
fprintf(fileID,'test_curiosity_8_red: %6.6f secs\n', toc);
fclose(fileID);

tic;
benchmark_alpha_beta_parallel
fileID = fopen('do_all_log.txt','a');
fprintf(fileID,'benchmark_alpha_beta_parallel: %6.6f secs\n', toc);
fclose(fileID);

tic;
average_percentage_unreduced
fileID = fopen('do_all_log.txt','a');
fprintf(fileID,'average_percentage_unreduced: %6.6f secs\n', toc);
fclose(fileID);

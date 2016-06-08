% DO_ALL.M
% Perform all analyses and create all plots for paper


tic;
explore_vr_ensemble_complexity
fileID = fopen('do_all_log.txt','w');
fprintf(fileID,'explore_vr_ensemble_complexity: %6.6f secs\n', toc);
fclose(fileID);

tic;
test_speed_standard_reduction
fileID = fopen('do_all_log.txt','a');
fprintf(fileID,'test_speed_standard_reduction: %6.6f secs\n', toc);
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


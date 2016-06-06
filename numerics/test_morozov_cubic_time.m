% Morozov takes cubic time
% Note: make sure that you give matlab enough heap space to work with
init;

figure_dir = './figures/';

%% Test Morozov cubic time

% Matrix for order 11 is too heavy for dense reduction
orders = 3:11;

stream_sizes = zeros(size(orders));
reduce_times_sparse = zeros(size(orders));
%reduce_times_dense = zeros(size(orders));

for i = 1:length(orders)
    order = orders(i);
    fprintf('Working on order %d...\t', order);
    
    time_order = tic;
    stream = examples.MorozovCubicTimeExample.getMorozovCubicTimeExample(order);
    stream_sizes(i) = stream.getSize();
    D = BoundaryMatrix(stream);
    
    time_sparse = tic;
    [~, ~] = D.do_classical_reduction_optimised();
    reduce_times_sparse(i) = toc(time_sparse);
    
    %time_dense = tic;
    %[~, ~] = D.do_classical_reduction_optimised_dense();
    %reduce_times_dense(i) = toc(time_dense);
    
    fprintf('done in %g sec!\n', toc(time_order));
end

%%

n_cube = stream_sizes.^3;

h_sparse = loglog(stream_sizes, reduce_times_sparse, 'r--');
set(gcf, 'color', [1 1 1]);
xlabel('Number of simplices (N)');
ylabel('Time (sec)');
title('Running time for Morozov complex');

hold on;
%h_dense = loglog(stream_sizes, reduce_times_dense, 'b--');
stream_ref = 2*stream_sizes(5:end);
h_ref = loglog(stream_ref, 1e1*(stream_ref/max(stream_ref)).^3, '-');
hold off;
legend([h_sparse, h_ref], 'Classic reduction algorithm', 't = C*N^3');

fileName = strcat(figure_dir, 'morozov_cubic_time');
print('-depsc', strcat(fileName, '.eps'));



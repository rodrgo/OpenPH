% Morozov takes cubic time
% Note: make sure that you give matlab enough heap space to work with
init;

import edu.stanford.math.plex4.*;

figure_dir = './figures/';
figure_tag = 'speed_morozov';

% Plot params
LW = 'LineWidth';
MS = 'MarkerSize';
markers = '+o*.xsd^v><ph';

% Algorithms to test
algorithms = {'std', 'twist', 'alpha_beta_std', 'rho_std'};

% Matrix dense?
as_dense = true;

% Make labels for plotting`
algorithms_labels = algorithms;
for i = 1:length(algorithms)
    algorithms_labels{i} = strrep(algorithms{i}, '_', '\_');
end

% Matrix for order 11 is too heavy for dense reduction
orders = 3:12;
stream_sizes = zeros(size(orders));
time_algorithms = zeros(length(algorithms), length(orders));

for i = 1:length(orders)
    order = orders(i);
    fprintf('Working on order %d...\t', order);
    
    time_order = tic;
    stream = examples.MorozovCubicTimeExample.getMorozovCubicTimeExample(order);
    stream_sizes(i) = stream.getSize();

    for l = 1:length(algorithms)
        algorithm = algorithms{l};
        [~, t, D] = reduce_stream(stream, 'unreduced', algorithm, as_dense);
        time_algorithms(l, i) = 1000*t;
    end

    fprintf('done in %g sec!\n', toc(time_order));
end

%%
% Plot results
handles = [];
set(gcf, 'color', [1 1 1]);
set(gca, 'Fontname', 'setTimes', 'Fontsize', 15);

for ii = 1:size(time_algorithms, 1)
    x = stream_sizes;
    y = time_algorithms(ii, :);
    handles(end + 1) = loglog(x, y, [markers(ii) '-'], LW, 1.5, MS, 10);
    hold on;
end

%stream_ref = 2*stream_sizes(5:end);
%handles(end + 1) = loglog(stream_ref, 1e1*(stream_ref/max(stream_ref)).^3, '-');
stream_ref = stream_sizes(8:10);
C = 1e-10;
handles(end + 1) = loglog(stream_ref, C*(stream_ref).^3, '-');
algorithms_labels{end + 1} = 't = C*m^3'; 

xlabel('m');
ylabel('time (ms)');
legend(handles, algorithms_labels, 'Location', 'NorthWest');
title('Running time for Morozov complex');
hold off;

file_name = strcat(figure_tag, '.eps');
file_path = strcat(figure_dir, file_name);
print('-depsc', file_path);
eps_to_pdf(file_path);



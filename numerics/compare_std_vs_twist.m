% test_speed_standard_reduction.m
init;

figure_dir = './figures/';

% Plot params
LW = 'LineWidth';
MS = 'MarkerSize';
markers = '+o*.xsd^v><ph';

% Explore complexity of Vietoris-Rips complexes
% We do not include random_torus.
vr_complexes = {'random_figure_8', ...
    'random_gaussian', 'random_trefoil_knot'};

% Algorithms to test
algo_cores = {'alpha_beta', 'rho', 'c8'};

% Build algorithms
algorithms = {'std_dense', 'twist_dense'};
for l = 1:length(algo_cores)
    algorithms{end + 1} = [algo_cores{l} '_' 'std_dense'];
    algorithms{end + 1} = [algo_cores{l} '_' 'twist_dense'];
end

% Make labels for plotting`
algorithms_labels = algorithms;
for i = 1:length(algorithms)
    algorithms_labels{i} = strrep(algorithms{i}, '_', '\_');
end

% Fixed complex parameters
max_dimension = 200;
num_divisions = 20;

% Variable complex parameters
max_filtration_values = 1:7;

% Number of complexes per parameter
num_samples = 10;

% Type of plots
plot_types = {'m', 'nnz'};

for i = 1:length(vr_complexes)

    complex = vr_complexes{i};

    % D.m, nnz(D), time in ms, algorithm
    time_data = [];

    for j = 1:length(max_filtration_values)

        mfv = max_filtration_values(j);
        fprintf('%s max_filtr_value: %s\n', complex, num2str(mfv));

        for k = 1:num_samples

            fprintf('\tSample %d/%d\n', k, num_samples);
            stream = example_factory(complex, max_dimension, mfv, num_divisions);

            for l = 1:length(algorithms)
                algorithm = algorithms{l};
                D = BoundaryMatrix(stream, 'unreduced');
                [lows, t] = reduce_matrix(D, algorithm);
                row = [D.m, nnz(D.matrix), 1000*t, l];
                time_data = [time_data; row];
            end

        end

    end

    % Relative speed results
    figure_tag = 'rel_speed_std_vs_twist';
    for p = 1:length(plot_types)
        x_var = plot_types{p};
        handles = [];
        set(gcf, 'color', [1 1 1]);
        set(gca, 'Fontname', 'setTimes', 'Fontsize', 15);
        rel_labels = {};

        % We want to beat Twist algorithm, so compare everything to it.
        for l = 1:floor(length(algorithms)/2)
            algo_num_idx = time_data(:, end) == (2*l-1);
            algo_den_idx = time_data(:, end) == (2*l);
            x = time_data(algo_num_idx, p);
            y = time_data(algo_num_idx, 3)./time_data(algo_den_idx, 3);
            handles(end + 1) = plot(x, y, 'x');
            rel_labels{end + 1} = [algorithms_labels{2*l-1}, '/', algorithms_labels{2*l}];
            hold on;
        end
        plot(xlim, ones(size(xlim)), '--k');
        hold on;

        xlabel(x_var);
        ylabel('relative time');
        legend(handles, rel_labels, 'Location', 'NorthEast');
        params_tag = sprintf('max_dim = %d, num_div = %d', max_dimension, num_divisions);
        params_tag = strrep(params_tag, '_', '\_');
        title_str = [strrep(complex, '_', '\_'), ', ', params_tag];
        title(title_str);
        hold off;

        file_name = [complex, '_', x_var, '_', figure_tag, '.eps'];
        file_path = strcat(figure_dir, file_name);
        print('-depsc', file_path);
        eps_to_pdf(file_path);
    end

end


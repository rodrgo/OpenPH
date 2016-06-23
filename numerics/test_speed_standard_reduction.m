% test_speed_standard_reduction.m
init;

figure_dir = './figures/';
figure_tag = 'speed_algos';

% Plot params
LW = 'LineWidth';
MS = 'MarkerSize';
markers = '+o*.xsd^v><ph';

% Explore complexity of Vietoris-Rips complexes
vr_complexes = {'random_torus', 'random_figure_8', ...
    'random_gaussian', 'random_trefoil_knot'};

% Algorithms to test
algorithms = {'std_sparse', 'std_dense', ...
    'twist_sparse', 'twist_dense', ... 
    'alpha_beta_sparse', 'alpha_beta_dense', ...
    'rho_sparse', 'rho_dense'};

% Make labels for plotting`
algorithms_labels = algorithms;
for i = 1:length(algorithms)
    algorithms_labels{i} = strrep(algorithms{i}, '_', '\_');
end

% Fixed complex parameters
max_dimension = 5;
num_divisions = 12;

% Variable complex parameters
max_filtration_values = 1:7;

% Number of complexes per parameter
num_samples = 10;

for i = 1:length(vr_complexes)

    complex = vr_complexes{i};
    avg_m = zeros(1, length(max_filtration_values));
    times_complex = zeros(length(algorithms), length(avg_m));

    for j = 1:length(max_filtration_values)

        mfv = max_filtration_values(j);
        fprintf('%s max_filtr_value: %s\n', complex, num2str(mfv));
        time_algorithms = zeros(length(algorithms), num_samples);
        ms = zeros(1, num_samples);

        for k = 1:num_samples

            fprintf('\tSample %d/%d\n', k, num_samples);
            stream = example_factory(complex, max_dimension, mfv, num_divisions);

            for l = 1:length(algorithms)
                algorithm = algorithms{l};
                D = BoundaryMatrix(stream, 'unreduced');
                [lows, t] = reduce_matrix(D, algorithm);
                time_algorithms(l, k) = 1000*t;
            end

            ms(k) = D.m;

        end

        avg_m(j) = round(mean(ms));
        times_complex(:, j) = mean(time_algorithms, 2);

    end

    % Plot results
    handles = [];
    set(gcf, 'color', [1 1 1]);
    set(gca, 'Fontname', 'setTimes', 'Fontsize', 15);

    for ii = 1:size(times_complex, 1)
        x = avg_m;
        y = times_complex(ii, :);
        handles(end + 1) = semilogy(x, y, [markers(ii) '-'], LW, 1.5, MS, 10);
        hold on;
    end

    xlabel('Average m');
    ylabel('time ms');
    legend(handles, algorithms_labels, 'Location', 'NorthWest');
    params_tag = sprintf('(max\_dim, num\_div, mfv)=(%d,%d,%d)', max_dimension, num_divisions, mfv);
    title_str = [strrep(complex, '_', '\_') ', ' params_tag];
    title(title_str);
    hold off;

    file_name = strcat(complex, '_', figure_tag, '.eps');
    file_path = strcat(figure_dir, file_name);
    print('-depsc', file_path);
    eps_to_pdf(file_path);
    fprintf('done in %g sec!\n', toc);

end


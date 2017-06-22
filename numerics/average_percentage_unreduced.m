% test_speed_standard_reduction.m
init;

figure_dir = './figures/';
experiment_tag = 'average_percentage_unreduced';

LW = 'LineWidth';
MS = 'MarkerSize';
markers = '+o*xsd^v><ph.';

% Explore complexity of Vietoris-Rips complexes
vr_complexes = {'random_gaussian'};

% Algorithms to test
algorithms = {'alpha_beta_parallel'};

% Matrix dense?
as_dense = true;

color_list = create_color_palette(length(algorithms));

% Homology mode
homology_mode = 'reduced'

% Complex parameters
max_dim = 5;
num_divs= 10;
max_filtration_value = 5;

% Number of complexes per parameter
num_samples = 10;

% num_points
num_points = 15;

time_init = tic;

for i = 1:length(vr_complexes)

    complex = vr_complexes{i};
    complex_tag = strrep(complex, '_', '\_');

    % ---------------------
    % One figure per complex and per metric
    % ---------------------

    % pecentage_reduced

    percentage_unreduced = {};
    average_per_iteration = zeros(1, 7); % 6 iterations 

    % labels_ind

    label_ind = [];

    mfv = max_filtration_value;
    fprintf('\n\t%s @ (max_dim, num_divs, mfv) = (%d, %d, %d)\n',...
        complex, max_dim, num_divs, mfv);

    for k = 1:num_samples

        stream = example_factory(complex, max_dim, mfv, num_divs, num_points);
        T = BoundaryMatrixFactory(stream, homology_mode, 'testing', as_dense);
        lows_test = reduce_matrix(T, 'testing');
        fprintf('\t\tSample %d/%d\tm = %d\n', k, num_samples, T.m);

        for l = 1:length(algorithms)

            algo = algorithms{l};

            D = BoundaryMatrixFactory(stream, homology_mode, algo, as_dense);
            [lows, t] = reduce_matrix(D, algo);

            fprintf('\t\t\t%s... ', algo);

            if k == 1
                label_ind(end + 1) = true;
            else
                label_ind(end + 1) = false;
            end

            % --------------
            % Extract metrics
            % --------------

            metrics = D.metrics;
            x = 1:metrics.iters; 

            if isequal(algo, 'alpha_beta_parallel')
                style = '-+';
            else
                style = '--';
            end

            % --------------
            % percentage_unreduced
            % --------------

            percentage_unreduced{end + 1} = {x, metrics.percentage_unreduced(x)};

            % Interested in average after 1, 3, 5 iterations
            v = metrics.percentage_unreduced(1:7)
            average_per_iteration = average_per_iteration + v;

            % Assert output is correct

            assert(all(lows == lows_test), 'Output incorrect!');
            fprintf('\t\tsuccess in %g secs!\n', t);

        end % end algorithms

    end % end num_samples

    ind = find(label_ind);

    % Compute average

    average_per_iteration = average_per_iteration/num_samples

    % num_column_adds cumulative

    figure_tag = strcat(experiment_tag, '-', complex, '-', 'num_col_adds_cumulative', '.eps');
    filepath = fullfile(figure_dir, figure_tag);

end % end vr_complexes

% --------------
% End
% --------------

time_total = toc(time_init);
fprintf('All tests finished successfully in %s secs :)\n', num2str(time_total));


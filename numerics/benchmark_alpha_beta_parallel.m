% test_speed_standard_reduction.m
init;

figure_dir = './figures/';
experiment_tag = 'benchmark_alpha_beta_parallel';

LW = 'LineWidth';
MS = 'MarkerSize';
markers = '+o*xsd^v><ph.';

% Explore complexity of Vietoris-Rips complexes
vr_complexes = {'house', 'random_figure_8', ...
                'random_torus', 'sphere_product', ...
                'icosahedron', 'random_trefoil_knot', ...
                'random_gaussian', 'morozov'};

% Algorithms to test
algorithms = {'alpha_beta_parallel', 'std', 'twist', 'alpha_beta_std', 'alpha_beta_twist'};

color_list = create_color_palette(length(algorithms));

% Matrix dense?
as_dense = true;

% Complex parameters
max_dim = 5;
num_divs= 10;
max_filtration_value = 5;

% Number of complexes per parameter
num_samples = 3;

% num_points
num_points = 15;

time_init = tic;

for i = 1:length(vr_complexes)

    complex = vr_complexes{i};

    complex_tag = strrep(complex, '_', '\_');

    % ---------------------
    % One figure per complex and per metric
    % ---------------------

    % num_column_adds

    figure(1);
    set(gcf, 'color', [1 1 1]);
    set(gca, 'Fontname', 'setTimes', 'Fontsize', 15);
    handles_1 = [];
    labels_1 = {};

    % num_entry_adds

    figure(2);
    set(gcf, 'color', [1 1 1]);
    set(gca, 'Fontname', 'setTimes', 'Fontsize', 15);
    handles_2 = [];
    labels_2 = {};

    % pecentage_reduced

    figure(3);
    set(gcf, 'color', [1 1 1]);
    set(gca, 'Fontname', 'setTimes', 'Fontsize', 15);
    handles_3 = [];
    labels_3 = {};

    % num_column_adds

    figure(4);
    set(gcf, 'color', [1 1 1]);
    set(gca, 'Fontname', 'setTimes', 'Fontsize', 15);
    handles_4 = [];
    labels_4 = {};

    % num_entry_adds

    figure(5);
    set(gcf, 'color', [1 1 1]);
    set(gca, 'Fontname', 'setTimes', 'Fontsize', 15);
    handles_5 = [];
    labels_5 = {};

    % Essential estimation

    figure(6);
    set(gcf, 'color', [1 1 1]);
    set(gca, 'Fontname', 'setTimes', 'Fontsize', 15);
    handles_6 = [];
    labels_6 = {};

    % Lowstar L1

    figure(7);
    set(gcf, 'color', [1 1 1]);
    set(gca, 'Fontname', 'setTimes', 'Fontsize', 15);
    handles_7 = [];
    labels_7 = {};

    % labels_ind

    label_ind = [];

    mfv = max_filtration_value;
    fprintf('\n\t%s @ (max_dim, num_divs, mfv) = (%d, %d, %d)\n',...
        complex, max_dim, num_divs, mfv);

    for k = 1:num_samples

        stream = example_factory(complex, max_dim, mfv, num_divs, num_points);
        [lows_test, ~, T]= reduce_stream(stream, 'testing', as_dense);

        fprintf('\t\tSample %d/%d\tm = %d\n', k, num_samples, T.m);

        for l = 1:length(algorithms)

            algo = algorithms{l};
            [lows, t, D] = reduce_stream(stream, algo, as_dense, lows_test);
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
            % num_column_adds 
            % --------------

            figure(1);
            labels_1{end + 1} = strrep(algo, '_', '\_');

            y = metrics.num_column_adds(x);
            handles_1(end + 1) = loglog(x, y, style, 'Color', color_list{l});
            hold on;

            % Arrow
            ll=ceil(length(x)*(k-1/2)/num_samples);
            txt=['\leftarrow m=' num2str(T.m)];
            h=text(x(ll), y(ll), txt, 'HorizontalAlignment', 'left');
            set(h, 'Color', color_list{l});

            % --------------
            % num_entry_adds 
            % --------------

            figure(2);
            labels_2{end + 1} = strrep(algo, '_', '\_');

            y = metrics.num_entry_adds(x);
            handles_2(end + 1) = loglog(x, y, style, 'Color', color_list{l});
            hold on;

            % Arrow
            ll=ceil(length(x)*(k-1/2)/num_samples);
            txt=['\leftarrow m=' num2str(T.m)];
            h=text(x(ll), y(ll), txt, 'HorizontalAlignment', 'left');
            set(h, 'Color', color_list{l});

            % --------------
            % percentage_unreduced
            % --------------

            figure(3)
            labels_3{end + 1} = strrep(algo, '_', '\_');

            y = metrics.percentage_unreduced(x);
            handles_3(end + 1) = loglog(x, y, style, 'Color', color_list{l});
            hold on;

            % Arrow
            ll=ceil(length(x)*(k-1/2)/num_samples);
            txt=['\leftarrow m=' num2str(T.m)];
            h=text(x(ll), y(ll), txt, 'HorizontalAlignment', 'left');
            set(h, 'Color', color_list{l});

            % --------------
            % num_column_adds cumulative 
            % --------------

            figure(4);
            labels_4{end + 1} = strrep(algo, '_', '\_');

            y = cumsum(metrics.num_column_adds(x));
            handles_4(end + 1) = loglog(x, y, style, 'Color', color_list{l});
            hold on;

            % Arrow
            ll=ceil(length(x)*(k-1/2)/num_samples);
            txt=['\leftarrow m=' num2str(T.m)];
            h=text(x(ll), y(ll), txt, 'HorizontalAlignment', 'left');
            set(h, 'Color', color_list{l});

            % --------------
            % num_entry_adds cumulative 
            % --------------

            figure(5);
            labels_5{end + 1} = strrep(algo, '_', '\_');

            y = cumsum(metrics.num_entry_adds(x));
            handles_5(end + 1) = loglog(x, y, style, 'Color', color_list{l});
            hold on;

            % Arrow
            ll=ceil(length(x)*(k-1/2)/num_samples);
            txt=['\leftarrow m=' num2str(T.m)];
            h=text(x(ll), y(ll), txt, 'HorizontalAlignment', 'left');
            set(h, 'Color', color_list{l});

            % --------------
            % Essential estimation (precision)
            % --------------

            figure(6);
            labels_6{end + 1} = strrep(algo, '_', '\_');

            y = metrics.essential_precision(x);
            handles_6(end + 1) = loglog(x, y, style, 'Color', color_list{l});
            hold on;

            % Arrow
            ll=ceil(length(x)*(k-1/2)/num_samples);
            txt=['\leftarrow m=' num2str(T.m)];
            h=text(x(ll), y(ll), txt, 'HorizontalAlignment', 'left');
            set(h, 'Color', color_list{l});

            % --------------
            % Lowstar distance (L1)
            % --------------

            figure(7);
            labels_7{end + 1} = strrep(algo, '_', '\_');

            y = metrics.lowstar.l1(x);
            handles_7(end + 1) = loglog(x, y, style, 'Color', color_list{l});
            hold on;

            % Arrow
            ll=ceil(length(x)*(k-1/2)/num_samples);
            txt=['\leftarrow m=' num2str(T.m)];
            h=text(x(ll), y(ll), txt, 'HorizontalAlignment', 'left');
            set(h, 'Color', color_list{l});

            % Assert output is correct

            assert(all(lows == lows_test), 'Output incorrect!');
            fprintf('\t\tsuccess in %g secs!\n', t);

        end % end algorithms

    end % end num_samples

    ind = find(label_ind);

    % num_column_adds

    figure_tag = strcat(experiment_tag, '-', complex, '-', 'num_col_adds', '.eps');
    filepath = fullfile(figure_dir, figure_tag);

    figure(1);
    hold off;
    legend(handles_1(ind), labels_1(ind));
    xlabel('iteration');
    ylabel('column additions');
    title({'Number of column additions', complex_tag});

    print('-depsc', filepath);
    eps_to_pdf(filepath);
    close(1);

    % num_entry_adds

    figure_tag = strcat(experiment_tag, '-', complex, '-', 'num_entry_adds', '.eps');
    filepath = fullfile(figure_dir, figure_tag);

    figure(2);
    hold off;
    legend(handles_2(ind), labels_2(ind));
    xlabel('iteration');
    ylabel('entries changed in column additions');
    title({'Number of entries changed in column additions', complex_tag});

    print('-depsc', filepath);
    eps_to_pdf(filepath);
    close(2);

    % percentage_reduced

    figure_tag = strcat(experiment_tag, '-', complex, '-', 'percentage_unreduced', '.eps');
    filepath = fullfile(figure_dir, figure_tag);

    figure(3);
    hold off;
    legend(handles_3(ind), labels_3(ind));
    xlabel('iteration');
    ylabel('% of unreduced columns');
    title({'Percentage of unreduced columns', complex_tag});

    print('-depsc', filepath);
    eps_to_pdf(filepath);
    close(3);

    % num_column_adds cumulative

    figure_tag = strcat(experiment_tag, '-', complex, '-', 'num_col_adds_cumulative', '.eps');
    filepath = fullfile(figure_dir, figure_tag);

    figure(4);
    hold off;
    legend(handles_4(ind), labels_4(ind));
    xlabel('iteration');
    ylabel('column additions (cumulative)');
    title({'Number of column additions (cumulative)', complex_tag});

    print('-depsc', filepath);
    eps_to_pdf(filepath);
    close(4);

    % num_entry_adds cumulative

    figure_tag = strcat(experiment_tag, '-', complex, '-', 'num_entry_adds_cumulative', '.eps');
    filepath = fullfile(figure_dir, figure_tag);

    figure(5);
    hold off;
    legend(handles_5(ind), labels_5(ind));
    xlabel('iteration');
    ylabel('entries changed in column additions (cumulative)');
    title({'Number of entries changed in column additions (cumulative)', complex_tag});

    print('-depsc', filepath);
    eps_to_pdf(filepath);
    close(5);

    % essential estimation

    figure_tag = strcat(experiment_tag, '-', complex, '-', 'essential_estimation_precision', '.eps');
    filepath = fullfile(figure_dir, figure_tag);

    figure(6);
    hold off;
    legend(handles_6(ind), labels_6(ind));
    xlabel('iteration');
    ylabel('Precision TP/(TP + FP)');
    title({'Positive Predictive Value (Precision) of Essential Estimation', complex_tag});

    print('-depsc', filepath);
    eps_to_pdf(filepath);
    close(6);

    % Lowstar L1 distance

    figure_tag = strcat(experiment_tag, '-', complex, '-', 'lowstar_l1_distance', '.eps');
    filepath = fullfile(figure_dir, figure_tag);

    figure(7);
    hold off;
    legend(handles_7(ind), labels_7(ind));
    xlabel('iteration');
    ylabel('||low - lowstar||_1');
    title({'L1 distance between lowstar and low at each iteration', complex_tag});

    print('-depsc', filepath);
    eps_to_pdf(filepath);
    close(7);

end % end vr_complexes

% --------------
% End
% --------------

time_total = toc(time_init);
fprintf('All tests finished successfully in %s secs :)\n', num2str(time_total));


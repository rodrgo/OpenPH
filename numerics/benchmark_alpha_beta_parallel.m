% test_speed_standard_reduction.m
init;

figure_dir = './figures/';
experiment_tag = 'lookBack_benchmark_alpha_beta_parallel';

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

% Homology mode
homology_modes = {'reduced', 'unreduced'};
homology_modes = {'reduced'};

% Complex parameters
max_dim = 5;
num_divs= 10;
max_filtration_value = 5;

% Number of complexes per parameter
num_samples = 3;

% num_points
num_points = 15;

time_init = tic;

for h = 1:length(homology_modes)

    homology_mode = homology_modes{h};
    fprintf('homology_mode = %s\n\n', homology_mode);

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

        % labels_ind

        label_ind = [];

        mfv = max_filtration_value;
        fprintf('\n\t%s @ (max_dim, num_divs, mfv) = (%d, %d, %d)\n',...
            complex, max_dim, num_divs, mfv);

        for k = 1:num_samples

            stream = example_factory(complex, max_dim, mfv, num_divs, num_points);
            lows_test = reduce_stream(stream, homology_mode, 'testing', as_dense);

            fprintf('\t\tSample %d/%d\tm = %d\n', k, num_samples, T.m);

            for l = 1:length(algorithms)

                algo = algorithms{l};
                [lows, t] = reduce_stream(stream, homology_mode, algo, as_dense);
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
                handles_1(end + 1) = plot(x, metrics.num_column_adds(x), style, 'Color', color_list{l});
                hold on;

                % --------------
                % num_entry_adds 
                % --------------

                figure(2);
                labels_2{end + 1} = strrep(algo, '_', '\_');
                handles_2(end + 1) = semilogy(x, metrics.num_entry_adds(x), style, 'Color', color_list{l});
                hold on;

                % --------------
                % percentage_unreduced
                % --------------

                figure(3)
                labels_3{end + 1} = strrep(algo, '_', '\_');
                handles_3(end + 1) = plot(x, metrics.percentage_unreduced(x), style, 'Color', color_list{l});
                hold on;

                % --------------
                % num_column_adds cumulative 
                % --------------

                figure(4);
                labels_4{end + 1} = strrep(algo, '_', '\_');
                handles_4(end + 1) = plot(x, cumsum(metrics.num_column_adds(x)), style, 'Color', color_list{l});
                hold on;

                % --------------
                % num_entry_adds cumulative 
                % --------------

                figure(5);
                labels_5{end + 1} = strrep(algo, '_', '\_');
                handles_5(end + 1) = semilogy(x, cumsum(metrics.num_entry_adds(x)), style, 'Color', color_list{l});
                hold on;

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

    end % end vr_complexes

end % end homology_mode

% --------------
% End
% --------------

time_total = toc(time_init);
fprintf('All tests finished successfully in %s secs :)\n', num2str(time_total));


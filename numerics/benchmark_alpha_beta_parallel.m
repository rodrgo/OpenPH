% test_speed_standard_reduction.m
init;

figure_dir = './figures/';
experiment_tag = 'benchmark_alpha_beta_parallel';

LW = 'LineWidth';
MS = 'MarkerSize';
markers = '+o*xsd^v><ph.';

% Font size
fs = [];
fs.title = 20;
fs.legend = 17;
fs.axis = 20;
fs.ticks = 20;


% Explore complexity of Vietoris-Rips complexes
vr_complexes = {'house', 'random_figure_8', ...
                'random_torus', 'sphere_product', ...
                'icosahedron', 'random_trefoil_knot', ...
                'random_gaussian', 'morozov'};

vr_complexes = {'random_gaussian', 'random_figure_8', ...
                'random_trefoil_knot', ...
                'sphere_product'};

% Algorithms to test
%algorithms = {'alpha_beta_parallel', 'std', 'twist', 'alpha_beta_std', 'alpha_beta_twist'};
algorithms = {'alpha_beta_parallel', 'std', 'twist'};

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

% ---------------------
% Create placeholders for tables
% ---------------------

% Table experiments:
num_experiments = 3;

num_complexes = length(vr_complexes);
num_algos = length(algorithms); % num of algorithms

% Essential estimation
ess_percentiles = [0.1, 0.5, 0.9, 0.95, 1];
tensor_ess = zeros(num_complexes, num_samples, num_algos, length(ess_percentiles));

% Percentage of unreduced columns
unreduced_percentiles = [0.9, 0.5, 0.1, 0.05, 0.01];
tensor_unreduced = zeros(num_complexes, num_samples, num_algos, length(unreduced_percentiles));

% L1 distance
l1_percentiles = [1e-1, 1e-2, 1e-3, 1e-4, 0];
tensor_l1 = zeros(num_complexes, num_samples, num_algos, length(l1_percentiles));

% Total column operation difference
tensor_col_ops = zeros(num_complexes, num_samples, num_algos); % last entry stores total number of column operations

for i = 1:length(vr_complexes)

    complex = vr_complexes{i};

    complex_tag = strrep(complex, '_', '\_');

    % ---------------------
    % One figure per complex and per metric
    % ---------------------

    % num_column_adds

    figure(1);
    set(gcf, 'color', [1 1 1]);
    set(gca, 'Fontname', 'setTimes', 'Fontsize', 18);
    handles_1 = [];
    labels_1 = {};

    % num_entry_adds

    figure(2);
    set(gcf, 'color', [1 1 1]);
    set(gca, 'Fontname', 'setTimes', 'Fontsize', 18);
    handles_2 = [];
    labels_2 = {};

    % pecentage_reduced

    figure(3);
    set(gcf, 'color', [1 1 1]);
    set(gca, 'Fontname', 'setTimes', 'Fontsize', 18);
    handles_3 = [];
    labels_3 = {};

    % num_column_adds

    figure(4);
    set(gcf, 'color', [1 1 1]);
    set(gca, 'Fontname', 'setTimes', 'Fontsize', 18);
    handles_4 = [];
    labels_4 = {};

    % num_entry_adds

    figure(5);
    set(gcf, 'color', [1 1 1]);
    set(gca, 'Fontname', 'setTimes', 'Fontsize', 18);
    handles_5 = [];
    labels_5 = {};

    % Essential estimation

    figure(6);
    set(gcf, 'color', [1 1 1]);
    set(gca, 'Fontname', 'setTimes', 'Fontsize', 18);
    handles_6 = [];
    labels_6 = {};

    % Lowstar L1

    figure(7);
    set(gcf, 'color', [1 1 1]);
    set(gca, 'Fontname', 'setTimes', 'Fontsize', 18);
    handles_7 = [];
    labels_7 = {};

    % labels_ind

    label_ind = [];

    mfv = max_filtration_value;
    fprintf('\n\t%s @ (max_dim, num_divs, mfv) = (%d, %d, %d)\n',...
        complex, max_dim, num_divs, mfv);

    for k = 1:num_samples

        [stream, complex_cell_info] = example_factory(complex, max_dim, mfv, num_divs, num_points);
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
            labels_1{end + 1} = strrep(name_change(algo), '_', '\_');

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
            labels_2{end + 1} = strrep(name_change(algo), '_', '\_');

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
            labels_3{end + 1} = strrep(name_change(algo), '_', '\_');

            y = metrics.percentage_unreduced(x);
            handles_3(end + 1) = loglog(x, y, style, 'Color', color_list{l});
            hold on;

            % Arrow
            ll=ceil(length(x)*(k-1/2)/num_samples);
            txt=['\leftarrow m=' num2str(T.m)];
            h=text(x(ll), y(ll), txt, 'HorizontalAlignment', 'left');
            set(h, 'Color', color_list{l});

            % Store tensor info

            % Number of iterations for 90 percent reduction
            for pp = 1:length(unreduced_percentiles)
            	tensor_unreduced(i, k, l, pp) = find(y - unreduced_percentiles(pp) <= 0 , 1, 'first');
            end

            % --------------
            % num_column_adds cumulative 
            % --------------

            figure(4);
            labels_4{end + 1} = strrep(name_change(algo), '_', '\_');

            y = cumsum(metrics.num_column_adds(x));
            handles_4(end + 1) = loglog(x, y, style, 'Color', color_list{l});
            hold on;

            % Arrow
            ll=ceil(length(x)*(k-1/2)/num_samples);
            txt=['\leftarrow m=' num2str(T.m)];
            h=text(x(ll), y(ll), txt, 'HorizontalAlignment', 'left');
            set(h, 'Color', color_list{l});

            % Store tensor info
            tensor_col_ops(i, k, l) = y(x(end));

            % --------------
            % num_entry_adds cumulative 
            % --------------

            figure(5);
            labels_5{end + 1} = strrep(name_change(algo), '_', '\_');

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
            labels_6{end + 1} = strrep(name_change(algo), '_', '\_');

            y = metrics.essential_precision(x);
            handles_6(end + 1) = loglog(x, y, style, 'Color', color_list{l});
            hold on;

            % Arrow
            ll=ceil(length(x)*(k-1/2)/num_samples);
            txt=['\leftarrow m=' num2str(T.m)];
            h=text(x(ll), y(ll), txt, 'HorizontalAlignment', 'left');
            set(h, 'Color', color_list{l});

            % Store tensor info
            for pp = 1:length(ess_percentiles)
                ep = find(y - ess_percentiles(pp) >= 0 , 1, 'first');
            	tensor_ess(i, k, l, pp) = ep;
            end

            % --------------
            % Lowstar distance (L1)
            % --------------

            figure(7);
            labels_7{end + 1} = strrep(name_change(algo), '_', '\_');

            y = metrics.lowstar.l1(x);
            handles_7(end + 1) = loglog(x, y, style, 'Color', color_list{l});
            hold on;

            % Arrow
            ll=ceil(length(x)*(k-1/2)/num_samples);
            txt=['\leftarrow m=' num2str(T.m)];
            h=text(x(ll), y(ll), txt, 'HorizontalAlignment', 'left');
            set(h, 'Color', color_list{l});

            % Store tensor info
            for pp = 1:length(l1_percentiles)
            	tensor_l1(i, k, l, pp) = find(y - l1_percentiles(pp) <= 0 , 1, 'first');
            end

            % Assert output is correct

            assert(all(lows == lows_test), 'Output incorrect!');
            fprintf('\t\tsuccess in %g secs!\n', t);

        end % end algorithms

    end % end num_samples

    ind = find(label_ind);

    complex_str = complex_cell_info{1};
    params_str = complex_cell_info{2}; 
    if length(params_str) == 0;
        second_line_str = [complex_str];
    else
        second_line_str = [complex_str, ': ', params_str];
    end

    % num_column_adds

    figure_tag = strcat(experiment_tag, '-', complex, '-', 'num_col_adds', '.eps');
    filepath = fullfile(figure_dir, figure_tag);

    figure(1);
    hold off;
    legend(handles_1(ind), labels_1(ind), 'Location', 'northeast', 'FontSize', fs.legend);
    xlabel('iteration', 'FontSize', fs.axis);
    ylabel('count of column additions', 'FontSize', fs.axis);
			% Tick size
			xt = get(gca, 'XTick');
			set(gca, 'FontSize', fs.ticks);

			xt = get(gca, 'YTick');
			set(gca, 'FontSize', fs.ticks);

    title(second_line_str, 'FontSize', fs.title);
    %title({'Number of column additions', second_line_str});

    print('-depsc', filepath);
    eps_to_pdf(filepath);
    close(1);

    % num_entry_adds

    figure_tag = strcat(experiment_tag, '-', complex, '-', 'num_entry_adds', '.eps');
    filepath = fullfile(figure_dir, figure_tag);

    figure(2);
    hold off;
    legend(handles_2(ind), labels_2(ind), 'Location', 'southeast', 'FontSize', fs.legend);
    xlabel('iteration', 'FontSize', fs.axis);
    ylabel('count of XOR operations', 'FontSize', fs.axis);
			% Tick size
			xt = get(gca, 'XTick');
			set(gca, 'FontSize', fs.ticks);

			xt = get(gca, 'YTick');
			set(gca, 'FontSize', fs.ticks);

    title(second_line_str, 'FontSize', fs.title);
    %title({'Number of entries changed in column additions', second_line_str});

    print('-depsc', filepath);
    eps_to_pdf(filepath);
    close(2);

    % percentage_reduced

    figure_tag = strcat(experiment_tag, '-', complex, '-', 'percentage_unreduced', '.eps');
    filepath = fullfile(figure_dir, figure_tag);

    figure(3);
    hold off;
    legend(handles_3(ind), labels_3(ind), 'FontSize', fs.legend);
    xlabel('iteration', 'FontSize', fs.axis);
    ylabel('% of unreduced columns', 'FontSize', fs.axis);
			% Tick size
			xt = get(gca, 'XTick');
			set(gca, 'FontSize', fs.ticks);

			xt = get(gca, 'YTick');
			set(gca, 'FontSize', fs.ticks);

    title(second_line_str, 'FontSize', fs.title);
    %title({'Percentage of unreduced columns', second_line_str});

    print('-depsc', filepath);
    eps_to_pdf(filepath);
    close(3);

    % num_column_adds cumulative

    figure_tag = strcat(experiment_tag, '-', complex, '-', 'num_col_adds_cumulative', '.eps');
    filepath = fullfile(figure_dir, figure_tag);

    figure(4);
    hold off;
    legend(handles_4(ind), labels_4(ind), 'Location', 'southeast', 'FontSize', fs.legend);
    xlabel('iteration', 'FontSize', fs.axis);
    ylabel('count of column additions (cumulative)', 'FontSize', fs.axis);
			% Tick size
			xt = get(gca, 'XTick');
			set(gca, 'FontSize', fs.ticks);

			xt = get(gca, 'YTick');
			set(gca, 'FontSize', fs.ticks);

    title(second_line_str, 'FontSize', fs.title);
    %title({'Number of column additions (cumulative)', second_line_str});

    print('-depsc', filepath);
    eps_to_pdf(filepath);
    close(4);

    % num_entry_adds cumulative

    figure_tag = strcat(experiment_tag, '-', complex, '-', 'num_entry_adds_cumulative', '.eps');
    filepath = fullfile(figure_dir, figure_tag);

    figure(5);
    hold off;
    legend(handles_5(ind), labels_5(ind), 'Location', 'southeast', 'FontSize', fs.legend);
    xlabel('iteration', 'FontSize', fs.axis);
    ylabel('count of XOR operations (cumulative)', 'FontSize', fs.axis);
			% Tick size
			xt = get(gca, 'XTick');
			set(gca, 'FontSize', fs.ticks);

			xt = get(gca, 'YTick');
			set(gca, 'FontSize', fs.ticks);

    title(second_line_str, 'FontSize', fs.title);
    %title({'Number of entries changed in column additions (cumulative)', second_line_str});

    print('-depsc', filepath);
    eps_to_pdf(filepath);
    close(5);

    % essential estimation

    figure_tag = strcat(experiment_tag, '-', complex, '-', 'essential_estimation_precision', '.eps');
    filepath = fullfile(figure_dir, figure_tag);

    figure(6);
    hold off;
    legend(handles_6(ind), labels_6(ind), 'FontSize', fs.legend);
    xlabel('iteration', 'FontSize', fs.axis);
    ylabel('Positive Predictive Value', 'FontSize', fs.axis);
			% Tick size
			xt = get(gca, 'XTick');
			set(gca, 'FontSize', fs.ticks);

			xt = get(gca, 'YTick');
			set(gca, 'FontSize', fs.ticks);

    title(second_line_str, 'FontSize', fs.title);
    %title({'Positive Predictive Value of Essential Estimation', second_line_str});

    print('-depsc', filepath);
    eps_to_pdf(filepath);
    close(6);

    % Lowstar L1 distance

    figure_tag = strcat(experiment_tag, '-', complex, '-', 'lowstar_l1_distance', '.eps');
    filepath = fullfile(figure_dir, figure_tag);

    figure(7);
    hold off;
    legend(handles_7(ind), labels_7(ind), 'Location', 'southwest', 'FontSize', fs.legend);
    xlabel('iteration', 'FontSize', fs.axis);
    ylabel('$\frac{\|low - low^*\|_1}{\|low^*\|_1}$', 'Interpreter', 'LaTex', 'FontSize', fs.axis);
			% Tick size
			xt = get(gca, 'XTick');
			set(gca, 'FontSize', fs.ticks);

			xt = get(gca, 'YTick');
			set(gca, 'FontSize', fs.ticks);

    title(second_line_str, 'FontSize', fs.title);
    %title({'||low - low^*||_1/||low^*||_1 per iteration', second_line_str});

    print('-depsc', filepath);
    eps_to_pdf(filepath);
    close(7);

end % end vr_complexes

% --------------
% Print stats tables
% --------------

fileId = fopen('tables.tex', 'w');

% table i
for i = 1:num_experiments
    if i == 1
        levels_x = unreduced_percentiles;
        tensor_x = tensor_unreduced;
    elseif i == 2
        levels_x = ess_percentiles;
        tensor_x = tensor_ess;
    else
        levels_x = l1_percentiles;
        tensor_x = tensor_l1;
    end
    fprintf(fileId,'\n\n\n\n');

    fprintf(fileId,'\\begin{small}\n');
    fprintf(fileId,'\\begin{table*}\n');
    fprintf(fileId,'\\centering\n');
    fprintf(fileId,'\\begin{tabular}{l');
    fprintf(fileId,'||');
    for j = 1:num_complexes
        for l = 1:num_algos
            fprintf(fileId,'c');
        end
        if j < num_complexes
    	    fprintf(fileId,'|');
        end
    end

    fprintf(fileId,'}\n');
    fprintf(fileId,'\\toprule\n');
    if i == 1
        fprintf(fileId,'\\multirow{2}{*}{Proportion} &\n');
    elseif i == 2
        fprintf(fileId,'\\multirow{2}{*}{Precision} &\n');
    else
        fprintf(fileId,'\\multirow{2}{*}{$\\frac{\\|\\low - \\lowstar\\|_1}{\\|\\lowstar\\|_1}$} &\n');
    end
    for j = 1:num_complexes
        fprintf(fileId,'\\multicolumn{%d}{c}{%s}', num_algos, name_change_table(vr_complexes{j}));
        if j < num_complexes
            fprintf(fileId,'&\n');
        else
            fprintf(fileId,'\\\\\n');
        end
    end
    for j = 1:num_complexes
        for l = 1:num_algos
            fprintf(fileId,'& {%s} ', name_change(algorithms{l}));
        end
    end
    fprintf(fileId,'\\\\\n');
    fprintf(fileId,'\\midrule\n');

    % row 'it'
    for pp = 1:length(levels_x)
        if i == 1 || i == 2
            fprintf(fileId,'%1.2f', levels_x(pp));
        else
            fprintf(fileId,'%g', levels_x(pp));
        end
        for j = 1:num_complexes
            for k = 1:num_algos
                fprintf(fileId,' & %d', tensor_x(j, 1, k, pp)); % Only do it for alpha_beta_parallel
            end
        end
        fprintf(fileId,'\\\\\n');
    end

    fprintf(fileId,'\\bottomrule\n');
    fprintf(fileId,'\\end{tabular}\n');
    if i == 1
        fprintf(fileId,'\\caption{Iterations to unreduced percentage}\n');
	fprintf(fileId,'\\label{tab:iterations_unreduced}\n');
    elseif i == 2
        fprintf(fileId,'\\caption{Iterations to essential-estimation precision}\n');
	fprintf(fileId,'\\label{tab:iterations_essential}\n');
    else
        fprintf(fileId,'\\caption{Iterations to relative $\\ell_1$-error}\n');
	fprintf(fileId,'\\label{tab:iterations_l1_error}\n');
    end
    fprintf(fileId,'\\end{table*}\n');
    fprintf(fileId,'\\end{small}\n');

end

% --------------
% Print savings table
% --------------

fprintf(fileId,'\n\n\n\n');
fprintf(fileId,'\\begin{small}\n');
fprintf(fileId,'\\begin{table*}\n');
fprintf(fileId,'\\centering\n');
fprintf(fileId,'\\begin{tabular}{l');
fprintf(fileId,'||');
for j = 1:num_complexes
    for l = 1:num_algos
        if strcmp(algorithms{l}, 'alpha_beta_parallel') == 0
            fprintf(fileId,'c');
        end
    end
    if j < num_complexes
        fprintf(fileId,'|');
    end
end
fprintf(fileId,'}\n');
fprintf(fileId,'\\toprule\n');
fprintf(fileId,'\\multirow{2}{*}{Sample} &\n');
for j = 1:num_complexes
    fprintf(fileId,'\\multicolumn{%d}{c}{%s}', num_algos-1, name_change_table(vr_complexes{j}));
    if j < num_complexes
        fprintf(fileId,'&\n');
    else
        fprintf(fileId,'\\\\\n');
    end
end
for j = 1:num_complexes
    for l = 1:num_algos
        if strcmp(algorithms{l}, 'alpha_beta_parallel') == 0
            fprintf(fileId,'& {%s} ', name_change(algorithms{l}));
        end
    end
end
fprintf(fileId,'\\\\\n');
fprintf(fileId,'\\midrule\n');

% row 'it'
for s = 1:num_samples
    fprintf(fileId,'%d', s);
    for j = 1:num_complexes
        for k = 1:num_algos
            if strcmp(algorithms{k}, 'alpha_beta_parallel') == 0
                ratio = tensor_col_ops(j, s, 1)/tensor_col_ops(j, s, k);
                fprintf(fileId,' & %1.2f', ratio);
            end
        end
    end
    fprintf(fileId,'\\\\\n');
end

fprintf(fileId,'\\bottomrule\n');
fprintf(fileId,'\\end{tabular}\n');
fprintf(fileId,'\\caption{Ratio of total column operations}\n');
fprintf(fileId,'\\label{tab:ratio_cumsum_operations}\n');
fprintf(fileId,'\\end{table*}\n');
fprintf(fileId,'\\end{small}\n');

fclose(fileId);

% --------------
% End
% --------------

time_total = toc(time_init);
fprintf('All tests finished successfully in %s secs :)\n', num2str(time_total));


% test_speed_standard_reduction.m
init;
plot_init;

% Complex parameters
cparams               = [];
cparams.max_dim       = 5;
cparams.num_divs      = 10;
cparams.max_filtr_val = 5;
cparams.num_points    = 15;

shapes = {'random_gaussian', ...
          'random_figure_8', ...
          'random_trefoil_knot'};

algos  = {'standard', ...
          'twist', ...
          'ph_row', ...
          'pms'};

time_init   = tic;
perc        = [];
perc.ess    = [0.1, 0.5, 0.9, 0.95, 1];
perc.unred  = [0.9, 0.5, 0.1, 0.05, 0.01];
perc.lone   = [1e-1, 1e-2, 1e-3, 1e-4, 0];
perc.linf   = [1e-1, 1e-2, 1e-3, 1e-4, 0];

num_shapes   = length(shapes);
num_algos    = length(algos);
num_samples  = 3;

COL_WIDTH    = 8;

tensor       = [];
tensor.ess   = zeros(num_shapes, num_samples, num_algos, length(perc.ess));
tensor.unred = zeros(num_shapes, num_samples, num_algos, length(perc.unred));
tensor.lone  = zeros(num_shapes, num_samples, num_algos, length(perc.lone));
tensor.linf  = zeros(num_shapes, num_samples, num_algos, length(perc.linf));

shapes_map  = containers.Map;

for i = 1:length(shapes)

    fprintf('\n\t%s: (max_dim, num_div, mfv) = (%d, %d, %d)\n',...
        shapes{i}, cparams.max_dim, cparams.num_divs, cparams.max_filtr_val);

    samples_map = containers.Map;
    for k = 1:num_samples

        [stream, complex_info] = complex_factory(shapes{i}, cparams);
        [r, c, m] = stream2cmo(stream);
        m = stream.getSize();
        low_true = get_true_low(r, c, m, COL_WIDTH);

        fprintf('\t\tSample %d/%d\tm = %d\n', k, num_samples, m);
        algos_map = containers.Map;
        for l = 1:length(algos)

            fprintf('\t\t\t%s... ', algos{l});
            [OUT, t] = openph(r, c, m, algos{l}, COL_WIDTH, low_true);
            fprintf('\t\tsuccess in %g secs!\n', t);
            tensor = populate_tensor(OUT, tensor, perc, i, k, l);
            algos_map(algos{l}) = OUT;

        end % end algorithms

        sample_struct = [];
        sample_struct.m = m;
        sample_struct.complex_info = complex_info;
        sample_struct.algos_map = algos_map;
        samples_map(num2str(k)) = sample_struct;

    end % end num_samples

    shapes_map(shapes{i}) = samples_map;

end % end shapes

benchmark_pms_plot(shapes_map, 'time', 'err_linf');
benchmark_pms_plot(shapes_map, 'time', 'err_lone');
benchmark_pms_plot(shapes_map, 'time', 'err_redu');
benchmark_pms_plot(shapes_map, 'time', 'err_ess');

benchmark_pms_plot(shapes_map, 'iters', 'err_linf');
benchmark_pms_plot(shapes_map, 'iters', 'err_lone');
benchmark_pms_plot(shapes_map, 'iters', 'err_redu');
benchmark_pms_plot(shapes_map, 'iters', 'err_ess');

% ---------------------
% Create tables
% ---------------------

create_table(shapes, algos, 'unreduced', perc.unred, tensor.unred);
create_table(shapes, algos, 'essential', perc.ess, tensor.ess);
create_table(shapes, algos, 'lone', perc.lone, tensor.lone);
create_table(shapes, algos, 'linf', perc.linf, tensor.linf);

% ---------------------
% Function definitions: populate_tensor
% ---------------------

function tensor = populate_tensor(OUT, tensor, perc, i, k, l)
    % ess
    y = OUT.err_ess(1:OUT.num_iters);
    for pp = 1:length(perc.ess)
        tensor.ess(i, k, l, pp) = find(y - perc.ess(pp) >= 0 , 1, 'first');
    end

    % unred
    y = OUT.err_redu(1:OUT.num_iters);
    for pp = 1:length(perc.unred)
        tensor.unred(i, k, l, pp) = find(y - perc.unred(pp) <= 0 , 1, 'first');
    end

    % lone
    y = OUT.err_lone(1:OUT.num_iters);
    for pp = 1:length(perc.lone)
        tensor.lone(i, k, l, pp) = find(y - perc.lone(pp) <= 0 , 1, 'first');
    end

    % linf
    y = OUT.err_linf(1:OUT.num_iters);
    for pp = 1:length(perc.linf)
        tensor.linf(i, k, l, pp) = find(y - perc.linf(pp) <= 0 , 1, 'first');
    end
end

% ---------------------
% Function definitions: benchmark_pms_plot
% ---------------------

function benchmark_pms_plot(shapes_map, x_unit, y_metric)
    plot_init;

    EXPERIMENT_TAG  = 'benchmark_pms';

    shapes = keys(shapes_map);
    for i = 1:length(shapes)
        shape = shapes{i};
        samples_map = shapes_map(shape);

        % label_ind
        label_ind = [];

        % Create figure
        figure(1);
        set(gcf, 'color', [1 1 1]);
        set(gca, 'Fontname', 'setTimes', 'Fontsize', 18);
        handles = [];
        labels  = {};

        samples = keys(samples_map);
        num_samples = length(samples);
        for k = 1:length(samples)
            sample = samples{k};
            sample_struct = samples_map(sample);

            algos_map = sample_struct.algos_map;
            m = sample_struct.m;
            complex_info = sample_struct.complex_info;

            algos = keys(algos_map);
            colors = create_color_palette(length(algos));
            for l = 1:length(algos)
                algo = algos{l};
                OUT = algos_map(algo);

                if k == 1
                    label_ind(end + 1) = true;
                else
                    label_ind(end + 1) = false;
                end

                if isequal(algo, 'pms')
                    style = '-+';
                else
                    style = '--';
                end

                % --------------
                % x_unit (iters)
                % --------------

                z = 1:OUT.num_iters;
                y = OUT.(y_metric)(z);

                if strcmp(x_unit, 'iters')
                    x = z;
                elseif strcmp(x_unit, 'time')
                    x = cumsum(OUT.time_track(z));
                else
                    error('not recognised');
                end

                % --------------
                % y_metric
                % --------------

                figure(1)
                labels{end + 1}  = strrep(algo, '_', '\_');
                handles(end + 1) = semilogy(x, y, style, 'Color', colors{l});
                hold on;

                % Arrow
                ll=ceil(length(x)*(k-1/2)/num_samples);
                txt=['\leftarrow m=' num2str(m)];
                h=text(x(ll), y(ll), txt, 'HorizontalAlignment', 'left');
                set(h, 'Color', colors{l});

            end
        end

        ind = find(label_ind);

        complex_str = complex_info{1};
        params_str = complex_info{2}; 
        if length(params_str) == 0;
            second_line_str = [complex_str];
        else
            second_line_str = [complex_str, ': ', params_str];
        end

        % figure

        figure_tag = strcat(EXPERIMENT_TAG, '-', shape, '-', x_unit, '-', y_metric, '.eps');
        filepath = fullfile(FIGURE_DIR, figure_tag);

        figure(1);
        hold off;
        legend(handles(ind), labels(ind), 'FontSize', fs.legend);


        if strcmp(x_unit, 'iters')
            xlab = 'Iteration';
        elseif strcmp(x_unit, 'time')
            xlab = 'Time (ms)';
        else
            error('not recognised');
        end
        xlabel(xlab, 'FontSize', fs.axis);

        if strcmp(y_metric, 'err_linf')
            ylabel('$\frac{\|low - low^*\|_\infty}{\|low^*\|_\infty}$', 'Interpreter', 'LaTex', 'FontSize', fs.axis);
        elseif strcmp(y_metric, 'err_lone')
            ylabel('$\frac{\|low - low^*\|_1}{\|low^*\|_1}$', 'Interpreter', 'LaTex', 'FontSize', fs.axis);
        elseif strcmp(y_metric, 'err_redu')
            ylabel('% of unreduced columns', 'FontSize', fs.axis);
        elseif strcmp(y_metric, 'err_ess')
            xlab = 'Time (ms)';
            ylabel('Precision of Ess estimation', 'FontSize', fs.axis);
        else
            error('not recognised');
        end

        % Tick size
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', fs.ticks);

        xt = get(gca, 'YTick');
        set(gca, 'FontSize', fs.ticks);

        title(second_line_str, 'FontSize', fs.title);

        print('-depsc', filepath);
        eps_to_pdf(filepath);
        close(1);

    end
end

% ---------------------
% Function definitions: create_table
% ---------------------

function create_table(shapes, algos, table_type, levels_x, tensor_x)
    plot_init;
    % table_type:
    %   + 'unreduced'
    %   + 'essential'
    %   + 'l1-error'
    %   + 'linf-error'

    EXPERIMENT_TAG = 'table';
    table_tag = strcat(EXPERIMENT_TAG, '_', table_type, '.tex');
    table_path = fullfile(FIGURE_DIR, table_tag);

    fileId = fopen(table_path, 'w');
    table_header(fileId);

    num_shapes = length(shapes);
    num_algos = length(algos);

    for j = 1:num_shapes
        for l = 1:num_algos
            fprintf(fileId,'c');
        end
        if j < num_shapes
            fprintf(fileId,'|');
        end
    end

    fprintf(fileId,'}\n');
    fprintf(fileId,'\\toprule\n');

    table_multirow(fileId, table_type);

    for j = 1:num_shapes
        fprintf(fileId,'\\multicolumn{%d}{c}{%s}', num_algos, rename_shape(shapes{j}));
        if j < num_shapes
            fprintf(fileId,'&\n');
        else
            fprintf(fileId,'\\\\\n');
        end
    end
    for j = 1:num_shapes
        for l = 1:num_algos
            fprintf(fileId,'& {%s} ', rename_alg(algos{l}));
        end
    end
    fprintf(fileId,'\\\\\n');
    fprintf(fileId,'\\midrule\n');

    % row 'it'
    for pp = 1:length(levels_x)
        if strcmp(table_type, 'unreduced') || strcmp(table_type, 'essential')
            fprintf(fileId,'%1.2f', levels_x(pp));
        else
            fprintf(fileId,'%g', levels_x(pp));
        end
        for j = 1:num_shapes
            for k = 1:num_algos
                fprintf(fileId,' & %d', tensor_x(j, 1, k, pp));
            end
        end
        fprintf(fileId,'\\\\\n');
    end

    fprintf(fileId,'\\bottomrule\n');
    fprintf(fileId,'\\end{tabular}\n');

end

function name = rename_shape(shape)
    name = '';
    if strcmp(shape, 'random_gaussian') == 1
       name = 'Gaussian';
    elseif strcmp(shape, 'random_figure_8') == 1
       name = 'Figure-8';
    elseif strcmp(shape, 'random_trefoil_knot') == 1
       name = 'Trefoil-Knot';
    elseif strcmp(shape, 'random_torus') == 1
       name = 'Random-Torus';
    elseif strcmp(shape, 'random_sphere_product') == 1
       name = 'Sphere-Product';
    else
       error('In rename_shape, shape not recognised');
    end
end

function name = rename_alg(alg)
    name = '';
    if strcmp(alg, 'standard') == 1
       name = 'std';
    elseif strcmp(alg, 'twist') == 1
       name = 'twist';
    elseif strcmp(alg, 'ph_row') == 1
       name = 'phRow';
    elseif strcmp(alg, 'standard_parallel') == 1
       name = 'std-parallel';
    elseif strcmp(alg, 'twist_parallel') == 1
       name = 'twist-parallel';
    elseif strcmp(alg, 'ph_row_parallel') == 1
       name = 'phRow-parallel';
    elseif strcmp(alg, 'pms') == 1
       name = 'pms';
    else
       error('In rename_alg, alg not recognised');
    end
end

function table_header(fileId)
    fprintf(fileId,'\\begin{tabular}{l');
    fprintf(fileId,'||');
end

function table_caption(fileId, table_type)
    if strcmp(table_type, 'unreduced')
        fprintf(fileId,'\\caption{Iterations to unreduced percentage}\n');
        fprintf(fileId,'\\label{tab:iterations_unreduced}\n');
    elseif strcmp(table_type, 'essential')
        fprintf(fileId,'\\caption{Iterations to essential-estimation precision}\n');
        fprintf(fileId,'\\label{tab:iterations_essential}\n');
    elseif strcmp(table_type, 'lone')
        fprintf(fileId,'\\caption{Iterations to relative $\\ell_1$-error}\n');
        fprintf(fileId,'\\label{tab:iterations_l1_error}\n');
    elseif strcmp(table_type, 'linf')
        fprintf(fileId,'\\caption{Iterations to relative $\\ell_\\infty$-error}\n');
        fprintf(fileId,'\\label{tab:iterations_linf_error}\n');
    else
        error('table_type not recognised');
    end
end

function table_multirow(fileId, table_type)
    if strcmp(table_type, 'unreduced')
        fprintf(fileId,'\\multirow{2}{*}{Proportion} &\n');
    elseif strcmp(table_type, 'essential')
        fprintf(fileId,'\\multirow{2}{*}{Precision} &\n');
    elseif strcmp(table_type, 'lone')
        fprintf(fileId,'\\multirow{2}{*}{$\\frac{\\|\\low - \\lowstar\\|_1}{\\|\\lowstar\\|_1}$} &\n');
    elseif strcmp(table_type, 'linf')
        fprintf(fileId,'\\multirow{2}{*}{$\\frac{\\|\\low - \\lowstar\\|_\\infty}{\\|\\lowstar\\|_\\infty}$} &\n');
    else
        error('table_type not recognised');
    end
end


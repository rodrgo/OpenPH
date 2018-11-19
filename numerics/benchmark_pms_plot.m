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
            handles(end + 1) = loglog(x, y, style, 'Color', colors{l});
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

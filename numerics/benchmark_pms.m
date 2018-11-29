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


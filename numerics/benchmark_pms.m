% test_speed_standard_reduction.m
init;

FIGURE_DIR      = '../figures/';
EXPERIMENT_TAG  = 'benchmark_pms';

% Complex parameters
cparams               = [];
cparams.max_dim       = 5;
cparams.num_divs      = 10;
cparams.max_filtr_val = 5;
cparams.num_points    = 15;

num_samples = 3;
num_samples = 1;

% Explore complexity of Vietoris-Rips complexes
shapes = {'random_gaussian', ...
          'random_figure_8', ...
          'random_trefoil_knot'};

shapes = {'random_gaussian', ...
          'random_figure_8'};

algos = {'standard_parallel', 'twist_parallel', 'ph_row_parallel', 'pms'};

% Matrix dense?
as_dense  = true;
time_init = tic;

shapes_map = containers.Map;
for i = 1:length(shapes)

    fprintf('\n\t%s: (max_dim, num_div, mfv) = (%d, %d, %d)\n',...
        shapes{i}, cparams.max_dim, cparams.num_divs, cparams.max_filtr_val);

    samples_map = containers.Map;
    for k = 1:num_samples

        [stream, complex_info] = complex_factory(shapes{i}, cparams);
        low_true = reduce_stream(stream, 'testing', as_dense);
        m = length(low_true);

        fprintf('\t\tSample %d/%d\tm = %d\n', k, num_samples, m);

        algos_map = containers.Map;
        for l = 1:length(algos)

            fprintf('\t\t\t%s... ', algos{l});

            D = BoundaryMatrix(stream);
            [OUT, t] = cuda_wrapper(D, algos{l}, low_true, 7);
            
            OUT.m = m;
            OUT.complex_info = complex_info;

            algos_map(algos{l}) = OUT;

            assert(all(OUT.low == low_true), 'Output incorrect!');
            fprintf('\t\tsuccess in %g secs!\n', t);

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


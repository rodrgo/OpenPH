% test_speed_standard_reduction.m
init;

% Explore complexity of Vietoris-Rips complexes
vr_complexes = {'house', 'random_figure_8', ...
                'random_torus', 'sphere_product', ...
                'icosahedron', 'random_trefoil_knot', ...
                'random_gaussian', 'morozov'};

% Algorithms to test
algorithms = {'alpha_beta_parallel'};

% Matrix dense?
as_dense = true;

% Complex parameters
max_dim = 5;
num_divs= 10;
max_filtration_values = [5];

% Number of complexes per parameter
num_samples = 10;

time_init = tic;

for i = 1:length(vr_complexes)

    complex = vr_complexes{i};

    for j = 1:length(max_filtration_values)

        mfv = max_filtration_values(j);
        fprintf('\n\t%s @ (max_dim, num_divs, mfv) = (%d, %d, %d)\n',...
            complex, max_dim, num_divs, mfv);

        for k = 1:num_samples

            stream = example_factory(complex, max_dim, mfv, num_divs);
            T = BoundaryMatrix(stream);
            lows_test = reduce_matrix(T, 'testing');
            fprintf('\t\tSample %d/%d\tm = %d\n', k, num_samples, T.m);

            for l = 1:length(algorithms)

                algo = algorithms{l};
                [lows, t, D] = reduce_stream(stream, algo, as_dense);

                %D.metrics.percentage_unreduced(1:D.metrics.iters)

                fprintf('\t\t\t%s... ', algo);
                assert(all(lows == lows_test), 'Output incorrect!');
                fprintf('\t\tsuccess in %g secs!\n', t);

            end

        end

    end

end

time_total = toc(time_init);
fprintf('All tests finished successfully in %s secs :)\n', num2str(time_total));


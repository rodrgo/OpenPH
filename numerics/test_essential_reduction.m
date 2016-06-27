% test_essential_reduction.m

init;

figure_dir = './figures/';
figure_tag = 'test_algos';

% Plot params
LW = 'LineWidth';
MS = 'MarkerSize';
markers = '+o*.xsd^v><ph';

% Explore complexity of Vietoris-Rips complexes
vr_complexes = {'house', 'random_figure_8', ...
                'random_torus', 'sphere_product', ...
                'icosahedron', 'random_trefoil_knot', ...
                'random_gaussian', 'morozov'};

% Algorithms to test
algorithms = {'std_dense', 'twist_dense', 'alpha_beta_dense', 'rho_dense'};
algorithms = {'std_dense'};

% Homology mode
homology_modes = {'reduced', 'unreduced'};

% Complex parameters
max_dim = 5;
num_divs= 10;
max_filtration_values = [5];

% Number of complexes per parameter
num_samples = 1;

time_init = tic;

for h = 1:length(homology_modes)

    homology_mode = homology_modes{h};
    fprintf('homology_mode = %s\n\n', homology_mode);
    for i = 1:length(vr_complexes)

        complex = vr_complexes{i};

        for j = 1:length(max_filtration_values)

            mfv = max_filtration_values(j);
            fprintf('\n\t%s @ (max_dim, num_divs, mfv) = (%d, %d, %d)\n',...
                complex, max_dim, num_divs, mfv);

            for k = 1:num_samples

                stream = example_factory(complex, max_dim, mfv, num_divs);
                fprintf('\t\tSample %d/%d\tm = %d\n', k, num_samples, stream.getSize());

                for l = 1:length(algorithms)

                    algo = algorithms{l};
                    D = BoundaryMatrix(stream, homology_mode);

                    [~, t_red] = reduce_matrix(D, algo);
                    ph_info = D.get_persistence_info();
                    relevant_ess = ph_info{2};

                    D = BoundaryMatrix(stream, homology_mode);
                    [selected_ess, t_ess] = essential_red(D);

                    isSubset = all(ismember(relevant_ess, selected_ess));
                    assert(isSubset, 'relevant_ess are not in selected_ess');

                    % Compute precision and recall 
                    tp = intersect(relevant_ess, selected_ess);
                    fp = setdiff(selected_ess, relevant_ess);
                    fn = setdiff(relevant_ess, tp);

                    num_tp = length(tp);
                    num_fp = length(fp);
                    num_fn = length(fn);

                    precision = num_tp/(num_tp + num_fp); 
                    recall = num_tp/(num_tp + num_fn);

                    fprintf('\t\t\t%s:\n', algo);
                    fprintf('\t\t\t\t\t t_ess/t_red=%g\n', t_ess/t_red);
                    fprintf('\t\t\t\t\t tp=%d, fp=%d, fn=%d\n', num_tp, num_fp, num_fn);
                    fprintf('\t\t\t\t\t precision=%g, recall=%g\n', precision, recall);

                end

            end

        end

    end

end

time_total = toc(time_init);
fprintf('All tests finished successfully in %s secs :)\n', num2str(time_total));


% alpha_beta_curves.m
% Alpha-Beta curves on Vietoris-Rips complexes
% Note: make sure that you give matlab enough heap space to work with
init;

figure_dir = './figures/';
figure_tag = 'alpha_beta_mask';

% Create Vietoris-Rips complexes
% We do not include random_torus
vr_complexes = {'house', 'random_figure_8', ...
                'sphere_product', ...
                'icosahedron', 'random_trefoil_knot', ...
                'random_gaussian'};

max_dimension = 200;
max_filtration_value = 7;
num_divisions = 20;

for i = 1:length(vr_complexes)
    tic;

    % create example
    complex_name = vr_complexes{i};
    fprintf('Working on %s...\t', vr_complexes{i});
    [stream, title_str] = example_factory(complex_name, max_dimension, max_filtration_value, num_divisions);

    D = BoundaryMatrix(stream);
    D.create_alpha();
    D.create_beta();
    masks = {D.get_alpha_mask(), D.get_beta_mask()};

    file_name = strcat(complex_name, '_', figure_tag, '.eps');
    file_path = strcat(figure_dir, file_name);
    plot_matrix(D.matrix, title_str, file_path, masks);

    fprintf('done in %g sec!\n', toc);
end


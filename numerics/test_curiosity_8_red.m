% test_curiosity_8_red.m

init;

figure_dir = './figures/';
figure_tag = 'c8';

% Plot params
LW = 'LineWidth';
MS = 'MarkerSize';
markers = '+o*xsd^v><ph.';

% Explore complexity of Vietoris-Rips complexes
% Do not include house, random_torus, morozov
vr_complexes = {'house', 'random_figure_8', ...
                'sphere_product', ...
                'icosahedron', 'random_trefoil_knot', ...
                'random_gaussian'};

%vr_complexes = {'house'};

% First we plot the boundary matrices with the curiosity_8 mask

max_dim = 5;
num_divs = 5;
mfv = 10;

for i = 1:length(vr_complexes)
    complex = vr_complexes{i};
    complex_label = strrep(complex, '_', '\_');
    fprintf('%s mask...\n', complex);

    stream = example_factory(complex, max_dim, mfv, num_divs);

    % Our way to create the Matlab's boundary matrix
    D = ReductionC8(stream, 'unreduced');
    matrix = D.matrix;

    alpha_mask = D.get_alpha_mask();
    beta_mask = D.get_beta_mask();
    c8_masks = D.get_curiosity_8_mask();

    masks = {alpha_mask, beta_mask, c8_masks{2}, c8_masks{1}, c8_masks{3}};

    file_path = [figure_dir complex '_' figure_tag '.eps'];
    plot_matrix(matrix, complex_label, file_path, masks);
    print('-depsc', file_path);
    eps_to_pdf(file_path);

end

% Homology mode
homology_mode = 'unreduced';

% Complex parameters
max_dimension = [10];
num_divisions= [15];
max_filtration_values = [11];

params = [];
for i = max_dimension
    for ii = num_divisions
        for iii = max_filtration_values
            row = [i, ii, iii];
            params = [params; row];
        end
    end
end

% Number of complexes per parameter
num_samples = 1;

time_init = tic;

% data
data = [];

for i = 1:length(vr_complexes)

    complex = vr_complexes{i};

    for j = 1:size(params, 1)

        max_dim = params(j, 1);
        num_divs = params(j, 2); 
        mfv = params(j, 3); 
        
        fprintf('\n\t%s @ (max_dim, num_divs, mfv) = (%d, %d, %d)\n',...
            complex, max_dim, num_divs, mfv);

        for k = 1:num_samples

            stream = example_factory(complex, max_dim, mfv, num_divs);
            fprintf('\t\tSample %d/%d\tm = %d\n', k, num_samples, stream.getSize());

            % Get groundtruth
            D = BoundaryMatrix(stream, homology_mode);
            [lows, t] = reduce_matrix(D, 'std');
            num_pos = nnz(D.classes == 1);

            % Measure positives
            D = ReductionC8(stream, homology_mode);
            D.init();
            alpha_beta_pos = zeros(1, D.m);
            D.alpha_beta_reduce();
            alpha_beta_pos(D.classes == 1) = 1;

            c8_pos = zeros(1, D.m);
            D.curiosity_8_clearing();
            c8_pos(D.classes == 1 & ~alpha_beta_pos) = 1;

            num_alpha_beta_pos = nnz(alpha_beta_pos);
            num_c8_pos = nnz(c8_pos);

            row = [i, D.m, num_pos, num_alpha_beta_pos, num_c8_pos];
            data = [data; row];

        end

    end

end

% Plot data
set(gcf, 'color', [1 1 1]);
set(gca, 'Fontname', 'setTimes', 'Fontsize', 15);

handles = [];
labels = {};
for i = 1:length(vr_complexes)
    % Get indices data
    complex_idx = data(:, 1) == i;
    data_complex = data(complex_idx, :);

    x = data_complex(:, 4)./data_complex(:, 3);
    y = data_complex(:, 5)./data_complex(:, 3);

    % Plot
    complex_label = strrep(vr_complexes{i}, '_', '\_');
    m = unique(data_complex(:, 2));
    assert(length(m) == 1, 'm is not unique');

    labels{end + 1} = [complex_label '; m=' num2str(m)];
    handles(end + 1) = plot(x, y, markers(i), LW, 1, MS, 5);
    hold on;
end

xlabel('\alpha-\beta\_positives/total\_positives');
ylabel('c8\_positives/total\_positives');
legend(handles, labels, 'Location', 'NorthEast');

title_str = ['Proportion of total positives identified'];
title(title_str);
hold off;

file_name = ['positive_ratio_' figure_tag, '.eps'];
file_path = strcat(figure_dir, file_name);
print('-depsc', file_path);
eps_to_pdf(file_path);

% Report experiment's time to completion

time_total = toc(time_init);
fprintf('All tests finished successfully in %s secs :)\n', num2str(time_total));


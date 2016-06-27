% test_vr_ensembles.m

init;

figure_dir = './figures/';
figure_tag = 'vr_ensembles';

% Create Vietoris-Rips complexes
% Excluding random_torus due to memory overflows at low mfv
vr_complexes = {'house', 'random_figure_8', ...
                'sphere_product', 'icosahedron', ...
                'random_trefoil_knot', 'random_gaussian'};

vr_complex_labels = cell(size(vr_complexes));
for i = 1:length(vr_complexes)
  vr_complex_labels{i} = strrep(vr_complexes{i}, '_', '\_');
end

% Plot params
LW = 'LineWidth';
MS = 'MarkerSize';
vr_markers = '+o*.xsd^v><ph';

% VR ensemble parameters

response_labels = {'time', 'm', 'sparsity'};
y_labels = {'time (ms)', 'm', 'sparsity: nnz/(m(m-1)/2)'};

% max_dimension and num_division_list are fixed
max_dim = 200;
num_div = 20;

max_dimensions = [max_dim];
max_filtration_values = 1:15;
num_division_list = [num_div];

parameter_labels = {'max_dim', 'max_dist', 'num_divisions'};
parameters = {max_dimensions, max_filtration_values, num_division_list}; 

% Get data
data = [];
for i = 1:length(max_dimensions)
    max_dim = max_dimensions(i);

    for j = 1:length(max_filtration_values)
        mfv = max_filtration_values(j);

        for k = 1:length(num_division_list)
            num_div = num_division_list(k);

            for h = 1:length(vr_complexes)
                complex_name = vr_complexes{h};

                t0 = tic;
                stream = example_factory(complex_name, max_dim, mfv, num_div);
                D = BoundaryMatrix(stream, 'unreduced');
                t = 1000*toc(t0);

                params = [max_dim, mfv, num_div];
		% Sparsity is defined as nnz(D)/(m^2 - (m*(m+1)/2)) = nnz(D)/(m(m-1)/2)
                response = [t, D.m, nnz(D.matrix)/(D.m*(D.m - 1)/2)];
                row = [params, response, h];
                data = [data; row];
                fprintf('%s\t%d\t%d\t%d\n', complex_name, max_dim, mfv, num_div);

            end
        end
    end
end

% Plot data
% Keep max_dimension and num_division_list fixed
idx_0 = (data(:, 1) == max_dim) & (data(:, 3) == num_div);

% x parameter is max_dist
x_param = 2;

for r = 1:length(response_labels)
    response_label = response_labels{r};

    handles = [];
    set(gcf, 'color', [1 1 1]);
    set(gca, 'Fontname', 'setTimes', 'Fontsize', 15);

    for h = 1:length(vr_complexes)
        complex_name = vr_complexes{h};
        idx = idx_0 & (data(:, end) == h);

        x_data = data(idx, x_param);
        y_data = data(idx, length(params) + r);

        handles(end + 1) = semilogy(x_data, y_data, '-.', LW, 1.5, MS, 10);
        hold on;
    end

    x_label = parameter_labels{x_param};
    y_label = y_labels{r};

    x_label_ = strrep(x_label, '_', '\_');
    y_label_ = strrep(y_label, '_', '\_');

    xlabel(x_label_);
    ylabel(y_label_);

    legend(handles, vr_complex_labels, 'Location', 'NorthEast');
    title_str = [x_label_, ' vs. ', response_label, ': ', ...
        'max\_dim = ',  num2str(max_dim),  ', ',  ...
        'num\_div = ', num2str(num_div)];
    title(title_str);
    hold off;

    % Save file
    file_name = [x_label, '_vs_', response_label, '_', figure_tag, '.eps'];
    file_path = strcat(figure_dir, file_name);
    print('-depsc', file_path);
    eps_to_pdf(file_path);
end



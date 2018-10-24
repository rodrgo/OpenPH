% test_vr_ensembles.m

init;

figure_dir = './figures/';
figure_tag = 'vr_ensembles_by_numpoints';

% Create Vietoris-Rips complexes
vr_complexes = {'random_figure_8', 'random_torus', ... 
                'sphere_product', 'random_trefoil_knot', ...
                'random_gaussian'};

vr_complex_labels = cell(size(vr_complexes));
for i = 1:length(vr_complexes)
  vr_complex_labels{i} = strrep(vr_complexes{i}, '_', '\_');
end

% Plot params
LW = 'LineWidth';
MS = 'MarkerSize';
vr_markers = '+o*.xsd^v><ph';

% max_dimension and num_division_list are fixed
max_dim = 200;
num_div = 20;
mfv = 15;

% Number of points

num_points_vect = 3:15;

% Get data
data = [];
for h = 1:length(vr_complexes)
    complex_name = vr_complexes{h};
    for i = 1:length(num_points_vect)
        num_points = num_points_vect(i);

        t0 = tic;
        stream = example_factory(complex_name, max_dim, mfv, num_div, num_points);
        D = BoundaryMatrix(stream);
        t = 1000*toc(t0);

        row = [h, num_points, max_dim, num_div, mfv, t, D.m, nnz(D.matrix)];
        data = [data; row];

        fprintf('%s\t%d\n', complex_name, num_points);
    end
end

%% Plot m vs time
handles = [];
labels = {};
set(gcf, 'color', [1 1 1]);
set(gca, 'Fontname', 'setTimes', 'Fontsize', 15);

for h = 1:length(vr_complexes)
    labels{end + 1} = [vr_complex_labels{h} ' params'];
    complex_data = data(data(:, 1) == h, :);
    x = complex_data(:, 7);
    y = complex_data(:, 6);
    handles(end + 1) = semilogy(x, y, '-x', LW, 1.5, MS, 10);
    hold on;
end

xlabel('m'); ylabel('time (ms)');
legend(handles, labels, 'Location', 'SouthEast');
title(['m vs. time (ms)']);
hold off;

file_name = ['m_vs_time' '_', figure_tag, '.eps'];
file_path = strcat(figure_dir, file_name);
print('-depsc', file_path);
eps_to_pdf(file_path);

%% Plot nnz vs time
handles = [];
labels = {};
set(gcf, 'color', [1 1 1]);
set(gca, 'Fontname', 'setTimes', 'Fontsize', 15);

for h = 1:length(vr_complexes)
    labels{end + 1} = [vr_complex_labels{h} ' params'];
    complex_data = data(data(:, 1) == h, :);
    x = complex_data(:, 8);
    y = complex_data(:, 6);
    handles(end + 1) = semilogy(x, y, '-x', LW, 1.5, MS, 10);
    hold on;
end

xlabel('nnz'); ylabel('time (ms)');
legend(handles, labels, 'Location', 'SouthEast');
title(['nnz vs. time (ms)']);
hold off;

file_name = ['nnz_vs_time' '_', figure_tag, '.eps'];
file_path = strcat(figure_dir, file_name);
print('-depsc', file_path);
eps_to_pdf(file_path);

%% Plot num_points vs m 
handles = [];
labels = {};
set(gcf, 'color', [1 1 1]);
set(gca, 'Fontname', 'setTimes', 'Fontsize', 15);

for h = 1:length(vr_complexes)
    labels{end + 1} = [vr_complex_labels{h} ' params'];
    complex_data = data(data(:, 1) == h, :);
    x = complex_data(:, 2);
    y = complex_data(:, 7);
    handles(end + 1) = semilogy(x, y, '-x', LW, 1.5, MS, 10);
    hold on;
end

xlabel('num\_points'); ylabel('m');
legend(handles, labels, 'Location', 'SouthEast');
title(['num\_points vs. m']);
hold off;

file_name = ['numpoints_vs_m' '_', figure_tag, '.eps'];
file_path = strcat(figure_dir, file_name);
print('-depsc', file_path);
eps_to_pdf(file_path);

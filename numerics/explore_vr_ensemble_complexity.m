% explore_ensemble_complexity.m
init;


figure_dir = './figures/';
figure_tag = 'speed_standard_reduction';

% Create Vietoris-Rips complexes
vr_complexes = {'house', 'random_figure_8', ...
                'random_torus', 'sphere_product', ...
                'icosahedron', 'random_trefoil_knot', ...
                'random_gaussian'};

vr_complex_labels = cell(size(vr_complexes));
for i = 1:length(vr_complexes)
  vr_complex_labels{i} = strrep(vr_complexes{i}, '_', '\_');
end

vr_markers = '+o*.xsd^v><ph';

max_dimension_list = [1, 3, 10];
max_filtration_value_list = [1, 10, 20];
num_divisions_list = [10, 20, 30];

data_time = [];
data_dimension = [];
data_sparsity = [];

for i = 1:length(vr_complexes)
  complex_name = vr_complexes{i};
  fprintf('Working on %s...\t', complex_name);

  for j = 1:length(max_dimension_list)
    max_dimension = max_dimension_list(j);

    for k = 1:length(max_filtration_value_list)
      max_filtration_value = max_filtration_value_list(k);

      for l = 1:length(num_divisions_list)
        num_divisions = num_divisions_list(l);

        row = [i, max_dimension, max_filtration_value, ...
          num_divisions, 0];

        t0 = tic;

        [stream, title_str] = example_factory(complex_name, ...
          max_dimension, max_filtration_value, num_divisions);
        D = BoundaryMatrix(stream, 'plain');

        row(end) = 1000*toc(t0);
        data_time = [data_time; row];

        row(end) = D.m;
        data_dimension = [data_dimension; row];

        row(end) = nnz(D.matrix)/((D.m)^2);
        data_sparsity = [data_sparsity; row];

      end
    end
  end
  fprintf('\n');
end

% Merge datasets

data = {data_time, data_dimension, data_sparsity};
data_labels = {'time', 'm', 'sparsity'};
matrix_labels = {'complex', 'max_dim', 'max_dist', 'num_divisions'};

% Plot results
% Plot (i,j) = data_labels{i} x matrix_labels{j}
LW = 'LineWidth';
MS = 'MarkerSize';

grid_size = length(data)*length(matrix_labels);

%%
% One complex per plot
%%

for k = 1:length(vr_complexes)
  for i = 1:length(data_labels)
    y_label = data_labels{i};
    data_matrix = data{i};
    complex_data = data_matrix(data_matrix(:,1) == k, :);
    for j = 2:length(matrix_labels)
      x_label = matrix_labels{j};

      % Factor data on j
      num_factors = size(complex_data, 2) - 1;
      factors_cols = (1:num_factors ~= j) & (1:num_factors ~= 1);
      factors_j = unique(complex_data(:, factors_cols), 'rows');
      label_keys = matrix_labels(factors_cols);
      factor_labels = {};
      for f = 1:size(factors_j, 1)
        label_values = factors_j(f, :);
        label = '';
        for ff = 1:length(label_keys)
          label = [label label_keys{ff} ':' num2str(label_values(ff)) ' '];
        end
        factor_labels{f} = strrep(label, '_', '\_');
      end

      handles = [];
      set(gcf, 'color', [1 1 1]);
      set(gca, 'Fontname', 'setTimes', 'Fontsize', 15);
      for l = 1:size(factors_j, 1)
        [~, factor_ind] = ismember(complex_data(:, factors_cols), factors_j(l, :), 'rows');
        x = complex_data(factor_ind == 1, j);
        y = complex_data(factor_ind == 1, end);
        handles(end + 1) = semilogy(x, y, [vr_markers(l) '-'], LW, 1.5, MS, 10);
        hold on;
      end

      x_label_f = strrep(x_label, '_', '\_');
      y_label_f = strrep(y_label, '_', '\_');
      xlabel(x_label_f);
      ylabel(y_label_f);
      legend(handles, factor_labels, 'Location', 'NorthEast');
      title([vr_complex_labels{k} ': ' x_label_f, ' vs. ', y_label_f]);
      hold off;

      % Save file
      file_name = ['ensemble_complexity_' vr_complexes{k} '_' x_label '_vs_' y_label '.eps'];
      file_path = strcat(figure_dir, file_name);
      print('-depsc', file_path);
    end
  end
end

%%
% Plot all data

if 1 == 0
  for i = 1:length(data_labels)
    y_label = data_labels{i};
    data_matrix = data{i};
    for j = 2:length(matrix_labels)
      x_label = matrix_labels{j};
      handles = [];
      set(gcf, 'color', [1 1 1]);
      set(gca, 'Fontname', 'setTimes', 'Fontsize', 15);
      for k = 1:length(vr_complexes)
        complex_data = data_matrix(data_matrix(:,1) == k, :);
        x = complex_data(:, j);
        y = complex_data(:, end);
        handles(end + 1) = semilogy(x, y, vr_markers(k), LW, 1.5, MS, 10);
        hold on;
      end
      x_label_f = strrep(x_label, '_', '\_');
      y_label_f = strrep(y_label, '_', '\_');
      xlabel(x_label_f);
      ylabel(y_label_f);
      legend(handles, vr_complex_labels, 'Location', 'NorthEast');
      title([x_label_f, ' vs. ', y_label_f]);
      hold off;
      % Save file
      file_name = ['ensemble_complexity_' 'all_' x_label '_vs_' y_label '.eps'];
      file_path = strcat(figure_dir, file_name);
      print('-depsc', file_path);
    end
  end
end

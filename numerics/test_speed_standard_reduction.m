% test_speed_standard_reduction.m
init;

figure_dir = './figures/';
figure_tag = 'speed_std_red';

% Explore complexity of Vietoris-Rips complexes
vr_complexes = {'random_torus', 'random_figure_8', ...
  'random_gaussian', 'random_trefoil_knot'};

% Fixed complex parameters
max_dimension = 5;
num_divisions = 10;

% Variable complex parameters
max_filtration_values = 5:15;

% Number of complexes per parameter
num_samples = 10;
num_samples = 1;

% Time

avg_m = zeros(1, length(max_filtration_values));
for i = 1:length(vr_complexes)
  % create example
  complex_name = vr_complexes{i};

  times_complex = zeros(3, length(avg_m));
  times_row_labels = {'std_red_sparse', 'std_red_sparse_opt', 'std_red_dense_opt'};

  for j = 1:length(max_filtration_values)
    max_filtration_value = max_filtration_values(j);
    time_sr = zeros(1, num_samples);
    time_sparse = zeros(1, num_samples);
    time_dense = zeros(1, num_samples);
    m_complex = zeros(1, num_samples);
    fprintf('%s filtration value %s\n',...
      vr_complexes{i}, num2str(max_filtration_value));

    for k = 1:num_samples

      fprintf('\tSample %s/%s\n',...
        num2str(k), num2str(num_samples));

      [stream, title_str] = example_factory(complex_name, max_dimension, max_filtration_value, num_divisions);

      D = BoundaryMatrix(stream, 'plain');

      m_complex(k) = D.m;

      [lows, t] = D.standard_reduction();
      t_sr(k) = 1000*t;

      [lows, t] = D.standard_reduction_sparse_opt();
      t_sparse(k) = 1000*t;

      [lows, t] = D.standard_reduction_dense_opt();
      t_dense(k) = 1000*t;

    end
    avg_m(j) = mean(m_complex);
    times_complex(1, j) = mean(t_sr);
    times_complex(2, j) = mean(t_sparse);
    times_complex(3, j) = mean(t_dense);
  end
  % Plot results
  handles = [];
  set(gcf, 'color', [1 1 1]);
  set(gca, 'Fontname', 'setTimes', 'Fontsize', 15);
  for ii = 1:size(times_complex, 1)
    x = floor(avg_m);
    y = times_complex(ii, :);
    handles(end + 1) = semilogy(x, y, [times_row_labels{ii} '-'], LW, 1.5, MS, 10);
    hold on;
  end
  hold off;
  file_name = strcat(complex_name, '_', figure_tag, '.eps');
  file_path = strcat(figure_dir, file_name);
  print('-depsc', file_path);
  fprintf('done in %g sec!\n', toc);
end


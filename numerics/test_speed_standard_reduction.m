% test_speed_standard_reduction.m
init;

figure_dir = './figures/';
figure_tag = 'speed_std_red';

% Plot params
LW = 'LineWidth';
MS = 'MarkerSize';
markers = '+o*.xsd^v><ph';

% Explore complexity of Vietoris-Rips complexes
vr_complexes = {'random_torus', 'random_figure_8', ...
  'random_gaussian', 'random_trefoil_knot'};

% Fixed complex parameters
max_dimension = 5;
num_divisions = 10;

% Variable complex parameters
max_filtration_values = 1:5;

% Number of complexes per parameter
num_samples = 10;

% Time

avg_m = zeros(1, length(max_filtration_values));
for i = 1:length(vr_complexes)
  % create example
  complex_name = vr_complexes{i};

  times_complex = zeros(4, length(avg_m));
  times_row_labels = {'std_sparse', 'std_sparse_opt', 'std_dense_opt', 'twist_sparse'};
  for l = 1:length(times_row_labels)
    times_row_labels{l} = strrep(times_row_labels{l}, '_', '\_');
  end

  for j = 1:length(max_filtration_values)
    max_filtration_value = max_filtration_values(j);

    time_sr_sparse = zeros(1, num_samples);
    time_sr_sparse_opt = zeros(1, num_samples);
    time_sr_dense_opt = zeros(1, num_samples);
    time_twist_sparse = zeros(1, num_samples);

    m_complex = zeros(1, num_samples);
    fprintf('%s filtration value %s\n',...
      vr_complexes{i}, num2str(max_filtration_value));

    for k = 1:num_samples

      fprintf('\tSample %s/%s\n',...
        num2str(k), num2str(num_samples));

      [stream, title_str] = example_factory(complex_name, max_dimension, max_filtration_value, num_divisions);

      D = BoundaryMatrix(stream, 'plain');

      m_complex(k) = D.m;

      [lows, t] = D.standard_reduction_sparse();
      time_sr_sparse(k) = 1000*t;

      [lows_sso, t] = D.standard_reduction_sparse_opt();
      time_sr_sparse_opt(k) = 1000*t;

      [lows_sdo, t] = D.standard_reduction_dense_opt();
      time_sr_dense_opt(k) = 1000*t;

      [lows_ts, t] = D.twist_reduction_sparse();
      time_twist_sparse(k) = 1000*t;

      fprintf('\t std_sparse_opt: %d\t std_dense_opt: %d\t twist_sparse: %d\n', ...
        sum(lows ~= lows_sso), sum(lows ~= lows_sdo), sum(lows ~= lows_ts));

    end
    avg_m(j) = mean(m_complex);
    times_complex(1, j) = mean(time_sr_sparse);
    times_complex(2, j) = mean(time_sr_sparse_opt);
    times_complex(3, j) = mean(time_sr_dense_opt);
    times_complex(4, j) = mean(time_twist_sparse);
  end

  % Plot results
  handles = [];
  set(gcf, 'color', [1 1 1]);
  set(gca, 'Fontname', 'setTimes', 'Fontsize', 15);
  for ii = 1:size(times_complex, 1)
    x = floor(avg_m);
    y = times_complex(ii, :);
    handles(end + 1) = semilogy(x, y, [markers(ii) '-'], LW, 1.5, MS, 10);
    hold on;
  end
  xlabel('Average m');
  ylabel('log(time) ms');
  legend(handles, times_row_labels, 'Location', 'NorthWest');
  title([strrep(complex_name, '_', '\_') ': ' 'reduction average time']);
  hold off;

  file_name = strcat(complex_name, '_', figure_tag, '.eps');
  file_path = strcat(figure_dir, file_name);
  print('-depsc', file_path);
  eps_to_pdf(file_path);
  fprintf('done in %g sec!\n', toc);
end


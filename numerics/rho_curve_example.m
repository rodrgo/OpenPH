% alpha_beta_curves.m
% Alpha-Beta curves on Vietoris-Rips complexes
% Note: make sure that you give matlab enough heap space to work with
init;

figure_dir = './figures/';
figure_tag = 'rho_curve_example';

% Create Vietoris-Rips complexes
% We do not include random_torus
complexes = {'morozov', 'house', 'random_figure_8', ...
                'sphere_product', ...
                'icosahedron', 'random_trefoil_knot', ...
                'random_gaussian'};

max_dimension = 3;
max_filtration_value = 4;
num_divisions = 10;

rho_curves = {};

for i = 1:length(complexes)
  tic;

  % create example
  complex_name = complexes{i};
  fprintf('Working on %s...\t', complexes{i});
  [stream, title_str] = example_factory(complex_name, max_dimension, max_filtration_value, num_divisions);

  D = BoundaryMatrix(stream, 'unreduced');

  D.create_rho();
  rho_curves{end + 1} = D.rho;

%  file_name = strcat(complex_name, '_', figure_tag, '.eps');
%  file_path = strcat(figure_dir, file_name);
%  plot_matrix(D.matrix, title_str, file_path, masks);

  fprintf('done in %g sec!\n', toc);
end

LW = 'LineWidth';
MS = 'MarkerSize';
set(gcf, 'color', [1 1 1]);
set(gca, 'Fontname', 'setTimes', 'Fontsize', 15);

handles = [];
rho_lengths = zeros(size(rho_curves));
for i = 1:length(rho_curves)
  rho_curve = rho_curves{i};
  rho_lengths(i) = length(rho_curve);
  handles(end + 1) = plot(1:length(rho_curve), rho_curve, LW, 1.5, MS, 10);
  hold on;
  zero_idx = find(rho_curve == 0);
  if length(zero_idx) > 1
    fprintf('The rho curve of %s has %d zeros!\n', complexes{i}, length(zero_idx));
    plot(zero_idx, zeros(size(zero_idx)), '*k', MS, 20);
    hold on;
  end
end
plot(1:max(rho_lengths), zeros(1, max(rho_lengths)), '--k');
for i = 1:length(complexes)
  complexes{i} = strrep(complexes{i}, '_', '\_');
end
xlabel('j in [m]');
ylabel('{rho}(j)');
legend(handles, complexes, 'Location', 'SouthWest');
title('rho curves');
hold off;

file_name = strcat(figure_tag, '.eps');
file_path = strcat(figure_dir, file_name);
print('-depsc', file_path);
eps_to_pdf(file_path);


% test_essential_reduction.m

init;

figure_dir = './figures/';
figure_tag = 'essential_estimation';

% Plot params
LW = 'LineWidth';
MS = 'MarkerSize';
markers = '+o*xsd^v><ph.';

% Explore complexity of Vietoris-Rips complexes
% Do not include house, random_torus, morozov
vr_complexes = {'random_figure_8', ...
                'sphere_product', ...
                'icosahedron', 'random_trefoil_knot', ...
                'random_gaussian'};

% Algorithms to test
algorithms = {'std', 'twist', 'alpha_beta_std', 'rho_std'};
algorithm = 'std';

% Homology mode
homology_mode = 'unreduced';

% Complex parameters
max_dimension = [3, 5, 10];
num_divisions= [5, 10, 15];
max_filtration_values = [3, 7, 11];

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
            [~, t_red, D] = reduce_stream(stream, homology_mode, algorithm, true);
            ph_info = D.get_persistence_info();
            relevant_ess = ph_info{2};

            num_rel_ess = length(relevant_ess);

            % Estimate essentials
            D = BoundaryMatrix(stream, homology_mode);
            [selected_ess, t_ess] = essential_std(D);

            num_sel_ess = length(selected_ess);

            % Assert relevant_ess \subset selected_ess
            isSubset = all(ismember(relevant_ess, selected_ess));
            assert(isSubset, 'relevant_ess are not in selected_ess');

            % Get data
            tp = intersect(relevant_ess, selected_ess);
            fp = setdiff(selected_ess, relevant_ess);
            fn = setdiff(relevant_ess, tp);

            num_tp = length(tp);
            num_fp = length(fp);
            num_fn = length(fn);

            row = [i, D.m, num_rel_ess, num_sel_ess, num_tp, num_fp, num_fn];
            data = [data; row];

            if false
                precision = num_tp/(num_tp + num_fp); 
                recall = num_tp/(num_tp + num_fn);

                fprintf('\t\t\t%s:\n', algorithm);
                fprintf('\t\t\t\t\t t_ess/t_red=%g\n', t_ess/t_red);
                fprintf('\t\t\t\t\t tp=%d, fp=%d, fn=%d\n', num_tp, num_fp, num_fn);
                fprintf('\t\t\t\t\t precision=%g, recall=%g\n', precision, recall);
            end

        end

    end

end

% Plot data
set(gcf, 'color', [1 1 1]);
set(gca, 'Fontname', 'setTimes', 'Fontsize', 15);

handles = [];
for i = 1:length(vr_complexes)
    % Get indices data
    complex_idx = data(:, 1) == i;

    % Compute precision
    ms = data(complex_idx, 2);
    num_rel_ess = data(complex_idx, 3);
    num_sel_ess = data(complex_idx, 4);
    tp = data(complex_idx, 5);
    fp = data(complex_idx, 6);
    precision = tp./(tp + fp);
    precision(num_rel_ess == 0) = 1;

    x = (tp + fp)./ms; 

    % Plot
    handles(end + 1) = plot(x, precision, markers(i), LW, 1, MS, 5);
    hold on;
end

vr_complex_labels = cell(size(vr_complexes));
for i = 1:length(vr_complex_labels)
    vr_complex_labels{i} = strrep(vr_complexes{i}, '_', '\_');
end

xlabel('(TP + FP)/m');
ylabel('Precision: TP/(TP + FP)');
legend(handles, vr_complex_labels, 'Location', 'SouthEast');

title_str = ['Estimation of essential classes in O(m)'];
title(title_str);
hold off;

file_name = [figure_tag, '.eps'];
file_path = strcat(figure_dir, file_name);
print('-depsc', file_path);
eps_to_pdf(file_path);

% Report experiment's time to completion

time_total = toc(time_init);
fprintf('All tests finished successfully in %s secs :)\n', num2str(time_total));


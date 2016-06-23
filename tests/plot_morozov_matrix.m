% TEST_MATLAB_BOUNDARY_MATRIX.M
% The results of these tests show that we can create the boundary matrix
% either by 
% 1. 
% 	streams.utility.StreamUtility.getCompressedBoundaryMatrix(stream);
%	D = BoundaryMatrix(stream)
%	matrix = D.matrix
%
% 2. 
% 	formal_sum = streams.utility.StreamUtility.createBoundaryMatrixAsDoubleSum(stream);
% 	matrix_converter = api.Plex4.createHomMatrixConverter(stream, stream);
% 	matrix = to_sparse_matlab_matrix(formal_sum, matrix_converter);
%
init;

figure_dir = './figures/';
tag = 'morozov';
complex = 'morozov'

order_complexity = 7;

for i = 5:order_complexity

    	stream = examples.MorozovCubicTimeExample.getMorozovCubicTimeExample(order_complexity);
	complex_name = [complex '_order_' num2str(i)];
	complex_name_label = strrep(complex_name, '_', '\_');

	% Our way to create the Matlab's boundary matrix
	D = BoundaryMatrix(stream, 'reduced');
	matrix_tda = D.matrix;

	% Spy boundary matrices
	figure;
	set(gcf, 'color', [1 1 1]);
	set(gca, 'xtick', [], 'ytick', [], 'XTickLabel', '', 'YTickLabel', '');

	spy_tda(matrix_tda);
	[m, n] = size(matrix_tda);
	xlabel(['(m, n) = (' num2str(m) ',' num2str(n) ')']);
	title([complex_name_label ': tda']);

	% Save plot
	file_path = [figure_dir complex_name '_' tag '.eps'];
	print('-depsc', file_path);
	eps_to_pdf(file_path);

	fprintf('done with %s!\n', complex_name);
end


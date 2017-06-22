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

import edu.stanford.math.plex4.*;

figure_dir = './figures/';
tag = 'bd_mat_comp';

% Create Vietoris-Rips complexes
vr_complexes = {'house', 'random_figure_8', ...
                'random_torus', 'sphere_product', ...
                'icosahedron', 'random_trefoil_knot', ...
                'morozov', 'random_gaussian'};

max_dimension = 5;
max_filtration_value = 5;
num_divisions = 10;

for i = 1:length(vr_complexes)
	complex_name = vr_complexes{i};
	complex_name_label = strrep(complex_name, '_', '\_');

	stream = example_factory(complex_name, ...
	max_dimension, max_filtration_value, num_divisions);

	% Our way to create the Matlab's boundary matrix
	D = BoundaryMatrix(stream, 'unreduced');
	matrix_tda = D.matrix;

	% Javaplex's way to create Matlab's boundary matrix
	% INPUTS:
	% 	formal_sum: an object of type DoubleSparseFormalSum<ObjectObjectPair<M, N>>
	% 	matrix_converter: the object of type DoubleMatrixConverter<M, N>
	%
	%formal_sum = streams.utility.StreamUtility.createBoundaryMatrixAsIntSum(stream);
	formal_sum = streams.utility.StreamUtility.createBoundaryMatrixAsDoubleSum(stream);
	matrix_converter = api.Plex4.createHomMatrixConverter(stream, stream);
	matrix_javaplex = to_sparse_matlab_matrix(formal_sum, matrix_converter);
	matrix_javaplex = matrix_javaplex';

	% Spy boundary matrices
	figure;
	set(gcf, 'color', [1 1 1]);
	set(gca, 'xtick', [], 'ytick', [], 'XTickLabel', '', 'YTickLabel', '');

	subplot(1,2,1);
	spy_tda(matrix_tda);
	[m, n] = size(matrix_tda);
	xlabel(['(m, n) = (' num2str(m) ',' num2str(n) ')']);
	title([complex_name_label ': tda']);

	subplot(1,2,2);
	spy_tda(matrix_javaplex);
	[m, n] = size(matrix_javaplex);
	xlabel(['(m, n) = (' num2str(m) ',' num2str(n) ')']);
	title([complex_name_label ': javaplex (transpose)']);

	fprintf('norm diff %s: %5.6f\n', complex_name, norm(abs(matrix_tda) - abs(matrix_javaplex), 1));
	% Save plot
	file_path = [figure_dir complex_name '_' tag '.eps'];
	print('-depsc', file_path);
	eps_to_pdf(file_path);

	fprintf('done with %s sec!\n', complex_name);
end


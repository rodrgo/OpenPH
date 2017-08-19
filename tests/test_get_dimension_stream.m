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

%Create Vietoris-Rips complexes
vr_complexes = {'house', 'random_figure_8', ...
                'random_torus', 'sphere_product', ...
                'icosahedron', 'random_trefoil_knot', ...
                'morozov', 'random_gaussian'};

max_dimension = 5;
max_filtration_value = 5;
num_divisions = 10;

for i = 1:length(vr_complexes)
    complex_name = vr_complexes{i};

    stream = example_factory(complex_name, max_dimension, max_filtration_value, num_divisions);

    % Our way to create the Matlab's boundary matrix
    D = BoundaryMatrix(stream);
    D.as_dense();

    D.simplex_dimensions

    % Get dimensions
end


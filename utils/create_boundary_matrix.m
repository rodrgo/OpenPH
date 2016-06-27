% CREATE_BOUNDARY_MATRIX.M

function D = create_boundary_matrix(stream)

import edu.stanford.math.plex4.*;

m = stream.getSize() + 1; %take dummy simplex index = 1 into account

ccs = streams.utility.StreamUtility.getCompressedBoundaryMatrix(stream);

rows_integer = ccs.get(0).toArray();
cols_integer = ccs.get(1).toArray();

rows = zeros(size(rows_integer));
cols = zeros(size(cols_integer));
vals = ones(size(cols_integer));

for i = 1:length(rows_integer)
    rows(i) = rows_integer(i);
    cols(i) = cols_integer(i);
end

matrix = sparse(rows, cols, vals, m, m);

nnz_count = sum(matrix);
centers = 0.5:1:(max(nnz_count) - 0.5);
element_count = hist(nnz_count, centers);

D.matrix = matrix;
D.element_count = @(i) element_count(i);
D.size = size(D.matrix);

end

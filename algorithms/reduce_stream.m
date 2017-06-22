% reduce_stream.m
%   Wraps BoundaryMatrix creation and reduction

function [lows, t, D] = reduce_stream(stream, homology_mode, algorithm, as_dense)

    D = BoundaryMatrixFactory(stream, homology_mode, algorithm, as_dense);
    [lows, t] = reduce_matrix(D, algorithm);

end

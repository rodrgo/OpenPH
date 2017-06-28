% reduce_stream.m
%   Wraps BoundaryMatrix creation and reduction

function [lows, t, D] = reduce_stream(stream, algorithm, as_dense)

    D = BoundaryMatrixFactory(stream, algorithm, as_dense);
    [lows, t] = reduce_matrix(D, algorithm);

end

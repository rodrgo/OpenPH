% reduce_stream.m
%   Wraps BoundaryMatrix creation and reduction

function [lows, t, D] = reduce_stream(stream, algorithm, as_dense, true_lowstar)

    if nargin > 3
        D = BoundaryMatrixFactory(stream, algorithm, as_dense, true_lowstar);
    else
        D = BoundaryMatrixFactory(stream, algorithm, as_dense);
    end
    
    [lows, t] = reduce_matrix(D, algorithm);

end

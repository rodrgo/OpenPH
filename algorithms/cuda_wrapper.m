% cuda_wrapper.M
% Standard reduction algorithm

function [lows, t] = cuda_wrapper(D, algorithm)
    % Matrix
    m = D.m;
    [r, c, v] = find(D.matrix);
    % Run algorithm
    t0 = tic;
    [low resRecord timeRecord] = ph(algorithm, int32(r), int32(c), int32(m));
    t = toc(t0);
    % Double lows
    lows = double(low);
end

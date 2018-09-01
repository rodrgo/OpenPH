% cuda_wrapper.M
% Standard reduction algorithm

function [lows, t] = cuda_wrapper(D, algorithm)
    % Matrix
    m = D.m;
    [r, c, v] = find(D.matrix);
    % with metrics
    if D.has_lowstar_tracking
        lows_true = D.true_lowstar;
    else
        lows_true = zeros(1, m);
    end
    % Run algorithm
    t0 = tic;
    [low resRecord timeRecord] = ph(algorithm, int32(r), int32(c), int32(m), int32(5), int32(lows_true));
    t = toc(t0);
    % Double lows
    lows = double(low);
end

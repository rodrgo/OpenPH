% cuda_wrapper.M
% Standard reduction algorithm

function [low, t, ess, err_linf, err_lone, err_redu, err_ess, time_track] = cuda_wrapper(D, algorithm)
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
    [o1 o2 o3 o4 o5 o6 o7] = ph(algorithm, int32(r), int32(c), int32(m), int32(5), int32(lows_true));
    t = toc(t0);
    % outputs
    low         = double(o1);
    ess         = double(o2);
    err_linf    = double(o3);
    err_lone    = double(o4);
    err_redu    = double(o5);
    err_ess     = double(o6);
    time_track  = double(o7);
    if false
        display(low)
        display(ess)
        display(err_linf)
        display(err_lone)
        display(err_redu)
        display(err_ess)
        display(time_track)
    end
end

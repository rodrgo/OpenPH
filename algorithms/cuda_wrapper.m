% cuda_wrapper.M

function [OUT, t] = cuda_wrapper(D, algo, low_true, col_width)
    % Matrix
    m = D.m;
    [r, c, ~] = find(D.matrix);
    % Run algorithm
    t0 = tic;
    [o1 o2 o3 o4 o5 o6 o7 o8] = ph(algo, int32(r), int32(c), int32(m), int32(col_width), int32(low_true));
    t = toc(t0);
    % outputs
    OUT             = [];
    OUT.low         = double(o1);
    OUT.ess         = double(o2);
    OUT.err_linf    = double(o3);
    OUT.err_lone    = double(o4);
    OUT.err_redu    = double(o5);
    OUT.err_ess     = double(o6);
    OUT.time_track  = double(o7);
    OUT.num_iters   = double(o8);
end

% ALPHA_BETA_STD.M
% 

function [lows, t] = alpha_beta_std(D)
    t0 = tic;
    % Initialise persistence vectors
    D.init();
    D.record_iteration();
    % Mark as lowstar all columns with alpha=beta
    D.alpha_beta_reduce();
    % Reduce remaining columns
    unreduced_cols = D.get_unreduced_cols();
    for j = unreduced_cols
        D.reduce_col(j);
        D.record_iteration();
    end
    % Extract lows
    lows = D.low;
    t = toc(t0);
end

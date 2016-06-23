% ALPHA_BETA_RED.M
% 

function [lows, t] = alpha_beta_red(D)
    t0 = tic;
    % Initialise persistence vectors
    D.init();
    % Mark as lowstar all columns with alpha=beta
    D.alpha_beta_reduce();
    % Reduce remaining columns
    unreduced_cols = D.get_unreduced_cols();
    for j = unreduced_cols
        D.reduce_col(j);
    end
    % Extract lows
    lows = D.low;
    t = toc(t0);
end

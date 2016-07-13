% ALPHA_BETA_TWIST.M
% 

function [lows, t] = alpha_beta_twist(D)
    t0 = tic;
    % Initialise persistence vectors
    D.init();

    % Mark as lowstar all columns with alpha=beta
    D.alpha_beta_reduce();

    % Reduce remaining columns
    twist_cols = D.get_twist_cols_unreduced();
    for j = twist_cols
        D.reduce_col_twist(j);
    end

    % Extract lows
    lows = D.low;
    t = toc(t0);
end

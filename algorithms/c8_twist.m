% C8_TWIST.M
%

function [lows, t] = c8_twist(D)
    t0 = tic;

    % Initialise persistence vectors
    D.init();

    % Do alpha/beta reduction
    D.alpha_beta_reduce();

    % Do curiosity_8_clearing
    D.curiosity_8_clearing();

    % Reduce remaining columns
    twist_cols = D.get_twist_cols_unreduced();
    for j = twist_cols
        D.reduce_col_twist(j);
    end

    % Extract lows
    lows = D.low;
    t = toc(t0);
end

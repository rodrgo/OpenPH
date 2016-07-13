% C8_STD.M
%

function [lows, t] = c8_std(D)
    t0 = tic;

    % Initialise persistence vectors
    D.init();

    % Do alpha/beta reduction
    D.alpha_beta_reduce();

    % Do curiosity_8_clearing
    D.curiosity_8_clearing();

    % Reduce remaining columns
    unreduced_cols = D.get_unreduced_cols();
    for j = unreduced_cols
        D.reduce_col(j);
    end

    % Extract lows
    lows = D.low;
    t = toc(t0);
end

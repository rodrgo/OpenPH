% RHO_STD.M
% 

function [lows, t] = rho_std(D)
    t0 = tic;
    % Initialise persistence vectors
    D.init();
    % Start clearing via rho curve
    D.create_rho();
    D.rho_clearing();
    % Do alpha-beta reduction in negative columns
    D.create_alpha();
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

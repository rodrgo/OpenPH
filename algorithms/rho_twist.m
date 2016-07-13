% RHO_TWIST.M
% 

function [lows, t] = rho_twist(D)
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
    twist_cols = D.get_twist_cols_unreduced();
    for j = twist_cols
        D.reduce_col_twist(j);
    end

    % Extract lows
    lows = D.low;
    t = toc(t0);
end

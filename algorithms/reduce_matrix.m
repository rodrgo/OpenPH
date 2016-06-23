% REDUCE_MATRIX.M

function [lows, t] = reduce_matrix(D, algorithm)

    % Vanilla algorithm for testing.
    if strcmp(algorithm, 'testing')
        R = D.matrix;
        [lows, t] = std_red_testing(R);
    % Standard algorithm
    elseif strcmp(algorithm, 'std_sparse') 
        [lows, t] = std_red(D);
    elseif strcmp(algorithm, 'std_dense') 
        D.as_dense();
        [lows, t] = std_red(D);
    % Twist algorithm
    elseif strcmp(algorithm, 'twist_sparse') 
        [lows, t] = twist_red(D);
    elseif strcmp(algorithm, 'twist_dense') 
        D.as_dense();
        [lows, t] = twist_red(D);
    % Alpha-beta algorithm
    elseif strcmp(algorithm, 'alpha_beta_sparse') 
        [lows, t] = alpha_beta_red(D);
    elseif strcmp(algorithm, 'alpha_beta_dense') 
        D.as_dense();
        [lows, t] = alpha_beta_red(D);
    % Rho-curve algorithm
    elseif strcmp(algorithm, 'rho_sparse') 
        [lows, t] = rho_red(D);
    elseif strcmp(algorithm, 'rho_dense') 
        D.as_dense();
        [lows, t] = rho_red(D);
    else
        assert(false, 'Algorithm not identified\n');
    end

end

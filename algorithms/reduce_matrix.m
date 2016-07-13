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
    %% Alpha-beta algorithm
    % Standard reduction
    elseif strcmp(algorithm, 'alpha_beta_std_sparse') 
        [lows, t] = alpha_beta_std(D);
    elseif strcmp(algorithm, 'alpha_beta_std_dense') 
        D.as_dense();
        [lows, t] = alpha_beta_std(D);
    % Twist reduction
    elseif strcmp(algorithm, 'alpha_beta_twist_sparse') 
        [lows, t] = alpha_beta_twist(D);
    elseif strcmp(algorithm, 'alpha_beta_twist_dense') 
        D.as_dense();
        [lows, t] = alpha_beta_twist(D);
    %% Rho-curve algorithm
    % Standard reduction
    elseif strcmp(algorithm, 'rho_std_sparse') 
        [lows, t] = rho_std(D);
    elseif strcmp(algorithm, 'rho_std_dense') 
        D.as_dense();
        [lows, t] = rho_std(D);
    %Twist algorithm
    elseif strcmp(algorithm, 'rho_twist_sparse') 
        [lows, t] = rho_twist(D);
    elseif strcmp(algorithm, 'rho_twist_dense') 
        D.as_dense();
        [lows, t] = rho_twist(D);
    % C8 algorithm
    % Standard reduction
    elseif strcmp(algorithm, 'c8_std_sparse') 
        [lows, t] = c8_std(D);
    elseif strcmp(algorithm, 'c8_std_dense') 
        D.as_dense();
        [lows, t] = c8_std(D);
    % C8 algorithm
    % Twist reduction
    elseif strcmp(algorithm, 'c8_twist_sparse') 
        [lows, t] = c8_twist(D);
    elseif strcmp(algorithm, 'c8_twist_dense') 
        D.as_dense();
        [lows, t] = c8_twist(D);
    else
        assert(false, 'Algorithm not identified\n');
    end
end

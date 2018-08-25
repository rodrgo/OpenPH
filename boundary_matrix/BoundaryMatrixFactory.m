% BoundaryMatrixFactory
%   Creates an instance of BoundaryMatrix with the
%   corresponding methods to reduce with each
%   algorithm.

function D = BoundaryMatrixFactory(stream, algorithm, as_dense, true_lowstar)

    % Default algorithms implemented by BoundaryMatrix
    default_algorithms = {'testing', 'std', 'twist', ...
        'alpha_beta_std', 'alpha_beta_twist'};

    if any(ismember(algorithm, default_algorithms))
        D = BoundaryMatrix(stream);
    elseif any(ismember(algorithm, {'rho_std', 'rho_twist'}))
        D = ReductionRho(stream);
    elseif any(ismember(algorithm, {'c8_std', 'c8_twist'}))
        D = ReductionC8(stream);
    elseif any(ismember(algorithm, {'alpha_beta_parallel'}))
        D = ReductionParallel(stream);
    elseif any(ismember(algorithm, {'ph_row'}))
        D = BoundaryMatrix(stream);
    else
        error('Algorithm not identified\n');
    end

    if nargin > 3 
        D.enable_lowstar_tracking(true_lowstar);
    end

    if as_dense
        D.as_dense();
    end

end

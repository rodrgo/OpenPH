% BoundaryMatrixFactory
%   Creates an instance of BoundaryMatrix with the
%   corresponding methods to reduce with each
%   algorithm.

function D = BoundaryMatrixFactory(stream, mode, algorithm, as_dense)

    % Default algorithms implemented by BoundaryMatrix
    default_algorithms = {'testing', 'std', 'twist', ...
        'alpha_beta_std', 'alpha_beta_twist'};

    if any(ismember(algorithm, default_algorithms))
        D = BoundaryMatrix(stream, mode);
    elseif any(ismember(algorithm, {'rho_std', 'rho_twist'}))
        D = ReductionRho(stream, mode);
    elseif any(ismember(algorithm, {'c8_std', 'c8_twist'}))
        D = ReductionC8(stream, mode);
    elseif any(ismember(algorithm, {'alpha_beta_parallel'}))
        D = ReductionParallel(stream, mode);
    else
        error('Algorithm not identified\n');
    end

    if as_dense
        D.as_dense();
    end

end

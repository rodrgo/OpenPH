% ESSENTIAL_STD.M
% Reduce Boundary matrix to find only essential simplices 

function [ess, t] = essential_std(D)
    t0 = tic;
    D.init();

    % Identify essential by Curiosity 7
    ess_candidates = ones(1, D.m);
    for j = 1:D.m
        if D.low(j) > 0
            ess_candidates(D.low(j)) = 0;
        end
    end

    % Remove negative columns that have a 'left'
    % If column has a 'left', it has 'beta' > 0
    D.create_beta();
    neg = D.beta > 0;
    ess_candidates(neg) = 0;

    % Refine ess_candidates by looking at alpha/beta info
    D.create_alpha();
    neg_pairs = (D.alpha == D.beta) & (D.beta > 0);
    pos_pairs = D.beta(neg_pairs);
    ess_candidates(pos_pairs) = 0;

    % Extract essentials
    ess = find(ess_candidates);
    t = toc(t0);
end

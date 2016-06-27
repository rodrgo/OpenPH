% ESSENTIAL_RED.M
% Reduce Boundary matrix to find only essential simplices 

function [ess, t] = essential_red(D)
    t0 = tic;
    D.init();
    ess_candidates = ones(1, D.m);
    for j = 1:D.m
        if D.low(j) > 0
            ess_candidates(D.low(j)) = 0;
        end
    end
    ess_candidates = find(ess_candidates);
    % Extract essentials
    ess = ess_candidates;
    t = toc(t0);
end

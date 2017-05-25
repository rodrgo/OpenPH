% TWIST_RED.M
% 

function [lows, t] = twist_red(D)
    t0 = tic;
    % Initialise persistence vectors
    D.init();
    D.record_iteration();
    % Start algorithm
    complex_dim = D.complex_dimension;
    for d = complex_dim:-1:1
        d_simplices = D.get_d_simplices(d);
        for j = d_simplices
            D.reduce_col_twist(j);
            D.record_iteration();
        end
    end
    % Extract lows
    lows = D.low;
    t = toc(t0);
end

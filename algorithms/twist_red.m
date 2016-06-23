% TWIST_RED.M
% 

function [lows, t] = twist_red(D)
    t0 = tic;
    % Initialise persistence vectors
    D.init();
    % Start algorithm
    complex_dim = D.complex_dimension;
    for d = complex_dim:-1:1
        d_simplices = D.get_d_simplices(d);
        for j = d_simplices
            D.reduce_col(j);
            if D.low(j) > 0
                i = D.low(j);
                D.clear_cols(i);
            end
        end
    end
    % Extract lows
    lows = D.low;
    t = toc(t0);
end

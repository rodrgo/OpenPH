% STD_RED.M
% Standard reduction algorithm

function [lows, t] = std_red(D)
    t0 = tic;
    % Initialise persistence vectors
    D.init();
    % Reduce
    for j = 1:D.m
        D.reduce_col(j);
    end
    % Extract lows
    lows = D.low;
    t = toc(t0);
end

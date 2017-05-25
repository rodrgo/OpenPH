% STD_RED.M
% Standard reduction algorithm

function [lows, t] = std_red(D)
    t0 = tic;
    % Initialise persistence vectors
    D.init();
    D.record_iteration();
    % Reduce
    for j = 1:D.m
        D.reduce_col(j);
	D.record_iteration();
    end
    % Extract lows
    lows = D.low;
    t = toc(t0);
end

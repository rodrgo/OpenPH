% PH_ROW.M
% Row reduction algorithm

function [lows, t] = ph_row(D)
    t0 = tic;
    % Initialise persistence vectors
    D.init();
    D.record_iteration();
    % Reduce
    for i = D.m:-1:1
        indices = find(D.low == i);
        for j = indices(2:end)
            D.left_to_right(indices(1), j);
	        D.record_iteration();
        end
    end
    % Extract lows
    lows = D.low;
    t = toc(t0);
end

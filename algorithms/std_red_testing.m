% STD_RED_TESTING.M
%

function [lows, t] = std_red_testing(R)
    tic;
    m = size(R, 2);
    lows = zeros(1, m);

    % Create lows vector
    for j = 1:m
        if any(R(:, j))
            l = find(R(:, j), 1, 'last');
        else
            l = 0;
        end
        lows(j) = l;
    end

    % Reduce columns
    for j = 1:m
        j0 = find(lows(1:(j-1)) ~= 0 & lows(1:(j-1)) == lows(j));
        while ~isempty(j0)
            R(:, j) = mod(R(:, j) + R(:, j0), 2);
            if any(R(:, j))
                lows(j) = find(R(:, j), 1, 'last');
            else
                lows(j) = 0;
            end
            j0 = find(lows(1:(j-1)) ~= 0 & lows(1:(j-1)) == lows(j));
        end
    end
    t = toc;
end

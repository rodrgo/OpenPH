
init;

folder = '../datasets/pointclouds';
fid = fopen('pointclouds_stream_success.log');
tline = fgetl(fid);

COL_WIDTH = 5;
while ischar(tline)
    if endsWith(tline, '.txt')
        disp(tline);
        fpath = [folder '/' tline];
        display(tline);

        d = 5; % max_dimension 
        n = 5; % num_steps
        p = 3; % max_filtration_value

        [r, c, m] = pointcloud_to_boundary_matrix(fpath, d, n, p);

        % TODO: Check dat2cmo and cmo2dat
        %[lows, t] = std_red_testing(r, c, m);
        [OUT, t] = openph(r, c, m, 'standard', COL_WIDTH, zeros(1,m));
        display(t);
        [OUT, t] = openph(r, c, m, 'pms', COL_WIDTH, zeros(1,m));
        display(t);
    end
end
fclose(fid);

%function [lows, t] = std_red_testing(r, c, m)
%    R = full(sparse(r, c, ones(size(r)), m, m));
%    tic;
%    m = size(R, 2);
%    lows = zeros(1, m);
%
%    % Create lows vector
%    for j = 1:m
%        if any(R(:, j))
%            l = find(R(:, j), 1, 'last');
%        else
%            l = 0;
%        end
%        lows(j) = l;
%    end
%
%    % Reduce columns
%    for j = 1:m
%        j0 = find(lows(1:(j-1)) ~= 0 & lows(1:(j-1)) == lows(j));
%        while ~isempty(j0)
%            R(:, j) = mod(R(:, j) + R(:, j0), 2);
%            if any(R(:, j))
%                lows(j) = find(R(:, j), 1, 'last');
%            else
%                lows(j) = 0;
%            end
%            j0 = find(lows(1:(j-1)) ~= 0 & lows(1:(j-1)) == lows(j));
%        end
%    end
%    t = toc;
%end
%

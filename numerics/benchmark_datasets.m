
init;

folder = '../datasets/pointclouds';
fid = fopen(fullfile('pointclouds_stream_success.list'), 'r');
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


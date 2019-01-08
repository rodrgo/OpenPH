
init;
listing = dir('../datasets/pointclouds');

COL_WIDTH = 5;

ts_phat = [];
ts_pms = [];
ms = [];
datasets = {};

for i = 1:length(listing)
    if endsWith(listing(i).name, '.txt')
        fpath = [listing(i).folder '/' listing(i).name];
        display(listing(i).name);
        d = 5; % max_dimension 
        n = 5; % num_steps
        p = 3; % max_filtration_value
        try
            % Read pointcloud
            [r, c, m] = pointcloud_to_boundary_matrix(fpath, d, n, p);
            fprintf(fid, [listing(i).name '\n']);
            display(sprintf('Conversion successful. m=%d', m));

            % PMS
            tic_pms = tic;
            [OUT, t] = openph(r, c, m, 'pms', COL_WIDTH, zeros(1,m));
            toc_pms = toc(tic_pms)
            display(sprintf('PMS success in %g', toc_pms));

            % PHAT
            tic_phat = tic;
            [lows_phat, t] = phat(r, c, m, PHAT_DIR, 'twist', false);
            toc_phat = toc(tic_phat)
            display(sprintf('PHAT success in %g', toc_phat));

            assert(sum(OUT.low ~= lows_phat) == 0);

            ts_pms(end+1) = toc_pms;
            ts_phat(end+1) = toc_phat;
            ms(end+1) = m;
            datasets{end+1} = listing(i).name;

        catch e
            display('Error');
            %e.message
            if(isa(e,'java.lang.OutOfMemoryError'))
                ex = e.ExceptionObject;
                assert(isjava(ex));
                ex.printStackTrace;
            end
        end
    end
end

% =============
% Create table
% =============

results_tag = strcat('table_datasets.tex');
table_path = fullfile(FIGURE_DIR, results_tag);

fid = fopen(table_path, 'w');

fprintf(fid,'\\begin{tabular}{l l || c c }\n');
fprintf(fid,'\\toprule\n');

fprintf(fid,'Dataset & $m$ & PMS & PHAT \\\\ \n');
fprintf(fid,'\\midrule\n');

for j = 1:length(datasets)
    fprintf(fid,'%s & %d & %4.3f & 4.3f \\\\ \n', datasets{j}, ms(j), ts_pms(j), ts_phat(j));
end 
fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}');

% fclose
fclose(fid);


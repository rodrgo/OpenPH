
init;
listing = dir('./datasets');

COL_WIDTH = 8;

ts_phat = [];
ts_pms = [];
ms = [];
datasets = {};

for i = 1:length(listing)
    if endsWith(listing(i).name, '.dat')
        fpath = [listing(i).folder '/' listing(i).name];
        display(listing(i).name);

        % Convert to CMO 
        % Read pointcloud
        [r, c, m] = dat2cmo(fpath);

        % PMS
        [OUT, toc_pms] = openph(r, c, m, 'pms', COL_WIDTH, zeros(1,m));
        display(sprintf('PMS success in %g', toc_pms));

        % PHAT
        [lows_phat, toc_phat] = phat(r, c, m, PHAT_DIR, 'twist', false);
        display(sprintf('PHAT success in %g', toc_phat));

        assert(sum(OUT.low ~= lows_phat) == 0);

        ts_pms(end+1) = toc_pms;
        ts_phat(end+1) = toc_phat;
        ms(end+1) = m;
        datasets{end+1} = listing(i).name;

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


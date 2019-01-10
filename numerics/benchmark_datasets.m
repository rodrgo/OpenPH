
init;
listing = dir('./datasets');

ts_phat = [];
ts_pms = [];
ts_pms_inner = [];
col_widths = [];
ms = [];
datasets = {};

% Warm up
for i = 1:length(listing)
    if endsWith(listing(i).name, '.dat')
        fpath = [listing(i).folder '/' listing(i).name];
        COL_WIDTH = 65;
        [r, c, m] = dat2cmo(fpath);
        [OUT, toc_pms] = openph(r(r>0), c(r>0), m, 'pms', COL_WIDTH, zeros(1,m));
        break
    end
end

% Continue

for i = 1:length(listing)
    if endsWith(listing(i).name, '.dat')
        fpath = [listing(i).folder '/' listing(i).name];
        display(listing(i).name);
        

        switch fpath
            case [listing(i).folder '/primal_explicit.complex.dat']
                COL_WIDTH = 10;
            case [listing(i).folder '/noise_3_16.complex.dat']
                COL_WIDTH = 30;
            case [listing(i).folder '/smooth_16.complex.dat']
                COL_WIDTH = 30;
            case [listing(i).folder '/ramp_3_16.complex.dat']
                COL_WIDTH = 10;
            case [listing(i).folder '/ramp_4_8.complex.dat']
                COL_WIDTH = 10;
            case [listing(i).folder '/noise_4_8.complex.dat']
                COL_WIDTH = 65;
            otherwise
                disp('other value')
        end

        % Convert to CMO 
        % Read pointcloud
        [r, c, m] = dat2cmo(fpath);

        % PMS
        [OUT, toc_pms] = openph(r(r>0), c(r>0), m, 'pms', COL_WIDTH, zeros(1,m));
        lows = OUT.low;
        display(sprintf('PMS: %g secs ', toc_pms));
        display(sprintf('PMS: %g ms (inner)', OUT.time_inner));

        % PHAT
        [lows_phat, toc_phat] = phat(r, c, m, PHAT_DIR, 'twist', false);
        display(sprintf('PHAT: %g secs', toc_phat));

        display(sprintf('Offsets = %d', sum(lows ~= lows_phat)));

        col_widths(end+1) = COL_WIDTH;
        ts_pms(end+1) = toc_pms;
        ts_pms_inner(end+1) = OUT.time_inner;
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

fprintf(fid,'\\begin{tabular}{l l || c | c c c }\n');
fprintf(fid,'\\toprule\n');

fprintf(fid,'Dataset & $m$ & PMS (width) & PMS (total) & PMS (inner) & PHAT \\\\ \n');
fprintf(fid,'\\midrule\n');

for j = 1:length(datasets)
    fprintf(fid,'%s & %d & %d & %d & %d & %d \\\\ \n', ...
        strrep(strrep(datasets{j}, '_', '-'), '.complex.dat', ''), ms(j), col_widths(j), ...
        round(ts_pms(j)*1000, 0), ...
        round(ts_pms_inner(j), 0), ...
        round(ts_phat(j)*1000, 0));
end 
fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}');

% fclose
fclose(fid);


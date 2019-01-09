init;

% Complex parameters
MAX_DIM = 20;
num_points = 5:19;

% Complex parameters
cparams               = [];
cparams.max_dim       = MAX_DIM;
cparams.num_divs      = 10;
cparams.max_filtr_val = 5;
cparams.num_points    = 5;

shape = 'random_gaussian';

phat_algos = {'twist', 'spectral_sequence', 'row'};

num_algos   = 1+2*length(phat_algos);
num_shapes  = 1;
num_samples = 3;
COL_WIDTH   = 8;

% Create shapes
ts  = zeros(num_algos, length(num_points));
nps = zeros(1, length(num_points));
ms  = zeros(1, length(num_points));
colors = create_color_palette(num_algos);

for j = 1:length(num_points)
    np = num_points(j);
    cparams.num_points = np;
    [stream, complex_info] = complex_factory(shape, cparams);
    [r, c, m] = stream2cmo(stream);

    nps(j) = np;
    ms(j) = m;

    %[lows, t] =  std_red_testing(full(sparse(r, c, ones(size(r)), m, m)));
    [OUT, time_pms] = openph(r, c, m, 'pms', COL_WIDTH, zeros(1,m));
    ts(1, j) = time_pms;
    lows = OUT.low;

    for k = 1:length(phat_algos)
        % Without dualization
        [lows_phat, t] = phat(r, c, m, PHAT_DIR, phat_algos{k}, false);
        ts(1+k, j) = t;
        assert(sum(lows ~= lows_phat) == 0);

        % With dualization
        [lows_phat, t] = phat(r, c, m, PHAT_DIR, phat_algos{k}, true);
        ts(1+length(phat_algos)+k, j) = t;
        assert(sum(lows ~= lows_phat) == 0);
    end

end 

% =============
% Create table
% =============

table_tag = strcat('table_phat.tex');
table_path = fullfile(FIGURE_DIR, table_tag);

fid = fopen(table_path, 'w');

fprintf(fid,'\\begin{tabular}{l l ||');
fprintf(fid,'S'); % pms
for k = 1:length(phat_algos)
    fprintf(fid,'SS');
end
fprintf(fid, '}\n');
fprintf(fid,'\\toprule\n');

fprintf(fid,'\\multirow{2}{*}{$N$} & \n');
fprintf(fid,'\\multirow{2}{*}{$m$} & \n');
fprintf(fid,'\\multicolumn{1}{c}{pms} & \n');
for k = 1:length(phat_algos)
    fprintf(fid,'\\multicolumn{2}{c}{%s}', strrep(phat_algos{k}, '_', '-'));
    if k < length(phat_algos)
        fprintf(fid,'& \n');
    else
        fprintf(fid,'\\\\\n');
    end
end

fprintf(fid,'& & '); % n & m & pms
for k = 1:length(phat_algos)
    fprintf(fid,' & {P} & {D} ');
end
fprintf(fid,'\\\\\n');
fprintf(fid,'\\midrule\n');

for j = 1:length(nps)
    fprintf(fid,'%d & %d', nps(j), ms(j));
    fprintf(fid,'& %d', round(ts(1,j)*1000, 0));
    for k = 1:length(phat_algos)
        fprintf(fid, '& %d & %d', round(ts(1+k,j)*1000, 0), round(ts(1+length(phat_algos)+k,j)*1000, 0));
    end
    fprintf(fid, '\\\\\n');
end 
fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}');

% fclose
fclose(fid);


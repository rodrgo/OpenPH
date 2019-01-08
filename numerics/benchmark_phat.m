
init;
plot_init;

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

phat_algos = {'twist', 'chunk', 'spectral_sequence', 'row'};

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

    test_tic = tic;
    %[OUT, t] = openph(r, c, m, 'pms', COL_WIDTH, zeros(1,m));
    [lows, t] =  std_red_testing(full(sparse(r, c, ones(size(r)), m, m)));
    ts(1, j) = toc(test_tic)
    %lows = OUT.low;

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
% Create figure
% =============

table_tag = strcat('table_phat.tex');
table_path = fullfile(FIGURE_DIR, table_tag);

fileId = fopen(table_path, 'w');

fprintf(fileId,'\\begin{tabular}{l l ||');
for k = 1:length(phat_algos)
    fprintf(fileId,'SS');
end
fprintf(fileId, '}\n');
fprintf(fileId,'\\toprule\n');

fprintf(fileId,'\\multirow{2}{*}{$N$} & \n');
fprintf(fileId,'\\multirow{2}{*}{$m$} & \n');
for k = 1:length(phat_algos)
    fprintf(fileId,'\\multicolumn{2}{c}{%s}', strrep(phat_algos{k}, '_', '-'));
    if k < length(phat_algos)
        fprintf(fileId,'& \n');
    else
        fprintf(fileId,'\\\\\n');
    end
end

fprintf(fileId,'& ');
for k = 1:length(phat_algos)
    fprintf(fileId,' & {P} & {D} ');
end
fprintf(fileId,'\\\\\n');
fprintf(fileId,'\\midrule\n');

for j = 1:length(nps)
    fprintf(fileId,'%d & %d', nps(j), ms(j));
    for k = 1:length(phat_algos)
        fprintf(fileId, '& %d & %d', round(ts(1+k,j)*1000, 0), round(ts(1+length(phat_algos)+k,j)*1000, 0));
    end
    fprintf(fileId, '\\\\\n');
end 
fprintf(fileId,'\\bottomrule\n');
fprintf(fileId,'\\end{tabular}');

% fclose
fclose(fileId);



init;
plot_init;

%folder = '../datasets/pointclouds';
%fid = fopen(fullfile('pointclouds_stream_success.list'), 'r');
%tline = fgetl(fid);

%[r, c, m] = dat2cmo(fullfile('single_triangle.dat'));
%low = phat(r, c, m, PHAT_DIR);
%low

%[lows, t] = std_red_testing(R);
%lows

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

num_algos   = 2;
num_shapes  = 1;
num_samples = 3;
COL_WIDTH   = 8;

% Create shapes
ts  = zeros(num_algos, length(num_points));
nps = zeros(1, length(num_points));
ms  = zeros(1, length(num_points));
colors = create_color_palette(num_algos);

phat_algos = {'standard', 'twist', 'chunk', 'chunk_sequential', 'spectral_sequence', 'row'};

for j = 1:length(num_points)
    np = num_points(j);
    cparams.num_points = np;
    [stream, complex_info] = complex_factory(shape, cparams);
    [r, c, m] = stream2cmo(stream);

    nps(j) = np;
    ms(j) = m;

    test_tic = tic;
    [OUT, t] = openph(r, c, m, 'pms', COL_WIDTH, zeros(1,m));
    ts(1, j) = toc(test_tic)
    lows = OUT.low;

    for k = 1:length(phat_algos)
        % Without dualization
        [lows_phat, t] = phat(r, c, m, PHAT_DIR, phat_algos{k}, false);
        ts(1+k, j) = t;

        % With dualization
        [lows_phat, t] = phat(r, c, m, PHAT_DIR, phat_algos{k}, true);
        ts(1+length(phat_algos)+k, j) = t;
    end

    assert(sum(lows ~= lows_phat) == 0);
end 


% =============
% Create figure
% =============

table_tag = strcat('phat_benchmark.tex');
table_path = fullfile(FIGURE_DIR, table_tag);

% fopen
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
    fprintf(fileId,'\\multicolumn{2}{c}{A (\%)}');
    if k < length(phat_algos)
        fprintf(fileId,'& \n');
    else
        fprintf(fileId,'\\\\\n');
    end
end

fprintf(fileId,'& \n');
for k = 1:length(phat_algos)
    fprintf(fileId,' & p & d ');
end
fprintf(fileId,'\\\\\n');
fprintf(fileId,'\\midrule\n');

for j = 1:length(nps)
    fprintf(fileId,'%d & %d', nps(j), ms(j));
    for k = 1:length(phat_algos)
        fprintf(fileId, '& %3.6f & %3.6f', ts(1+k,j), ts(1+length(phat_algos)+k,j));
    end
    fprintf(fileId, '\\\\\n');
end 
fprintf(fileId,'\\bottomrule\n');
fprintf(fileId,'\\{tabular}');

% fclose
fclose(fileId);


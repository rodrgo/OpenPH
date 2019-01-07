
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

for j = 1:length(num_points)
    np = num_points(j);
    cparams.num_points = np;
    [stream, complex_info] = complex_factory(shape, cparams);
    [r, c, m] = stream2cmo(stream);

    nps(j) = np;
    ms(j) = m;

    test_tic = tic;
    [OUT, t] = openph(r, c, m, 'pms', COL_WIDTH, zeros(1,m));
    %[lows, t] = std_red_testing(full(sparse(r, c, ones(size(r)), m, m)));
    ts(1, j) = toc(test_tic)
    lows = OUT.low;

    [lows_phat, t] = phat(r, c, m, PHAT_DIR);
    ts(2, j) = t

    %[OUT, t] = openph(r, c, m, algos{l}, COL_WIDTH, low_true);
    % t = sum(OUT.time_track(1:(OUT.num_iters)));

    assert(sum(lows ~= lows_phat) == 0);
end 


% =============
% Create figure
% =============

table_tag = strcat('phat_benchmark.tex');
table_path = fullfile(FIGURE_DIR, table_tag);

fileId = fopen(table_path, 'w');

fprintf(fileId,'\\begin{tabular}{l l || c | c}\n');
fprintf(fileId,'\\toprule\n');
fprintf(fileId,'{$N$} & {$m$} & {pms} & {phat}\n');
fprintf(fileId,'\\midrule\n');
for j = 1:length(nps)
    fprintf(fileId,'%d & %d & %3.6f & %3.6f\\\\\n', ...
        nps(j), ms(j), ts(1,j), ts(2,j));
end 
fprintf(fileId,'\\bottomrule\n');
fprintf(fileId,'\\{tabular}');

fclose(fileId);


%
%figure(1);
%set(gcf, 'color', [1 1 1]);
%set(gca, 'Fontname', 'setTimes', 'Fontsize', 18);
%handles = [];
%labels  = {};
%
%algos = {'pms', 'phat'};
%for l = 1:length(algos)
%    labels{end + 1}  = strrep(algos{l}, '_', '\_');
%    hold on;
%    handles(end + 1) = loglog(ms, ts(l,:), '-*', 'Color', colors{l});
%
%    % Arrow
%    ll=ceil(length(ms)*l/length(algos));
%    txt=['\leftarrow ' strrep(algos{l}, '_', '\_')];
%    h=text(ms(ll), ts(l,ll), txt, 'HorizontalAlignment', 'left');
%    set(h, 'Color', colors{l});
%end
%
%% Cosmetics
%xlabel('m', 'FontSize', fs.axis);
%ylabel('time (ms)', 'FontSize', fs.axis);
%
%%legend(handles, labels, 'FontSize', fs.legend);
%% Tick size
%xt = get(gca, 'XTick');
%set(gca, 'FontSize', fs.ticks);
%
%xt = get(gca, 'YTick');
%set(gca, 'FontSize', fs.ticks);
%
%title('PMS VS PHAT', 'FontSize', fs.title);
%filepath = fullfile(FIGURE_DIR, 'phat.eps');
%print('-depsc', filepath);
%eps_to_pdf(filepath);
%

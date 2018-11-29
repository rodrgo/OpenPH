% test_speed_standard_reduction.m
init;
plot_init;

% Complex parameters
MAX_DIM = 20;
num_points = 5:16

cparams               = [];
cparams.max_dim       = MAX_DIM;
cparams.num_divs      = 10;
cparams.max_filtr_val = 5;
cparams.num_points    = 5;

% Explore complexity of Vietoris-Rips complexes
shape = 'random_gaussian';

%algos = {'standard', 'twist', 'ph_row', 'standard_parallel', 'twist_parallel', 'ph_row_parallel', 'pms'};
algos = {'standard', 'twist', 'ph_row', 'pms'};

% Matrix dense?
as_dense  = true;
time_init = tic;

% Create figure
figure(1);
set(gcf, 'color', [1 1 1]);
set(gca, 'Fontname', 'setTimes', 'Fontsize', 18);
handles = [];
labels  = {};

COL_WIDTH = 8;

% Create shapes
ts  = zeros(length(algos), length(num_points));
nps = zeros(1, length(num_points));
ms  = zeros(1, length(num_points));
colors = create_color_palette(length(algos));

for j = 1:length(num_points)
    np = num_points(j);
    cparams.num_points = np;
    [stream, complex_info] = complex_factory(shape, cparams);
    [r, c, m] = stream2cmo(stream);
    low_true = get_true_low(r, c, m, COL_WIDTH);

    if false
        D = BoundaryMatrix(stream);
        [r_verify, c_verify] = find(D.matrix);
        assert(all(r == r_verify));
        assert(all(c == c_verify));
    end

    nps(j) = np;
    ms(j) = m;

    for l = 1:length(algos)

        [OUT, ~] = openph(r, c, m, algos{l}, COL_WIDTH, low_true);
        t = sum(OUT.time_track(1:(OUT.num_iters)));
        ts(l, j) = t;
        
        fprintf('%s\t%s: (np, m, t) = (%d, %d, %g)\n',...
            shape, algos{l}, cparams.num_points, m, t);

    end
    fprintf('\n');
end 

A = [ts; ms];
A = A';
A = sortrows(A, size(A, 2));
A = A';
ts = A(1:(end-1),:);
ms = A(end,:);

for l = 1:length(algos)
    labels{end + 1}  = strrep(algos{l}, '_', '\_');
    hold on;
    handles(end + 1) = loglog(ms, ts(l,:), '-*', 'Color', colors{l});

    % Arrow
    ll=ceil(length(ms)*l/length(algos));
    txt=['\leftarrow ' strrep(algos{l}, '_', '\_')];
    h=text(ms(ll), ts(l,ll), txt, 'HorizontalAlignment', 'left');
    set(h, 'Color', colors{l});
end

% Add scaling references 
colour = 'black';

if true
    slope = 2;
    ordinate = 0;
    x = [ms(end-2), ms(end)];
    y = exp(ordinate)*x.^(slope);
    y = y/(y(end))*max(max(ts))*1.1;
    labels{end + 1}  = strrep('Quadratic', '_', '\_');
    hold on;
    handles(end + 1) = loglog(x, y, '-', 'Color', colour);
    h=text(x(end), y(end), 'O(m^2)', 'HorizontalAlignment', 'right');
    set(h, 'Color', colour);
end

if true
    slope = 1;
    ordinate = 0;
    x = [ms(end-3), ms(end-2)];
    y = exp(ordinate)*x.^(slope);
    y = y/(y(end))*max(max(ts))*0.5;
    labels{end + 1} = strrep('Linear', '_', '\_');
    hold on;
    handles(end + 1) = loglog(x, y, '-', 'Color', colour);
    h=text(x(end), y(end), 'O(m)', 'HorizontalAlignment', 'right');
    set(h, 'Color', colour);
end

if true
    slope = 1;
    ordinate = 0;
    x = [ms(end-1), ms(end)];
    y = x;
    labels{end + 1} = strrep('Logarithmic', '_', '\_');
    hold on;
    handles(end + 1) = loglog(x, y, '-', 'Color', colour);
    h=text(x(end), y(end), 'O(log(m))', 'HorizontalAlignment', 'right');
    set(h, 'Color', colour);
end

% Cosmetics

xlabel('m', 'FontSize', fs.axis);
ylabel('time (ms)', 'FontSize', fs.axis);

%legend(handles, labels, 'FontSize', fs.legend);

title('Scaling of algorithms.', 'FontSize', fs.title);
filepath = fullfile(FIGURE_DIR, 'scalings.eps');
print('-depsc', filepath);
eps_to_pdf(filepath);



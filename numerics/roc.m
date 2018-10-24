% test_speed_standard_reduction.m
init;

% Complex parameters
cparams               = [];
cparams.max_dim       = 5;
cparams.num_divs      = 10;
cparams.max_filtr_val = 5;
cparams.num_points    = 5;

num_points = 10:21;

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

% Create shapes

ts = zeros(length(algos), length(num_points));

nps = zeros(1, length(num_points));
ms = zeros(1, length(num_points));
colors = create_color_palette(length(algos));
for j = 1:length(num_points)
    np = num_points(j);
    cparams.num_points = np;
    [stream, complex_info] = complex_factory(shape, cparams);
    D = BoundaryMatrix(stream);

    %low_true = reduce_stream(stream, 'testing', as_dense);
    OUT = cuda_wrapper(D, 'pms', zeros(1,D.m), 7);
    low_true = OUT.low;

    nps(j) = np;
    ms(j) = D.m;

    for l = 1:length(algos)
        fprintf('\t\t\t%s... ', algos{l});

        [OUT, ~] = cuda_wrapper(D, algos{l}, low_true, 7);
        assert(all(OUT.low == low_true), 'Output incorrect!');

        t = sum(OUT.time_track(1:OUT.num_iters));
        if (t ~= sum(OUT.time_track))
            display(t);
            display(sum(OUT.time_track));
        end
        ts(l, j) = t/1000;
        
        fprintf('\n\t%s: (np, m, t) = (%d, %d, %g)\n',...
            shape, cparams.num_points, D.m, t);

    end
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

x = ms(end-2):100:ms(end);
y = 3*(x.^2)/max(x.^2);
labels{end + 1}  = strrep('Quadratic', '_', '\_');
hold on;
handles(end + 1) = loglog(x, y, '-', 'Color', colour);
ll=ceil(length(x)*3/4);
h=text(x(ll), y(ll), 'O(m^2)', 'HorizontalAlignment', 'right');
set(h, 'Color', colour);

x = ms(end-2):100:ms(end);
y = 1*x/max(x);
labels{end + 1} = strrep('Linear', '_', '\_');
hold on;
handles(end + 1) = loglog(x, y, '-', 'Color', colour);
ll=ceil(length(x)*3/4);
h=text(x(ll), y(ll), 'O(m)', 'HorizontalAlignment', 'right');
set(h, 'Color', 'black');

% Cosmetics

xlabel('m', 'FontSize', fs.axis);
ylabel('time (sec)', 'FontSize', fs.axis);

%legend(handles, labels, 'FontSize', fs.legend);

title('Scaling of algorithms', 'FontSize', fs.title);
filepath = fullfile(FIGURE_DIR, 'scalings.eps');
print('-depsc', filepath);
eps_to_pdf(filepath);


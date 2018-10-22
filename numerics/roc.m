% test_speed_standard_reduction.m
init;

FIGURE_DIR      = '../figures/';
EXPERIMENT_TAG  = 'benchmark_pms';

% Font size
fs          = [];
fs.title    = 20;
fs.legend   = 17;
fs.axis     = 20;
fs.ticks    = 20;

% Complex parameters
cparams               = [];
cparams.max_dim       = 5;
cparams.num_divs      = 10;
cparams.max_filtr_val = 5;
cparams.num_points    = 15;

num_points = 10:19;

% Explore complexity of Vietoris-Rips complexes
shape = 'random_gaussian';

%algos = {'standard', 'twist', 'ph_row', 'standard_parallel', 'twist_parallel', 'ph_row_parallel', 'pms'};
algos = {'standard', 'twist', 'ph_row', 'pms'};

% Matrix dense?
as_dense  = true;
time_init = tic;

algos_map = containers.Map;

% Create figure
figure(1);
set(gcf, 'color', [1 1 1]);
set(gca, 'Fontname', 'setTimes', 'Fontsize', 18);
handles = [];
labels  = {};

colors = create_color_palette(length(algos));
for l = 1:length(algos)
    fprintf('\t\t\t%s... ', algos{l});

    ms = [];
    ts = [];
    nps = [];
    for np = num_points
        cparams.num_points = np;

        [stream, complex_info] = complex_factory(shape, cparams);
        low_true = reduce_stream(stream, 'testing', as_dense);
        D = BoundaryMatrix(stream);
        [OUT, ~] = cuda_wrapper(D, algos{l}, low_true, 7);

        t = sum(OUT.time_track(1:OUT.num_iters));
        if (t ~= sum(OUT.time_track))
            display(t);
            display(sum(OUT.time_track));
        end
        
        nps(end+1) = np;
        ms(end+1) = D.m;
        ts(end+1) = t/1000;
        assert(all(OUT.low == low_true), 'Output incorrect!');

        fprintf('\n\t%s: (np, m, t) = (%d, %d, %g)\n',...
            shape, cparams.num_points, D.m, t);

    end

    labels{end + 1}  = strrep(algos{l}, '_', '\_');
    hold on;
    handles(end + 1) = loglog(ms, ts, '--', 'Color', colors{l});

    % Arrow
    ll=ceil(length(ms)*l/length(algos));
    txt=['\leftarrow ' strrep(algos{l}, '_', '\_')];
    h=text(ms(ll), ts(ll), txt, 'HorizontalAlignment', 'left');
    set(h, 'Color', colors{l});

end % end algorithms

% Add scaling references 
x = ms(1:end/2);
colour = 'black';

y = 1e-7*x.^(3/2);
labels{end + 1}  = strrep('Superlinear', '_', '\_');
hold on;
handles(end + 1) = loglog(x, y, '--', 'Color', colour);
ll=ceil(length(x)*3/4);
h=text(x(ll), y(ll), 'O(m^(3/2))', 'HorizontalAlignment', 'right');
set(h, 'Color', colour);

y = 1e-7*x.^(1);
labels{end + 1} = strrep('Linear', '_', '\_');
hold on;
handles(end + 1) = loglog(x, y, '--', 'Color', colour);
ll=ceil(length(x)*3/4);
h=text(x(ll), y(ll), 'O(m)', 'HorizontalAlignment', 'right');
set(h, 'Color', 'black');

% Cosmetics

xlabel('m', 'FontSize', fs.axis);
ylabel('time (sec)', 'FontSize', fs.axis);

%legend(handles, labels, 'FontSize', fs.legend);

title('Scaling', 'FontSize', fs.title);
filepath = fullfile(FIGURE_DIR, 'scalings.eps');
print('-depsc', filepath);
eps_to_pdf(filepath);


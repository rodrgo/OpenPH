% SPY_TDA.M
% Visualise sparsity pattern for boundary matrices
% and return handle

function plot_handle = spy_tda(S, LineSpec, markersize)

	% Default variables

	if nargin < 3 || markersize == 0
		markersize = get(gcf, 'defaultlinemarkersize');
	end

	if nargin < 2
		LineSpec = '.k';
	end

	[~, color, marker] = colstyle(LineSpec);
	if isempty(marker)
		marker = '.';
	end
	if isempty(color)
		color = [0,0,0];
	end

	linestyle = 'none';

	% We must have three inputs

	[m, n] = size(S);
	[rows, cols] = find(S);

	if isempty(rows)
		rows = NaN;
		cols = NaN;
	end
	if isempty(S)
		marker = 'none';
	end
	plot_handle = plot(cols, rows, ...
		'marker', marker, ...
		'markersize', markersize, ...
		'linestyle', linestyle, ...
		'color', color);
	xlabel(['nnz = ' int2str(nnz(S))]);
	set(gca,'xlim', [0 n+1], 'ylim', [0 m+1], ...
		'ydir','reverse', ...
		'GridLineStyle','none', ...
		'plotboxaspectratio', [n+1 m+1 1]);

end

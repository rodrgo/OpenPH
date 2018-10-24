function palette_cell = create_color_palette(alg_list)

	numColors = 9;
	if length(alg_list) > 9
		numColors = length(alg_list);
	end

	palette = colorscale(numColors, 'hue', [1/100 1], 'saturation' , 1, 'value', 0.7);
	palette(1, :) = [0 0 0];
	palette(4, :) = palette(4, :) + [0 -0.2 0];
	idx = 1:numColors;
	idx(1:9) = [1 6 9 8 5 4 2 7 3];

	palette_cell = cell(length(idx), 1);
	for i = 1:length(idx)
		palette_cell{i} = palette(idx(i),:);
	end

end

function cols = colorscale(n, varargin)
    p = inputParser;
    p.addRequired('n', @isnumeric);
    p.addOptional('hue', [0.1 0.9], @(x) length(x) == 2 & min(x) >=0 & max(x) <= 1);
    p.addOptional('saturation', 0.5, @(x) length(x) == 1);
    p.addOptional('value', 0.8, @(x) length(x) == 1);
    p.parse(n, varargin{:});
    cols = hsv2rgb([transpose(linspace(p.Results.hue(1), p.Results.hue(2), p.Results.n)), ...
        repmat(p.Results.saturation, p.Results.n, 1), ...
        repmat(p.Results.value, n,1) ]);
end

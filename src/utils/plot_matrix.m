function h = plot_matrix(matrix, title_str, file_path, masks)

h = spy_tda(matrix);
set(gcf, 'color', [1 1 1]);
set(gca, 'xtick', [], 'ytick', [], 'XTickLabel', '', 'YTickLabel', '');
xlabel('');

% Plot masks (if any are given)
if nargin > 3
    handles = [];
    names = {};
    for i = 1:length(masks)
        mat = masks{i}{1};
        name = masks{i}{2};
        LineSpec = masks{i}{3};
        hold on;
        if size(mat, 1) > 150
            MS = 0.01;
        else
            MS = 0;
        end
        handles(end + 1) = spy_tda(mat, LineSpec, MS);
        names{end + 1} = name;
        xlabel('');
    end
    legend(handles, names, 'Location', 'SouthWest');
    hold off;
end
title(title_str);

% Save file
print('-depsc', file_path);
eps_to_pdf(file_path);

end

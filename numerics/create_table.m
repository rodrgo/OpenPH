
function create_table(shapes, algos, table_type, levels_x, tensor_x)
    init;
    % table_type:
    %   + 'unreduced'
    %   + 'essential'
    %   + 'l1-error'
    %   + 'linf-error'

    EXPERIMENT_TAG = 'table'
    table_tag = strcat(EXPERIMENT_TAG, '_', table_type, '.tex');
    tablepath = fullfile(FIGURE_DIR, figure_tag);

    fileId = fopen(table_path, 'w');
    table_header(fileId);

    num_shapes = length(shapes);
    num_algos = length(algos);

    for j = 1:num_shapes
        for l = 1:num_algos
            fprintf(fileId,'c');
        end
        if j < num_shapes
            fprintf(fileId,'|');
        end
    end

    fprintf(fileId,'}\n');
    fprintf(fileId,'\\toprule\n');

    table_multirow(fileId, table_type);

    for j = 1:num_shapes
        fprintf(fileId,'\\multicolumn{%d}{c}{%s}', num_algos, shapes{j});
        if j < num_shapes
            fprintf(fileId,'&\n');
        else
            fprintf(fileId,'\\\\\n');
        end
    end
    for j = 1:num_shapes
        for l = 1:num_algos
            fprintf(fileId,'& {%s} ', algos{l});
        end
    end
    fprintf(fileId,'\\\\\n');
    fprintf(fileId,'\\midrule\n');

    % row 'it'
    for pp = 1:length(levels_x)
        if strcmp(table_type, 'unreduced') || strcmp(table_type, 'essential')
            fprintf(fileId,'%1.2f', levels_x(pp));
        else
            fprintf(fileId,'%g', levels_x(pp));
        end
        for j = 1:num_shapes
            for k = 1:num_algos
                fprintf(fileId,' & %d', tensor_x(j, 1, k, pp));
            end
        end
        fprintf(fileId,'\\\\\n');
    end

    fprintf(fileId,'\\bottomrule\n');
    fprintf(fileId,'\\end{tabular}\n');

    table_caption(fileId, table_type);

    fprintf(fileId,'\\end{table*}\n');
    fprintf(fileId,'\\end{small}\n');

end

function table_header(fileId)
    fprintf(fileId,'\n\n\n\n');
    fprintf(fileId,'\\begin{small}\n');
    fprintf(fileId,'\\begin{table*}\n');
    fprintf(fileId,'\\centering\n');
    fprintf(fileId,'\\begin{tabular}{l');
    fprintf(fileId,'||');
end

function table_caption(fileId, table_type)
    if strcmp(table_type, 'unreduced')
        fprintf(fileId,'\\caption{Iterations to unreduced percentage}\n');
        fprintf(fileId,'\\label{tab:iterations_unreduced}\n');
    elseif strcmp(table_type, 'essential')
        fprintf(fileId,'\\caption{Iterations to essential-estimation precision}\n');
        fprintf(fileId,'\\label{tab:iterations_essential}\n');
    elseif strcmp(table_type, 'lone')
        fprintf(fileId,'\\caption{Iterations to relative $\\ell_1$-error}\n');
        fprintf(fileId,'\\label{tab:iterations_l1_error}\n');
    elseif strcmp(table_type, 'linf')
        fprintf(fileId,'\\caption{Iterations to relative $\\ell_\\infty$-error}\n');
        fprintf(fileId,'\\label{tab:iterations_linf_error}\n');
    else
        error('table_type not recognised');
    end
end

function table_multirow(fileId, table_type)
    if strcmp(table_type, 'unreduced')
        fprintf(fileId,'\\multirow{2}{*}{Proportion} &\n');
    elseif strcmp(table_type, 'essential')
        fprintf(fileId,'\\multirow{2}{*}{Precision} &\n');
    elseif strcmp(table_type, 'lone')
        fprintf(fileId,'\\multirow{2}{*}{$\\frac{\\|\\low - \\lowstar\\|_1}{\\|\\lowstar\\|_1}$} &\n');
    elseif strcmp(table_type, 'linf')
        fprintf(fileId,'\\multirow{2}{*}{$\\frac{\\|\\low - \\lowstar\\|_\\infty}{\\|\\lowstar\\|_\\infty}$} &\n');
    else
        error('table_type not recognised');
    end
end

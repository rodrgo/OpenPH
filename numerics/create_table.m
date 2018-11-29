
function create_table(shapes, algos, table_type, levels_x, tensor_x)
    plot_init;
    % table_type:
    %   + 'unreduced'
    %   + 'essential'
    %   + 'l1-error'
    %   + 'linf-error'

    EXPERIMENT_TAG = 'table';
    table_tag = strcat(EXPERIMENT_TAG, '_', table_type, '.tex');
    table_path = fullfile(FIGURE_DIR, table_tag);

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
        fprintf(fileId,'\\multicolumn{%d}{c}{%s}', num_algos, rename_shape(shapes{j}));
        if j < num_shapes
            fprintf(fileId,'&\n');
        else
            fprintf(fileId,'\\\\\n');
        end
    end
    for j = 1:num_shapes
        for l = 1:num_algos
            fprintf(fileId,'& {%s} ', rename_alg(algos{l}));
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

end

function name = rename_shape(shape)
    name = '';
    if strcmp(shape, 'random_gaussian') == 1
       name = 'Gaussian';
    elseif strcmp(shape, 'random_figure_8') == 1
       name = 'Figure-8';
    elseif strcmp(shape, 'random_trefoil_knot') == 1
       name = 'Trefoil-Knot';
    elseif strcmp(shape, 'random_torus') == 1
       name = 'Random-Torus';
    elseif strcmp(shape, 'random_sphere_product') == 1
       name = 'Sphere-Product';
    else
       error('In rename_shape, shape not recognised');
    end
end

function name = rename_alg(alg)
    name = '';
    if strcmp(alg, 'standard') == 1
       name = 'std';
    elseif strcmp(alg, 'twist') == 1
       name = 'twist';
    elseif strcmp(alg, 'ph_row') == 1
       name = 'phRow';
    elseif strcmp(alg, 'standard_parallel') == 1
       name = 'std-parallel';
    elseif strcmp(alg, 'twist_parallel') == 1
       name = 'twist-parallel';
    elseif strcmp(alg, 'ph_row_parallel') == 1
       name = 'phRow-parallel';
    elseif strcmp(alg, 'pms') == 1
       name = 'pms';
    else
       error('In rename_alg, alg not recognised');
    end
end

function table_header(fileId)
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

function create_table_time_ratio(table_path, shapes, algos)

    fileId = fopen(table_path, 'w');

    num_shapes = length(shapes);
    num_algos = length(algos);

    fprintf(fileId,'\n\n\n\n');
    fprintf(fileId,'\\begin{small}\n');
    fprintf(fileId,'\\begin{table*}\n');
    fprintf(fileId,'\\centering\n');
    fprintf(fileId,'\\begin{tabular}{l');
    fprintf(fileId,'||');
    for j = 1:num_shapes
        for l = 1:num_algos
            if strcmp(algos{l}, 'pms') == 0
                fprintf(fileId,'c');
            end
        end
        if j < num_shapes
            fprintf(fileId,'|');
        end
    end
    fprintf(fileId,'}\n');
    fprintf(fileId,'\\toprule\n');
    fprintf(fileId,'\\multirow{2}{*}{Sample} &\n');
    for j = 1:num_shapes
        fprintf(fileId,'\\multicolumn{%d}{c}{%s}', num_algos-1, vr_complexes{j});
        if j < num_shapes
            fprintf(fileId,'&\n');
        else
            fprintf(fileId,'\\\\\n');
        end
    end
    for j = 1:num_shapes
        for l = 1:num_algos
            if strcmp(algos{l}, 'pms') == 0
                fprintf(fileId,'& {%s} ', algos{l});
            end
        end
    end
    fprintf(fileId,'\\\\\n');
    fprintf(fileId,'\\midrule\n');

    % row 'it'
    for s = 1:num_samples
        fprintf(fileId,'%d', s);
        for j = 1:num_shapes
            for k = 1:num_algos
                if strcmp(algos{k}, 'pms') == 0
                    ratio = tensor_col_ops(j, s, 1)/tensor_col_ops(j, s, k);
                    fprintf(fileId,' & %1.2f', ratio);
                end
            end
        end
        fprintf(fileId,'\\\\\n');
    end

    fprintf(fileId,'\\bottomrule\n');
    fprintf(fileId,'\\end{tabular}\n');
    fprintf(fileId,'\\caption{Ratio of total column operations}\n');
    fprintf(fileId,'\\label{tab:ratio_cumsum_operations}\n');
    fprintf(fileId,'\\end{table*}\n');
    fprintf(fileId,'\\end{small}\n');

    fclose(fileId);
end

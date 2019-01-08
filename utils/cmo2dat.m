
function dat = cmo2dat(r, c, m)
    % Assume r,c,m are in row-column order
    % Assume if a column "j" is zero, we set c(idx)=j, r(idx)=-1
    % c(1) need not be 1
    % This assumption is needed by PMS
    dat = cell(m, 1);
    idx = 1;
    j_expected = 1;
    while idx < length(r)
        pos_start = idx;
        % Get segment of r for column j
        while idx <= length(r) && c(idx) == c(pos_start)
            idx = idx + 1;
        end
        pos_end = idx - 1;
        % Get dimension of the cell
        cell_dim = pos_end - pos_start;
        if cell_dim == 0 && r(pos_start) == -1
            cell_dim = -1;
        end
        % Get j
        j = c(pos_start);
        % Fill in missing columns with zeros
        while j > j_expected
            dat{j_expected} = [-1];
            j_expected = j_expected + 1;
        end
        % Insert to dat
        if cell_dim == -1
            dat{j} = [cell_dim];
        else
            dat{j} = [cell_dim; r(pos_start:pos_end)-1]; % PHAT uses 0-indexing
        end
        % j_expected
        j_expected = j_expected + 1;
    end
end


function [r, c, m] = stream2cmo(stream)
    import edu.stanford.math.plex4.*;
    ccs = streams.utility.StreamUtility.getCompressedBoundaryMatrix(stream);
    m = stream.getSize();
    r = double(cell2mat(cell(ccs.get(0).toArray()))); % rows
    c = double(cell2mat(cell(ccs.get(1).toArray()))); % cols
    % Sort to get in right format (do this in C later)
    [r, c] = sort_columns(r, c);
end

function [r, c] = sort_columns(r, c)
    % Assume c is a vector with non-decreasing values
    % Assume c(l) = k for all l in [i,j] and c(i-1) < k < c(j+1)
    % Then this function sorts each segment r(i:j)
    idx = 1;
    while (idx <= length(c))
        idx_left  = idx;
        idx_right = idx;
        while (idx_right <= length(c) & c(idx_right) == c(idx_left))
            idx_right = idx_right + 1;
        end
        idx_right = idx_right - 1;
        r(idx_left:idx_right) = sort(r(idx_left:idx_right));
        idx = idx_right + 1;
    end
end


function dat = cmo2dat(r, c, m)
    dat = cell(m, 1);
    j = 1;
    idx = 1;
    while j <= m
        pos_start = idx;
        while idx <= m && c(idx) == r(pos_start)
            idx = idx + 1;
        end
        pos_end = idx - 1;
        if j == c(pos_start)
            dat{j} = r(pos_start:pos_end);
        else
            while j < c(pos_start)
                dat{j} = 0 
                j = j + 1;
            end
        end
        j = j + 1;
    end
end


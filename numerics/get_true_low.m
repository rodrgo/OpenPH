function low_true = get_true_low(r, c, m, col_width)
    low_true = ph('pms', int32(r), int32(c), int32(m), int32(col_width), int32(zeros(1,m)));
    low_true = double(low_true);
end


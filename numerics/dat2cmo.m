
function [r, c, m] = dat2cmo(fname)
    % DAT file indexing starts at 0
    batch = 100;
    r = zeros(batch, 1);
    c = zeros(batch, 1);
    m = 0;
    idx = 0;
    % Opens file and reads line
    fid = fopen(fname);
    tline = fgetl(fid);
    while ischar(tline)
        s = strsplit(tline, ' ');
        if strcmp(s{1}, '#') == 0 & length(s{1}) > 0
            m = m + 1;
            if strcmp(s{1}, '0') == 0
                for i = 2:length(s)
                    idx = idx + 1;
                    r(idx) = str2num(s{i});
                    c(idx) = m-1;
                end
                if m >= length(r)
                    r = [r; zeros(batch, 1)];
                    c = [c; zeros(batch, 1)];
                end
            else
                m = m + 1;
                idx = idx + 1;
                r(idx) = m-1;
                c(idx) = m-1;
            end
        end
        tline = fgetl(fid);
    end
    if idx > 0
        r = r(1:idx);
        c = c(1:idx);
    else
        r = [];
        c = [];
    end
    fclose(fid);
end


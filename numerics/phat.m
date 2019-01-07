% PHAT.M
%   A wrapper for PHAT

function [low, t] = phat(r, c, m, PHAT_DIR, algorithm)
    algos = {'standard', 'twist', 'chunk', 'chunk_sequential', 'spectral_sequence', 'row'};
    if nargin < 5
        algorithm = 'twist';
    end
    assert(any(strcmp(algos, algorithm)));
    %"--standard, --twist, --chunk, --chunk_sequential, --spectral_sequence, --row 
    dat = cmo2dat(r, c, m);
    input_fname = print_dat(dat);
    %print_cmo(r, c, m);  % for debugging
    output_fname = generate_tempname();
    command = [PHAT_DIR '/phat ' '--' algorithm ' --ascii ' input_fname ' ' output_fname];
    t0 = tic;
    status = system(command);
    t = toc(t0);
    low = out2low(output_fname, m); 
    % Remove input_fname and output_fname
    delete(input_fname);
    delete(output_fname);
end

function print_cmo(r, c, m)
    fileID = fopen('cmo_tmp.debug','w');
    idx = 1;
    while idx <= length(c)
        fprintf(fileID, '%d, %d\n', c(idx), r(idx));
        idx = idx + 1;
    end
    fclose(fileID);
end

function fname = print_dat(dat)
    fname = generate_tempname();
    fileID = fopen(fname,'w');
    for i = 1:length(dat)
        for j = 1:length(dat{i}) 
            fprintf(fileID, '%d', dat{i}(j));
            if j < length(dat{i})
                fprintf(fileID, ' ');
            end
        end
        fprintf(fileID, '\n');
    end
    fclose(fileID);
end

function fname = generate_tempname()
    [~, tempfname] = fileparts(tempname);
    fname = [tempfname '.tmp'];
end

function low = out2low(fname, m)
    low = zeros(1, m);
    % Opens file and reads line
    fid = fopen(fname);
    tline = fgetl(fid);
    while ischar(tline)
        s = strsplit(tline, ' ');
        if length(s) == 2
            low(str2num(s{2})+1) = str2num(s{1})+1;
        end
        tline = fgetl(fid);
    end
    fclose(fid);
end

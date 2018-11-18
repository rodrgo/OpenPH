function tensor = populate_tensor(OUT, tensor, perc, i, k, l)
    % ess
    y = [OUT.ess(1:OUT.num_iters) 1];
    for pp = 1:length(perc.ess)
        tensor.ess(i, k, l, pp) = find(y - perc.ess(pp) >= 0 , 1, 'first');
    end

    % unred
    y = [OUT.err_redu(1:OUT.num_iters) 0];
    for pp = 1:length(perc.unred)
        tensor.unred(i, k, l, pp) = find(y - perc.unred(pp) <= 0 , 1, 'first');
    end

    % lone
    y = [OUT.err_lone(1:OUT.num_iters) 0];
    for pp = 1:length(perc.lone)
        tensor.lone(i, k, l, pp) = find(y - perc.lone(pp) <= 0 , 1, 'first');
    end

    % linf
    y = [OUT.err_linf(1:OUT.num_iters) 0];
    for pp = 1:length(perc.linf)
        tensor.linf(i, k, l, pp) = find(y - perc.linf(pp) <= 0 , 1, 'first');
    end
end

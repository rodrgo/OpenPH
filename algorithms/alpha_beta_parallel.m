% ALPHA_BETA_PARALLEL.M
% 
% Not designed for unreduced boundary matrix 

function [lows, t] = alpha_beta_parallel(D)
    t0 = tic;

    % Initialise persistence vectors
    D.init();
    D.record_iteration();

    % Set threading
    %   Can be either 'alpha_beta' or 'first_lows'
    D.set_threading('alpha_beta');

    verbose = false;
    if verbose
        fprintf('\nStart of algorithm:');
        D.print_ph_status();
    end

    % -------------
    % Phase 0
    % -------------

    % Get arglows and lowstars from alpha-beta reduce
    % These are initial pivots
    D.alpha_beta_reduce();

    while ~D.matrix_is_reduced()

        % -------------
        % Main iteration: Phase I
        % -------------

        % For each dimension, add to pivots first unreduced
        % column "j" such that low(j) == i
        D.mark_unique_unreduced_per_dimension();

        % Get known lowstars
        pivots_cell = D.get_pivots_cell();

        % -------------
        % Main iteration: Phase II
        % -------------

        D.reset_updated();

        % Parallel-reduce the matrix
        for l = 1:length(pivots_cell)
            pivot = pivots_cell{l};
            j0 = pivot.col;
            for j = pivot.neighbours
                % Given that "j0" is a true lowstar
                % only "j0" will hit "j"
                % Hence, only one thread will update
                % obj.low(j) at each iteration
                D.left_to_right(j0, j);
                D.mark_updated(j);
            end
        end

        % Try to infer nature of updated columns
        % We do this in two stages

        %   Check for alpha/beta information
        updated = D.get_updated();
        for j = updated
            D.alpha_beta_check(j);
        end

        % -------------
        % Main iteration: Phase III
        % -------------

        %   Run c8 check
        updated = D.get_updated();
        for j = updated
            D.c8_check(j);
        end

        if verbose
            D.print_ph_status();
            [1:length(D.low); D.simplex_dimensions; D.classes; D.low]
        end
        D.record_iteration();

    end

    % Set unmarked columns to essential
    D.set_unmarked();

    % Extract lows
    lows = D.low;
    t = toc(t0);
end


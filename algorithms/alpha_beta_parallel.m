% ALPHA_BETA_PARALLEL.M
% 
% Not designed for unreduced boundary matrix 

function [lows, t] = alpha_beta_parallel(D)
    t0 = tic;

    % Initialise persistence vectors
    D.init();
    D.record_iteration();

    D.create_candidate_pivots();

    verbose = false;
    if verbose
        fprintf('\nStart of algorithm:');
        D.print_ph_status();
    end

    while ~D.matrix_is_reduced()

        % Get arglows and lowstars from alpha-beta reduce
        D.alpha_beta_reduce();

        % Reduce first low
        % Sets as low the first unreduced column j such that
        % low(1:(j-1)) is reduced
        D.mark_first_low();

        % Get known lowstars
        lowstar_pivots = D.get_lowstar_pivots();

        % Reset candidate pivots here
        D.reset_candidate_pivots();
        D.reset_previous_lowstars();

        % Reduce
        can_update = ~all(cellfun(@(x) isempty(x.neighbours), lowstar_pivots));

        if can_update

            % -------------
            % Left-to-right operations
            % -------------

            % Parallel-reduce the matrix
            for l = 1:length(lowstar_pivots)
                pivot = lowstar_pivots{l};
                j0 = pivot.col;
                for j = pivot.neighbours
                    D.left_to_right(j0, j);
                end
            end

            % Update matrix invariants
            for l = 1:length(lowstar_pivots)
                pivot = lowstar_pivots{l};
                for j = pivot.neighbours
                    D.update_features(j);
                end
            end

            % -------------
            % Jared's suggestion 
            % -------------

            % Mark previous lowstars
            for l = 1:length(lowstar_pivots)
                pivot = lowstar_pivots{l};
                j0 = pivot.col;
                D.mark_previous_lowstars(j0);
            end

            % Update matrix invariants
            prev_lowstars = find(D.previous_lowstars);
            for l = 1:length(prev_lowstars)
                j0 = prev_lowstars(l);
                D.declare_previous_lowstars(j0);
            end

        else
            % set to unmarked was here before
        end
        if verbose
            fprintf('can_update = %d, iter = %d\n', can_update, iter);
            D.print_ph_status();
            %pause;
        end
        D.record_iteration();
    end

    % Set unmarked columns to essential
    D.set_unmarked();

    % Extract lows
    lows = D.low;
    t = toc(t0);
end


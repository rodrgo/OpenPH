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

    D.create_candidate_pivots();

    verbose = false;
    if verbose
        fprintf('\nStart of algorithm:');
        D.print_ph_status();
    end

    % Get arglows and lowstars from alpha-beta reduce
    D.alpha_beta_reduce();

    while ~D.matrix_is_reduced()

        % Reduce first low
        % Sets as low the first unreduced column j such that
        % low(1:(j-1)) is reduced
        % D.mark_first_low();

        % For each dimension, add to pivots first column
        % that has a low at a given position
        %D.mark_first_unreduced_per_dimension();
        D.mark_unique_unreduced_per_dimension();

        % Get known lowstars
        %   pivots = D.get_lowstar_pivots();
        pivots_cell = D.get_pivots_cell();

        % Reset candidate pivots here
        % D.reset_candidate_pivots();
        % D.reset_previous_lowstars();

        if length(pivots_cell) > 0

            % -------------
            % Left-to-right operations
            % -------------

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

            % Try to update information on
            % "unclassified pivots" and "updated columns"
            % We do this in two stages.
            % Update in two stages:
            %   First new alpha-beta and clearing
            %   Then other strategies 

            updated = D.get_updated();
            for j = updated
                D.alpha_beta_check(j);
            end

            updated = D.get_updated();
            for j = updated
                D.c8_check(j);
            end
             
%            % Try to find class of given markers
%            unreduced_pivots = find(obj.classes == 0 & obj.pivots > 0);
%            for j = unreduced_pivots
%                D.update_persistence_markers(j);
%            end
%
%            D.update_pivot_vectors();

%            % Update matrix invariants
%            for l = 1:length(pivots)
%                pivot = pivots{l};
%                for j = pivot.neighbours
%                    % Checks whether after parallel-reducing
%                    % we have enough information to declare
%                    % that "j" is reduced.
%                    D.update_features(j);
%                end
%            end

            % -------------
            % Jared's suggestion 
            % -------------

            % Mark lowstar strategy (move inside update features)
            % Mark previous lowstars
%            for l = 1:length(pivots)
%                pivot = pivots{l};
%                j0 = pivot.col;
%                D.mark_previous_lowstars(j0);
%            end

            % Update matrix invariants
%            prev_lowstars = find(D.previous_lowstars);
%            for l = 1:length(prev_lowstars)
%                j0 = prev_lowstars(l);
%                D.declare_previous_lowstars(j0);
%            end

        else
            % set to unmarked was here before
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


classdef ReductionParallel < BoundaryMatrix
    properties
        % Slow might be because we are not staging sufficient
        % candidate pivots or because some columns need a lot
        % of left to right oppertations to be reduced.
        % In this case it would be worth hitting more colums
        % at each iteration using some rule.

        % When updating in parallel we do not always have the
        % complete reduction information of columns 1:(j-1)
        % in thread "j".
        %
        % The parallel implementation first performs left-to-right
        % operations in parallel and then infers the nature of the
        % final column.
        % Sometimes there is enough information to tell whether
        % a given column is positive, negative, essential; however
        % some other times we can't guarantee that the column is
        % a true low, but we can benefit from using it to
        % reduce nevertheless. These are the candidate pivots. 
        candidate_pivots;
        previous_lowstars;

        % pivot_arglows(i) = min{j \in [m] : low(j) = i}
        % pivot_arglows(i) > 0 iff exists column with low(j) = i
        pivot_arglows; 
        pivot_arglows_needs_update;

        % Mark as true if column needs to be inspected for updating persistence markers
        columns_to_inspect;

        % pivot_lows(j) > 0 iff "j" is a pivot.
        % pivot_arglows(pivot_lows(j)) = j
        pivot_lows; 

        pivots_cell;

        updated;

        % Threading mode
        threading;

    end

    methods

        function obj = ReductionParallel(stream)
            obj = obj@BoundaryMatrix(stream);
            obj.create_pivot_vectors();
            obj.pivots_cell = cell(1, obj.m);
            obj.updated = false(1, obj.m);
            obj.threading = 'alpha_beta';
        end

        function set_threading(obj, threading)
            if ismember(threading, {'alpha_beta', 'first_low'})
                obj.threading = threading;
            else
                error('obj.threading not recognised');
            end
        end 

        function create_candidate_pivots(obj)
            obj.candidate_pivots = zeros(1, obj.m);
        end

        function create_pivot_vectors(obj)

            % pivot_arglows
            obj.pivot_arglows = zeros(obj.m, 1);
            for i = 1:obj.m
                pos = find(obj.low == i, 1, 'first');
                if ~isempty(pos)
                    obj.pivot_arglows(i) = pos;
                end
            end

            % pivots
            obj.pivot_lows = zeros(1, obj.m);
            for i = 1:obj.m
                if obj.pivot_arglows(i) ~= 0
                    obj.pivot_lows(obj.pivot_arglows(i)) = i;
                end
            end

            % pivot_arglows_needs_update
            obj.pivot_arglows_needs_update = zeros(obj.m, 1);

        end

        function mark_updated(obj, j)
            obj.updated(j) = true;
        end

        function mark_update_pivot_vectors(obj, j)
            % If a column "j" is in "neighbours"
            % it means that it was updated at this iteration.
            % Hence, there is a pivot column that operated
            % on it.
            % Need to update obj.pivots and obj.pivot_arglows
            % accordingly 
            %
            % obj.low(j) == 0
            %   We do nothing as this is clearly not a new pivot
            %
            % obj.low(j) > 0
            %   If obj.pivot_arglow(obj.low(j)) > j
            %       The current pivot for obj.low(j) is after j
            %       We don't know whether this condition holds for
            %       several columns.
            %       If so, we have to get the minimum, so mark column
            %       to recompute later.
            %   If obj.pivot_arglow(obj.low(j)) <= j
            %       obj.pivot_arglows is consistent, so we don't do
            %       anything in this case.
            if obj.pivot_arglows(obj.low(j)) > j
                obj.pivot_arglows_needs_update(obj.low(j)) = 1;
            end
        end

        function update_persistence_markers(obj, j)

            if obj.low(j) > 0
                is_lowstar = false;
                if obj.left(obj.low(j)) == j
                    is_lowstar = true;
                else
                    first_low = isempty(find(obj.low(1:(j-1)) == obj.low(j)));
                    if first_low
                        previous_reduced = all(obj.classes(1:(j-1)) ~= 0);
                        if previous_reduced
                            is_lowstar = true;
                        end
                    end
                end

                % -----------
                % If is_lowstar, do a twist clearing.
                % -----------
                if is_lowstar
                    obj.arglow(obj.low(j)) = j;
                    obj.mark_as_negative(j);
                    i = obj.low(j);
                    obj.clear_cols(i);
                end
            else
                obj.mark_as_positive(j);
            end

        end

        function update_pivot_vectors(obj)
            ind = find(obj.pivot_arglows_needs_update);
            for l = 1:length(ind)
                i = ind(l);
                % We do this because of collisions
                %   otherwise, do
                %   obj.pivot_arglows(obj.low(j)) = j;
                %
                % We don't need to do the search for all
                % rows. This can be made more efficient
                % by an atomic operation
                pos = find(obj.low == i, 1, 'first');
            end
            obj.pivot_arglows_needs_update(:) = 0;
        end

        function create_previous_lowstars(obj)
            obj.previous_lowstars = zeros(1, obj.m);
        end

        function reset_candidate_pivots(obj)
            obj.candidate_pivots(:) = 0;
        end

        function reset_previous_lowstars(obj)
            obj.previous_lowstars(:) = 0;
        end

        function mark_first_unreduced_per_dimension(obj)
            % Let cols_dim(d) = {columns of dimension "d"}
            % We look for j=min{col \in cols_dim(d) : col is unreduced}

            for d = 1:obj.complex_dimension
                % If this happens, then surely we haven't observed obj.(j)
                j = find(obj.classes == 0 & obj.simplex_dimensions == d, 1, 'first'); 
                if ~isempty(j) & obj.arglow(obj.low(j)) == 0
                    % We claim that isempty(j) == false
                    % if and only if the matrix is still unreduced.
                    % and dimension is not zero
                    assert(obj.low(j) > 0);
                    % If low(j) is positive and we have
                    % not observed it before, then this is
                    % a lowstar.
                    obj.arglow(obj.low(j)) = j;
                    obj.mark_as_negative(j);
                    i = obj.low(j);
                    % We call obj.clear_cols, though since
                    % i<j, i is already reduced.
                    obj.clear_cols(i);
                end
            end
        end

        function mark_first_low(obj)
            % We look for first column j such that
            % columns 1:(j-1) are reduced.
            % We update obj.classes so that if obj.classes ~=0
            % then column is reduced

            % Get first column j such that 1:(j-1) has
            % been reduced
            j = find(obj.classes == 0, 1, 'first'); 

            if ~isempty(j)
                % We claim that isempty(j) == false
                % if and only if the matrix is still unreduced.
                assert(obj.low(j) > 0);
                if obj.low(1:(j-1)) ~= obj.low(j)
                    % If low(j) is positive and we have
                    % not observed it before, then this is
                    % a lowstar.
                    obj.arglow(obj.low(j)) = j;
                    obj.mark_as_negative(j);
                    i = obj.low(j);
                    % We call obj.clear_cols, though since
                    % i<j, i is already reduced.
                    obj.clear_cols(i);
                end
            end
        end

        function mark_previous_lowstars(obj, j)

            % -------------
            % Look for lowstars in columns 1:(j-1)
            % -------------

            % For each "j", look for unreduced columns 1:(j-1)
            % and mark more obvious lowstars

            lows = obj.low(1:j);
            mm = max(lows);
            if mm > 0
                count = zeros(mm, 1);
                for l = 1:j
                    if lows(l) ~= 0
                        count(lows(l)) = count(lows(l)) + 1;
                    end
                end
                % These are the lowstars
                ind = find(count == 1 & obj.arglow(1:mm) == 0);
                % Mark these as lowstars
                for l = 1:length(ind)
                    col = find(lows == ind(l));
                    obj.previous_lowstars(col) = 1;
                end
            end

        end

        function declare_previous_lowstars(obj, j)
            obj.arglow(obj.low(j)) = j;
            obj.mark_as_negative(j);
            i = obj.low(j);
            obj.clear_cols(i);
        end

        function lowstars = find_nonzero_lowstar(obj)
            idx = obj.lowstar ~= -1;
            pos = find(idx);
            lowstars = obj.lowstar(pos);
        end

        function neighbours = get_row_neighbours(obj, i, j0)
            % Given i and j0 returns
            %   {j \in [m] : j > j0 & i == low(j)}
            neighbours = find(obj.low((j0+1):end) == i);
            if ~isempty(neighbours)
                neighbours = neighbours + j0;
                %obj.candidate_pivots(neighbours) = 0;
            end
        end

        function pivots_cell = get_pivots_cell(obj)
            if strcmp(obj.threading, 'alpha_beta')
                pivots_cell = obj.get_pivots_alphaBeta();
            elseif strcmp(obj.threading, 'firstLow')
                pivots_cell = obj.get_pivots_firstLow();
            else
                error('obj.threading not recognised');
            end
        end

        %function lowstar_pivots = get_lowstar_pivots(obj)
        function pivots_cell = get_pivots_alphaBeta(obj)
            % lowstar_pivots is a cell array of objects with
            %   pivot.col = column of lowstar
            %   pivot.row = row of lows
            %   pivot.neighbours = {j \in [m] : j > pivot.col & i == low(j)}

            % Reset pivots_cell
            obj.pivots_cell(:) = {[]};

            % Add pivots which are already a lowstar
            % and which have nonempty neighbours
            cols = obj.arglow(find(obj.arglow));
            cols = cols';
            assert(size(cols, 1) == 1);
            for j = cols
                i = obj.low(j); 
                neighbours = obj.get_row_neighbours(i, j);
                if ~isempty(neighbours)
                    pivot = [];
                    pivot.row = i;
                    pivot.col = j;
                    pivot.neighbours = neighbours;
                    obj.pivots_cell{j} = pivot;
                end
            end

            % Add candidate pivots
            % Not already 
            % Get indices of nonzero candidate_pivots
            cols = find(obj.candidate_pivots ~= 0); 
            for j = cols
                i = obj.low(j);
                neighbours = obj.get_row_neighbours(i, j);
                if ~isempty(neighbours)
                    pivot = [];
                    pivot.row = i;
                    pivot.col = j;
                    pivot.neighbours = neighbours;
                    obj.pivots_cell{j} = pivot;
                end
            end

            % Only return those that are not empty
            pivots_cell = obj.pivots_cell(cellfun(@(x) ~isempty(x), obj.pivots_cell));

        end

        function pivots_cell = get_pivots_firstLow(obj)
            % Recall
            %   pivot_arglows(i) = min{j \in [m] : low(j) = i}
            %   pivot_arglows(i) > 0 iff exists column with low(j) = i
            %   pivot_lows(j) > 0 iff "j" is a pivot.
            %   pivot_arglows(pivot_lows(j)) = j

            % Reset pivots_cell
            obj.pivots_cell(:) = {[]};

            % Get columns with nonzero pivot_lows
            cols = find(pivot_lows);
            for j = cols
                i = pivot_lows(j);
                neighbours = obj.get_row_neighbours(i, j);
                if ~isempty(neighbours)
                    pivot = [];
                    pivot.row = i;
                    pivot.col = j;
                    pivot.neighbours = neighbours;
                    obj.pivots_cell{j} = pivot; 
                end
            end 

            % Only return those that are not empty
            pivots_cell = obj.pivots_cell(cellfun(@(x) ~isempty(x), obj.pivots_cell));

        end

        function alpha_beta_check(obj, j)
            % Need to infer whether 'j' is a lowstar
            % This function checks:
            %   A. obj.low(j) == 0, then "j" is positive.
            %   B. If obj.left(obj.low(j)) == j, then "j" is reduced and negative.
            column_classified = false;
            is_lowstar = false;
            if obj.low(j) > 0
                if obj.left(obj.low(j)) == j
                    % If obj.left(obj.low(j)) == j, then
                    % obj.beta(obj.low(j)) == j, so the column is reduced
                    is_lowstar = true;
                    column_classified = true;
                end
            else
                obj.mark_as_positive(j);
                column_classified = true;
            end

            % If is_lowstar, do a twist clearing.
            if is_lowstar
                obj.arglow(obj.low(j)) = j;
                obj.mark_as_negative(j);
                i = obj.low(j);
                obj.clear_cols(i);
            end

            % Change updated vector for next stage
            if column_classified
                obj.updated(j) = false;
            end
        end

        function c8_check(obj, j)
            % Extended c8-check that also works with columns that do not necessarily
            % have obj.beta(j) == 0
            %
            % For any "j", let
            %   I = {i \in [j-1] : dim(i) == dim(j)-1}
            %
            % lowstar(j) \in S where
            %
            %    S = {i \in I : i \in [beta(j), low(j)]}                \setminus
            %        {i \in I : arglow(l) = i > 0}                      \setminus
            %        {i \in I : i = lowstar(l) for some l \in [j-1]}
            %
            % If obj.beta(j) == 0, then
            %   S = S \setminus {i \in I : beta(i) == 0}
            %
            % If obj.beta(j) >  0, then
            %   S = S \setminus {i \in I : beta(i) > 0}
            %

            dim_j = obj.simplex_dimensions(j);
            S = true(1, obj.m);
            % 
            S(1:(obj.beta(j)-1)) = false;
            S((obj.low(j)+1):end) = false;
            S(obj.simplex_dimensions ~= (dim_j-1)) = false;
            S(obj.arglow(obj.arglow > 0)) = false;
            S(obj.lowstar > 0) = false;
            S(obj.classes == -1) = false;

            if obj.beta(j) == 0
                % Column can be positive or negative
                S(obj.beta(j) > 0) = false;
                if nnz(S) == 0
                    obj.clear_cols(j);
                    fprintf('c8 success #1\n');
                end
            else
                % Column has to be negative
                S(obj.beta(j) > 0) = false;
                assert(nnz(S) > 0);
                if nnz(S) == 1
                    i = find(S);
                    obj.clear_cols(i);
                    fprintf('c8 success #2\n');
                end
            end

        end

        function ind = get_updated(obj)
            ind = find(obj.updated);
        end

        function update_features(obj, j)
            % Need to infer whether 'j' is a lowstar
            % Options:
            %   A. obj.low(j) == 0, then "j" is positive.
            %   B. If obj.left(obj.low(j)) == j, then "j" is reduced and negative.
            %   C. If all columns 1:(j-1) have been reduced and "j" and
            %      is the first column with given pivot
            %   D. (Using C8). Get candidate lows of unreduced columns and see if
            %       low for "j" doesn't overlap with these.
            if obj.low(j) > 0
                is_lowstar = false;
                if obj.left(obj.low(j)) == j
                    % If obj.left(obj.low(j)) == j, then
                    % obj.beta(obj.low(j)) == j, so the column is reduced
                    is_lowstar = true;
                else
                    c8_success = false;

                    if obj.beta(j) == 0 & false
                        % Use curiosity-8
                        % lowstar(j) \in {p \in [1, low(j)] : beta(p) == 0} \cup {0}
                        %loc = [find(obj.beta == 0 & 1:obj.m <= obj.low(j)), 0];

                        ind = 1:obj.low(j);
%                        obj.beta(ind) == 0;
%                        not obj.lowstar(ind) > 0;
%                        not obj.classes(ind) == -1;
                        loc = [0, find(obj.beta(ind) == 0 & obj.lowstar(ind) == 0 & obj.classes(ind) ~= -1)];
%                        %   Substract locations of known lowstars
%                        loc = setdiff(loc, obj.lowstar(obj.lowstar > 0));
%                        %   Substract locations of known negatives
%                        negative_cols = [find(obj.beta > 0), find(obj.classes == -1)];
%                        loc = setdiff(loc, negative_cols);
                        max_loc = max(loc);
                        if max_loc == 0
                            obj.clear_cols(j);
                            c8_success = true;
                        elseif false
                            % If all {l \in [j-1] : D(low(j), l) = 1}
                            % have low(l) > low(j), then this is a low  
                            % Check lefts and low 

                            % Once we have computed 'loc', we check look at x=left(low(j)) and at y=arglow(low(j)).
                            % If y exists, we look at the low(y, -1) and see if this belongs to loc. If it doesn't, then this is a lowstar.
                            % If this doesn't give information, look at 

                        elseif length(loc) == 1
                            if max_loc == obj.low(j)
                                is_lowstar = true;
                                c8_success = true;
                            else
                                % In this case the lowstar is somewhere
                                % above low(j), but we don't know what the
                                % final sparsity pattern of the matrix is
                                % We can't consider this columns a pivot.
                            end
                        elseif max_loc == obj.low(j) 
                            %obj.candidate_pivots(j) = obj.low(j);
                            obj.candidate_pivots(j) = 1;
                        end
                    end

                    if ~c8_success
                        % Only pursue this strategy if c8_success
                        % Check if "j" is the first column that has obj.low(j)
                        first_low = isempty(find(obj.low(1:(j-1)) == obj.low(j)));
                        if first_low
                            % If it is, check if all columns 1:(j-1) have
                            % already been reduced
                            previous_reduced = all(obj.classes(1:(j-1)) ~= 0);
                            if previous_reduced
                                is_lowstar = true;
                            else
                                % We need to mark those columns which
                                % we can't guarantee are lowstars due
                                % to the parallel implementation, but 
                                % which are the first lows observed
                                %
                                % This object does not update
                                % automatically
                                % obj.candidate_pivots(j) = obj.low(j);
                                obj.candidate_pivots(j) = 1;
                            end
                        end
                    end
                end

                % -----------
                % If is_lowstar, do a twist clearing.
                % -----------
                if is_lowstar
                    obj.arglow(obj.low(j)) = j;
                    obj.mark_as_negative(j);
                    i = obj.low(j);
                    obj.clear_cols(i);
                end
            else
                obj.mark_as_positive(j);
            end

        end


    end
end


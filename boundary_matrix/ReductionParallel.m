classdef ReductionParallel < BoundaryMatrix
    properties
        % When updating in parallel we do not always have the
        % complete reduction information of columns 1:(j-1)
        % in thread "j".
        
        % candidate_pivots is an auxiliary indicator vector to mark more pivots
        % than alpha/beta's
        candidate_pivots;

        % pivot_arglows(i) = min{j \in [m] : low(j) = i}
        % pivot_arglows(i) > 0 iff exists column with low(j) = i
        pivot_arglows; 

        % pivot_lows(j) > 0 iff "j" is a pivot.
        % pivot_arglows(pivot_lows(j)) = j
        pivot_lows; 

        pivots_cell;

        updated;

        % Dimensions index
        dimensions_index;

        % Threading mode
        threading;

    end

    methods

        % ===================
        % Constructor & algo properties
        % ===================

        function obj = ReductionParallel(stream)
            obj = obj@BoundaryMatrix(stream);
            obj.create_pivot_vectors();
            obj.pivots_cell = cell(1, obj.m);
            obj.updated = false(1, obj.m);
            obj.threading = 'alpha_beta';

            % Dimensions index
            obj.dimensions_index = cell(1, obj.complex_dimension);
            for d = 1:obj.complex_dimension
                obj.dimensions_index{d} = find(obj.simplex_dimensions == d);
            end
        end

        function set_threading(obj, threading)
            if ismember(threading, {'alpha_beta', 'first_low'})
                obj.threading = threading;
            else
                error('obj.threading not recognised');
            end
        end 

        % ===================
        % Create
        % ===================

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
        end

        % ===================
        % Get
        % ===================

        function reset_candidate_pivots(obj)
            obj.candidate_pivots(:) = 0;
        end

        function ind = get_updated(obj)
            ind = find(obj.updated);
        end

        function reset_updated(obj)
            obj.updated(:) = false;
        end

        function neighbours = get_row_neighbours(obj, i, j0)
            % Given i and j0 returns
            %   {j \in [m] : j > j0 & i == low(j)}
            neighbours = find(obj.low((j0+1):end) == i);
            if ~isempty(neighbours)
                neighbours = neighbours + j0;
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

        % ===================
        % Mark
        % ===================

        function mark_updated(obj, j)
            obj.updated(j) = true;
        end

        function mark_unique_unreduced_per_dimension(obj)
            % Let cols_dim(d) = {columns of dimension "d"}
            % We look for j=min{col \in cols_dim(d) : col is unreduced}

            % For each dimension
            visited = false(obj.m, 1);
            for d = obj.complex_dimension:-1:1
                % If this happens, then surely we haven't observed obj.(j)
                % j = find(obj.classes == 0 & obj.simplex_dimensions == d, 1, 'first'); 
                % if ~isempty(j) & obj.arglow(obj.low(j)) == 0
                visited(:) = false;
                ind = obj.dimensions_index{d};
                ceiling = 0;
                for j = ind
                    if obj.low(j) > 0
                        if ~visited(obj.low(j))
                            if obj.classes(j) == 0 && obj.low(j) > ceiling
                                % mark as lowstar
                                obj.arglow(obj.low(j)) = j;
                                obj.mark_as_negative(j);
                                i = obj.low(j);
                                % We call obj.clear_cols, though since
                                % i<j, i is already reduced.
                                obj.clear_cols(i);
                            end
                        else
                            ceiling = max(ceiling, obj.low(j));
                        end
                        % mark as visited
                        visited(obj.low(j)) = true;
                    end
                end
            end
        end

        % ===================
        % Check column nature
        % ===================

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
            % Extended c8-check that also works
            % with columns that do not necessarily
            % have obj.beta(j) == 0
            %
            % For any "j", let
            %   I = {i \in [m] : dim(i) == dim(j)-1}
            %
            % lowstar(j) \in S where
            %
            %    S = {i \in I : i \in [beta(j), low(j)]}                \setminus
            %        {i \in I : arglow(l) = i > 0}                      \setminus
            %        {i \in I : i = lowstar(l) for some l \in [j-1]}
            %   S = S \setminus {i \in I : beta(i) > 0}

            dim_j = obj.simplex_dimensions(j);
            S = true(1, obj.m);

            % Global-palette

            S(obj.simplex_dimensions ~= (dim_j-1)) = false;
            S(obj.beta > 0) = false;
            S(obj.lowstar > 0) = false;
            S(obj.arglow(obj.arglow > 0)) = false;
            S(obj.classes == -1) = false;

            % Local-palette

            S(1:(obj.beta(j)-1)) = false;
            S((obj.low(j)+1):end) = false;

            if obj.beta(j) == 0
                % Column can be positive or negative
                if nnz(S) == 0
                    obj.clear_cols(j);
                    %fprintf('success #1: %d\n', j);
                else
                    if nnz(S) == 1
                        %fprintf('almost: %d\n', j);
                    end
                end
            else
                % Column has to be negative
                assert(nnz(S) > 0);
                if nnz(S) == 1
                    i = find(S);
                    obj.clear_cols(i);
                    %fprintf('success #2: %d\n', i);
                end
            end

        end

    end
end


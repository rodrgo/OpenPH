classdef ReductionParallel < BoundaryMatrix
    properties

        % Candidate pivots (for alpha-beta-parallel reduction)
        candidate_pivots;
        previous_lowstars;

    end

    methods

        function obj = ReductionParallel(stream, homology_mode)
            obj = obj@BoundaryMatrix(stream, homology_mode);
        end

        function create_candidate_pivots(obj)
            obj.candidate_pivots = zeros(1, obj.m);
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

        function mark_first_low(obj)
            % We look for first non-essential column j such that
            % columns 1:(j-1) are unreduced.
            % We update obj.classes so that if obj.classes ~=0
            % then column is reduced

            j = find(obj.classes == 0, 1, 'first'); 

            if ~isempty(j)
                assert(obj.low(j) > 0);
                if obj.low(1:(j-1)) ~= obj.low(j)
                    obj.arglow(obj.low(j)) = j;
                    obj.mark_as_negative(j);
                    i = obj.low(j);
                    obj.clear_cols(i);
                end
            end

        end

        function update_features(obj, j)
            if obj.low(j) > 0
                is_lowstar = false;
                if obj.left(obj.low(j)) == j
                    is_lowstar = true;
                else
                    % Check if all columns 1:(j-1) have been reduced
                    previous_reduced = all(obj.classes(1:(j-1)) ~= 0);
                    first_low = isempty(find(obj.low(1:(j-1)) == obj.low(j)));
                    % If low is not a left, we can only guarantee that
                    % it is a lowstar if the previous columns
                    % have already been reduced (since we are
                    % updating in parallel)
                    if first_low
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
                            obj.candidate_pivots(j) = obj.low(j);
                        end
                    end
                end
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

        function neighbours = get_row_neighbours(obj, i, j0, mode)
            if strcmp(mode, 'all')
                neighbours = find(obj.matrix(i, (j0+1):end) == 1);
            elseif strcmp(mode, 'lowstar')
                neighbours = find(obj.low((j0+1):end) == i);
            else
                error('mode not recognised');
            end

            if ~isempty(neighbours)
                neighbours = neighbours + j0;
            end
        end

        function lowstar_pivots = get_lowstar_pivots(obj)
            % lowstar_pivots is a cell array of objects with
            %   pivot.col = column of lowstar
            %   pivot.row = row of lows
            %   pivot.neighbours = {j \in [m] : D_{i,j} = 1 for j >= pivot.col and i == pivot.row}

            lowstars = find(obj.arglow);
            lowstar_cols = obj.arglow(lowstars);
            lowstar_pivots = cell(1, length(lowstars));
            for l = 1:length(lowstar_pivots)
                i = lowstars(l);
                j = lowstar_cols(l);
                pivot = [];
                pivot.row = i;
                pivot.col = j;
                pivot.neighbours = obj.get_row_neighbours(i, j, 'lowstar');
                lowstar_pivots{l} = pivot; 
            end 

            % Add candidate pivots
            % Not already 
            % Get indices of nonzero candidate_pivots
            ind = find(obj.candidate_pivots ~= 0); 
            candidate_pivots = cell(1, length(ind));
            for l = 1:length(ind)
                j = ind(l);
                i = obj.candidate_pivots(j);
                pivot = [];
                pivot.row = i;
                pivot.col = j;
                pivot.neighbours = obj.get_row_neighbours(i, j, 'lowstar');
                candidate_pivots{l} = pivot;
            end

            % Add first lows

            lowstar_pivots = [lowstar_pivots, candidate_pivots];

        end

    end
end


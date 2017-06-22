
classdef ReductionC8 < BoundaryMatrix & handle

    methods
        function obj = ReductionC8(stream, homology_mode)
            obj = obj@BoundaryMatrix(stream, homology_mode);
        end

        % ===================
        % Reduce
        % ===================

        function curiosity_8_clearing(obj)
            % This function assumes that:
            %   obj.init() has been called
            %   obj.alpha_beta_reduce() has been done
            % We do not perform it again for performance reasons

            % Depending on the homology_mode we need to set
            % a dim_thresh parameter
            if strcmp(obj.homology_mode, 'reduced')
                dim_thresh = 1;
            elseif strcmp(obj.homology_mode, 'unreduced')
                dim_thresh = 2;
            end

            % We only consider columns "j" such that:
            %   1. Have beta(j) == 0,
            %   2. Have not been identified as positives by alpha_beta_reduction.
            beta_zero_cols = find((obj.beta == 0) & (obj.classes ~= 1));

            % We only consider rows "i" such that:
            %   1. After doing alpha_beta reduction arglow(i) == 0.
            %       That is, columns for which low* has not been identified.
            %   2. obj.beta(i) == 0 (as stated by the curiosity).
            rows = (obj.arglow == 0)' & (obj.beta == 0);
            for j = beta_zero_cols
                % Consider "i" st dim(i) = dim(j) - 1 
                dim_j = obj.initial_dimensions(j);
                if dim_j > dim_thresh
                    ind_bd = obj.initial_dimensions == dim_j - 1;
                    row_range = rows & ind_bd;
                else
                    row_range = rows;
                end
                % Clear column if row_range is empty
                % (since therefore low* has to be zero, see curiosity)
                if ~any(row_range(1:obj.alpha(j)))
                    obj.clear_cols(j);
                end
            end
        end

        % ===================
        % Visualisation
        % ===================

        function mask_cell = get_curiosity_8_mask(obj)
            % Do alpha beta reduction
            obj.init();
            obj.create_alpha();
            obj.create_beta();
            neg = (obj.alpha == obj.beta) & (obj.beta > 0);
            obj.mark_as_negative(neg);
            pos_pairs = obj.beta(neg);
            obj.clear_cols(pos_pairs);
            obj.arglow(pos_pairs) = find(neg);

            % Get an indicator vector with location of lowstars
            lowstar_ind = zeros(size(neg));
            lowstar_ind(pos_pairs) = true;

            % Create vector of indicator classes
            pos_markers = zeros(1, obj.m);
            pos_markers(pos_pairs) = 1;

            % Locations
            locs = cell(obj.m, 1);
            n = 0;
            for j = 1:obj.m
                if obj.beta(j) == 0 && obj.classes(j) ~= 1
                    % Get \{p \in [1, \alpha_j] : \beta_p = 0\}
                    row_range = (1:obj.m <= obj.alpha(j)) & obj.beta == 0;
                    % Consider only i st dim(i) = dim(j) - 1 
                    dim_j = obj.initial_dimensions(j);
                    % The dimension filter changes depending on
                    % whether the homology_mode is 'reduced' or 'unreduced'
                    if strcmp(obj.homology_mode, 'reduced')
                        dim_thresh = 1;
                    elseif strcmp(obj.homology_mode, 'unreduced')
                        dim_thresh = 2;
                    end
                    if dim_j > dim_thresh
                        ind_bd = obj.initial_dimensions == dim_j - 1;
                        row_range = row_range & ind_bd;
                    end
                    % Consier only those that are not already lowstars
                    row_range = row_range & (~lowstar_ind);
                    % Get locations
                    loc = find(row_range);
                    % If there are no locations, it must be a zero
                    % and this zero is not given by alpha-beta reduction
                    % so mark it explicitly
                    if isempty(loc)
                        pos_markers(j) = 2;
                    end
                else
                    loc = [];
                end
                locs{j} = loc;
                n = n + length(loc);
            end

            % Build these indices as a matrix
            rows = zeros(n, 1);
            cols = zeros(n, 1);
            l = 1;
            for j = 1:obj.m
                loc = locs{j};
                if ~isempty(loc)
                    idx = l:(l+length(loc)-1);
                    rows(idx) = loc;
                    cols(idx) = j;
                    l = l + length(loc);
                end
            end
            mask = sparse(rows, cols, ones(size(rows)), obj.m, obj.m);
            mask_cell_1 = {mask, 'c8 candidate', 'dk'};

            % Create two more masks with pos_markers data
            cols = find(pos_markers == 1);
            mask_pos_ab = sparse(ones(size(cols)), cols, ones(size(cols)), obj.m, obj.m); 
            mask_cell_2 = {mask_pos_ab, 'alpha/beta => positive', 'c<'};

            cols = find(pos_markers == 2);
            mask_pos_c8 = sparse(ones(size(cols)), cols, ones(size(cols)), obj.m, obj.m); 
            mask_cell_3 = {mask_pos_c8, 'c8 => positive', 'm>'};

            % Merge all masks
            mask_cell = {mask_cell_1, mask_cell_2, mask_cell_3};
        end

    end
end


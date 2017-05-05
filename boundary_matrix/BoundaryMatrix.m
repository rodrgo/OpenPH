classdef BoundaryMatrix < handle
    properties
        % Matrix homology mode
        % homology_mode is 'reduced' or 'unreduced'
        %
        % 'reduced': Add dummy simplex so that 0-Betti number
        %   counts components.
        % 'unreduced': Does not add dummy simplex.
        homology_mode;

        % Matrix properties
        m;
        matrix;
        complex_dimension;
        initial_dimensions;

        % Persistence information
        %
        % low(j) = max{i \in [m] : D_{i,j} = 1}
        % low converges to low* at the end of reduction.
        %
        % classes(j) = class of column j after reduction.
        %   + 1 iff j is positive
        %   - 1 iff j is negative
        %   Inf iff j is essential
        %   0   iff j is unreduced
        %
        % arglow(i) = \{j \in [m] : low*(j) = i\}
        %
        % These vectors must be maintained whenever the
        % the state of BoundaryMatrix is modified
        low;
        classes;
        arglow;

        % Persistence sub-information
        left;
        alpha;
        beta;
        rho;

        % Booleans
        has_low;
        has_arglow;
        has_classes;
        has_alpha;
        has_left;
        has_beta;
        has_rho;
    end

    methods
        function obj = BoundaryMatrix(stream, homology_mode)
            if nargin > 0

                % Default homology_mode is plain
                if nargin == 1
                    obj.homology_mode = 'unreduced';
                else
                    obj.homology_mode = homology_mode;
                end

                % Get boundary matrix for given stream from Javaplex
                import edu.stanford.math.plex4.*;
                if strcmp(obj.homology_mode, 'reduced')
                    ccs = streams.utility.StreamUtility.getCompressedBoundaryMatrixRedHom(stream);
                    obj.m = stream.getSize() + 1; % Add dummy simplex
                elseif strcmp(obj.homology_mode, 'unreduced')
                    ccs = streams.utility.StreamUtility.getCompressedBoundaryMatrix(stream);
                    obj.m = stream.getSize();
                else
                    assert(false, 'homology_mode not recognised');
                end

                % Create sparse boundary matrix
                rows = double(cell2mat(cell(ccs.get(0).toArray())));
                cols = double(cell2mat(cell(ccs.get(1).toArray())));
                vals = ones(size(rows));
                obj.matrix = sparse(rows, cols, vals, obj.m, obj.m);
                assert(istriu(obj.matrix));

                % Get dimension of the complex
                obj.initial_dimensions = sum(obj.matrix);
                obj.complex_dimension = max(obj.initial_dimensions);

                % Set vectors 
                obj.has_low = false;
                obj.has_arglow = false;
                obj.has_classes = false;
                obj.has_alpha = false;
                obj.has_left = false;
                obj.has_beta = false;
                obj.has_rho = false;

            end
        end

        %%%%%%%%%%%%%%%%%%%%% 
        % Query information
        %%%%%%%%%%%%%%%%%%%%% 

        function d_simplices = get_d_simplices(obj, d)
            d_simplices = find(obj.initial_dimensions == d);
        end

        function b = is_reduced(obj, j)
            b = ~(obj.low(j) > 0 && obj.arglow(obj.low(j)) ~= 0);
        end

        function unreduced_cols = get_unreduced_cols(obj)
            unreduced_cols = find(obj.classes == 0);
        end

        function ess = get_essential(obj)
            ess = find(obj.classes == Inf);
        end

        function twist_cols = get_twist_cols_unreduced(obj)
            % Get inidicator vector of unreduced columns
            unreduced_cols_idx = obj.classes == 0;

            % Copy initial_dimensions vector and mark with 0
            % all the columns that have already been reduced
            initial_dimensions_aux = obj.initial_dimensions;
            initial_dimensions_aux(~unreduced_cols_idx) = 0;

            twist_cols = zeros(1, nnz(initial_dimensions_aux));
            
            start = 1;
            complex_dim = obj.complex_dimension;
            for d = complex_dim:-1:1 
                d_simplices_idx = initial_dimensions_aux == d;
                stop = start + nnz(d_simplices_idx) - 1;
                twist_cols(start:stop) = find(d_simplices_idx);
                start = stop + 1;
            end

        end

        %%%%%%%%%%%%%%%%%%%%% 
        % Modify matrix
        %%%%%%%%%%%%%%%%%%%%% 
        
        function as_dense(obj)
            obj.matrix = full(obj.matrix);
        end

        function clear_cols(obj, idx)
            obj.matrix(:, idx) = 0;
            obj.low(idx) = 0;
            obj.mark_as_positive(idx);
        end

        function reduce_col_twist(obj, j)
            obj.reduce_col(j);
            if obj.low(j) > 0
                i = obj.low(j);
                obj.clear_cols(i);
            end
        end

        function reduce_col(obj, j)
            while (obj.low(j) > 0 && (obj.arglow(obj.low(j)) ~=0))
                j0 = obj.arglow(obj.low(j));
                obj.matrix(:, j) = mod(obj.matrix(:, j) + obj.matrix(:, j0), 2);
                obj.low(j) = obj.get_low(j);
            end
            obj.low(j) = obj.get_low(j);
            if obj.low(j) > 0
                obj.arglow(obj.low(j)) = j;
                obj.mark_as_negative(j);
            else
                obj.mark_as_positive(j);
            end
        end

        function alpha_beta_reduce(obj)
            % Create relevant vectors
            obj.create_alpha();
            obj.create_beta();

            % Pick obvious lows from alpha/beta pairs
            neg = (obj.alpha == obj.beta) & (obj.beta > 0);
            obj.mark_as_negative(neg);

            % Clear and mark associated positive columns
            pos_pairs = obj.beta(neg);
            obj.clear_cols(pos_pairs);

            % Mark arglows
            obj.arglow(pos_pairs) = find(neg);
        end

        function rho_clearing(obj)

            % Define range of action for rho_clearing
            max_rho_zero = find(obj.rho == 0, 1, 'last');

            if ~isempty(max_rho_zero)
                ind = 1:obj.m <= max_rho_zero;

                % Clear positive columns and mark them as positive
                pos_idx = (obj.beta == 0) & ind;
                obj.clear_cols(pos_idx);

                % Remaining indices in ind are negative
                % Use local alpha-beta reduction on remaining indices
                % to find lows.
                % TO DO:
                % If we are doing alpha-beta reduction this is redundant
                obj.create_alpha();
                neg = (obj.alpha == obj.beta) & (obj.beta > 0) & ind;
                obj.mark_as_negative(neg);
                obj.arglow(obj.low(find(neg))) = find(neg);
            end

            % Experiment: reduce with 'kappa' curve
            % Explore kappa curve for morozov
            if false
                kappa = zeros(1, obj.m);
                for j = 1:obj.m
                    if obj.low(j) > 0
                        l = obj.left(obj.low(j));
                        if (obj.low(l) == obj.low(j)) && (l < j)
                            kappa(j) = obj.low(j);
                        end
                    end
                end

                for j = 2:obj.m
                    kappa(j) = max(kappa(j-1), kappa(j));
                end

                unreduced_cols = obj.get_unreduced_cols();
                for j = unreduced_cols
                    fprintf('j=%d, low=%d, kappa=%d\n',j,obj.low(j),kappa(j));
                    if (obj.arglow(obj.low(j)) == 0) && (obj.low(j) > kappa(j))
                        obj.arglow(obj.low(j)) = j;
                        obj.mark_as_negative(j);
                    end
                end
                fprintf('Size kappa-red: %d\n', length(unreduced_cols));
            end

        end

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

        %%%%%%%%%%%%%%%%%%%%% 
        % Persistence information vectors
        %%%%%%%%%%%%%%%%%%%%% 

        % This function wraps the creation of
        % low, arglow, classes
        function init(obj)
            obj.create_low();
            obj.create_arglow();
            obj.create_classes();
        end

        function l = get_low(obj, j)
            if any(obj.matrix(:, j))
                l = find(obj.matrix(:, j), 1, 'last');
            else
                l = 0;
            end
        end

        function l = get_left(obj, i)
            if any(obj.matrix(i, :))
                l = find(obj.matrix(i, :), 1, 'first');
            else
                l = 0;
            end
        end

        function create_classes(obj)
            if ~obj.has_classes
                obj.classes = zeros(1, obj.m);
                obj.has_classes = true;
            end
        end

        function create_low(obj)
            if ~obj.has_low
                % Create low
                obj.low = zeros(1, obj.m);
                for j = 1:obj.m
                    obj.low(j) = obj.get_low(j);
                end
                % Mark low as created
                obj.has_low = true;
            end
        end

        function create_arglow(obj)
            if ~obj.has_arglow
                % Create arglow
                obj.arglow = zeros(obj.m, 1);
                % Mark arglow as created
                obj.has_arglow = true;
            end
        end


        function create_left(obj)
            if ~obj.has_left
                % Create left
                obj.left = zeros(obj.m, 1);
                for i = 1:obj.m
                    obj.left(i) = obj.get_left(i);
                end
                % Mark left as created
                obj.has_left = true;
            end
        end

        function create_alpha(obj)
            if ~obj.has_alpha
                % Requirements for alpha
                if ~obj.has_low
                    obj.create_low();
                end
                % Create alpha
                obj.alpha = zeros(1, obj.m);
                for j = 1:obj.m
                  obj.alpha(j) = obj.low(j);
                end
                % Mark alpha as created
                obj.has_alpha = true;
            end
        end

        function create_rho(obj)
            if ~obj.has_rho
                % Requirements for rho
                if ~obj.has_beta
                    obj.create_beta();
                end
                % Create rho
                obj.rho = zeros(1, obj.m);
                if strcmp(obj.homology_mode, 'reduced')
                    obj.rho(1) = (-1)^(0 + (obj.beta(1) == 0));
                elseif strcmp(obj.homology_mode, 'unreduced')
                    obj.rho(1) = 0;
                end

                for j = 2:obj.m
                    if obj.beta(j) > 0
                            obj.rho(j) = obj.rho(j-1) + 1;
                    else
                            obj.rho(j) = obj.rho(j-1) - 1;
                    end
                end
                % Mark rho as created
                obj.has_rho = true;
            end
        end

        function create_beta(obj)
            if ~obj.has_beta
                % Requirements for beta
                if ~obj.has_left
                    obj.create_left();
                end
                % Create beta
                obj.beta = zeros(1, obj.m);
                for j = 1:obj.m
                    lefts_in_col_j = find(obj.left == j);
                    if ~isempty(lefts_in_col_j)
                        obj.beta(j) = max(lefts_in_col_j);
                    end
                end
                % Mark beta as created
                obj.has_beta = true;
            end
        end

        %%%%%%%%%%%%%%%%%%%%% 
        % Persistence markers
        %%%%%%%%%%%%%%%%%%%%% 

        function ph_info = get_persistence_info(obj)
            neg = reshape(find(obj.low), [], 1);
            % Find pairs
            pairs = [reshape(obj.low(neg), [], 1) neg];
            % Get indicator vector of paired
            paired = zeros(obj.m, 1);
            paired(pairs(:, 1)) = 1;
            paired(pairs(:, 2)) = 1;
            % Get essentials
            essentials = find(paired == 0);
            ph_info = {pairs, essentials};
        end

        function mark_classes(obj)
            obj.mark_as_positive(obj.low == 0);
            obj.mark_as_negative(obj.low > 0);
            % Essentials are positive columns that are not
            % paired with any negative column.
            obj.mark_as_essential((obj.classes == 1)' & obj.arglow == 0);
        end

        % Positives can either be paired or essential
        function mark_as_positive(obj, idx)
            obj.classes(idx) = +1;
        end

        function mark_as_negative(obj, idx)
            obj.classes(idx) = -1;
        end

        function mark_as_essential(obj, idx)
            obj.classes(idx) = Inf;
        end

        %%%%%%%%%%%%%%%%%%%%% 
        % Visualisation masks
        %%%%%%%%%%%%%%%%%%%%% 

        function mask_cell = get_alpha_mask(obj)
            obj.create_alpha();
            idx = find(obj.alpha);
            mask = sparse(obj.alpha(idx), idx, ones(size(idx)), obj.m, obj.m);
            mask_cell = {mask, 'alpha', 'xr'};
        end

        function mask_cell = get_beta_mask(obj)
            obj.create_beta();
            idx = find(obj.beta);
            mask = sparse(obj.beta(idx), idx, ones(size(idx)), obj.m, obj.m);
            mask_cell = {mask, 'beta', '+b'};
        end

        function mask_cell = get_left_mask(obj)
            obj.create_left();
            idx = find(obj.left);
            mask = sparse(idx, obj.left(idx), ones(size(idx)), obj.m, obj.m);
            mask_cell = {mask, 'left', '.g'};
        end

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


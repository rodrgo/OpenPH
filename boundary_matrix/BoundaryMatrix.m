classdef BoundaryMatrix < handle
    properties

        % ---------------
        % Homology mode
        % ---------------
        %   'reduced': Add dummy simplex so that 0-Betti number
        %              counts components.
        %   'unreduced': Does not add dummy simplex.
        homology_mode;
        
        % ---------------
        % Matrix properties
        % ---------------

        % The boundary matrix
        matrix;

        % Number of columns/rows in matrix
        m;

        % Dimension of the complex
        complex_dimension;

        % Initial dimensions of columns
        initial_dimensions;

        % ---------------
        % Persistence information
        %
        %   These vectors must be maintained whenever the
        %   the state of BoundaryMatrix is modified
        % ---------------

        % low(j) = max{i \in [m] : matrix(i,j) = 1}
        %   low converges to low* at the end of reduction.
        low;
        has_low;

        % lowstar(j) is the lowstar of column j
        lowstar;
        has_lowstar;

        % arglow(i) = {j \in [m] : low*(j) = i}
        arglow;
        has_arglow;

        % classes(j) = class of column j after reduction.
        %   + 1 iff j is positive
        %   - 1 iff j is negative
        %   Inf iff j is essential
        %   0   iff j is unreduced
        classes;
        has_classes;

        % left(i) = min{j \in [m] : matrix(i, j) = 1}
        left;
        has_left;

        % alpha(j) = max{i \in [m] : matrix(i, j) = 1}
        alpha;
        has_alpha;

        % beta(j) = max{i \in [m] : left(i) = j}
        beta;
        has_beta;

        % ---------------
        % Convergence Metrics
        % ---------------

        metrics;


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

                % Booleans
                obj.has_low = false;
                obj.has_lowstar = false;
                obj.has_arglow = false;
                obj.has_classes = false;
                obj.has_alpha = false;
                obj.has_left = false;
                obj.has_beta = false;

                % ----------------
                % obj.metrics
                % ----------------

                obj.metrics = [];
                obj.metrics.next_iter = 1;
                obj.metrics.iters = obj.metrics.next_iter - 1;

                % iters.num_column_adds(k) = number of column adds at iteration k
                obj.metrics.num_column_adds= zeros(1, obj.m);

                % iters.num_entry_adds(k) = number of entries that change from
                %   0->1 or 1->0 at each iteration
                obj.metrics.num_entry_adds = zeros(1, obj.m);

                % iters.percentage_reduced(k) = percentage of reduced columns at
                %   each iteration
                obj.metrics.percentage_unreduced = ones(1, obj.m);

            end
        end

        % ===================
        % as_dense
        % ===================
        
        function as_dense(obj)
            obj.matrix = full(obj.matrix);
        end

        % ===================
        % Metrics
        % ===================

        function record_column_add(obj, left, right)
            pos = obj.metrics.next_iter;
            obj.metrics.num_column_adds(pos) = obj.metrics.num_column_adds(pos) + 1;
            num_flips = obj.m - nnz(obj.matrix(:, left) ~= 0 & obj.matrix(:, right) ~= 0);
            obj.metrics.num_entry_adds(pos) = obj.metrics.num_entry_adds(pos) + num_flips;
        end

        function record_iteration(obj)
            pos = obj.metrics.next_iter;
            obj.metrics.percentage_unreduced(pos) = 1 - nnz(obj.classes)/obj.m;
            obj.metrics.iters = obj.metrics.next_iter;
            obj.metrics.next_iter = obj.metrics.next_iter + 1;
            if obj.metrics.next_iter > length(obj.metrics.num_column_adds)
                pos = obj.metrics.next_iter;
                obj.metrics.num_column_adds(pos) = 0;
                obj.metrics.num_entry_adds(pos) = 0;
                obj.metrics.percentage_unreduced(pos) = 1;
            end
        end

        % ===================
        % Booleans
        % ===================

        function b = is_reduced(obj, j)
            b = ~(obj.low(j) > 0 && obj.arglow(obj.low(j)) ~= 0);
        end

        function b = matrix_is_reduced(obj)
            lows = obj.low(obj.low > 0);
            b = length(unique(lows)) == length(lows);
        end

        function b = is_negative(obj, j)
            if obj.has_beta && obj.has_classes
                b = obj.beta(j) > 0 || obj.classes(j) == -1;
            else
                error('Object has no beta/classes vector');
            end
        end

        function b = is_positive(obj, j)
            if obj.has_low && obj.has_classes
                b = obj.low(j) == 0 || obj.classes(j) == 1;
            else
                error('Object has no beta/classes vector');
            end
        end

        % ===================
        % Get/Find functions
        % ===================

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

        function idx = find_arglow(obj)
            idx = obj.arglow(obj.arglow ~= 0);
        end

        function d_simplices = get_d_simplices(obj, d)
            d_simplices = find(obj.initial_dimensions == d);
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

        % ===================
        % Reduce
        % ===================

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

        function col = left_to_right(obj, l, j)
            obj.record_column_add(l, j);
            obj.matrix(:, j) = mod(obj.matrix(:, j) + obj.matrix(:, l), 2);
            obj.low(j) = obj.get_low(j);
        end

        function reduce_col(obj, j)
            while (obj.low(j) > 0 && (obj.arglow(obj.low(j)) ~=0))
                j0 = obj.arglow(obj.low(j));
                obj.record_column_add(j0, j);
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
            neg_ind = (obj.alpha == obj.beta) & (obj.beta > 0);
            neg = find(neg_ind);
            obj.mark_as_negative(neg);

            % Clear and mark associated positive columns
            pos_pairs = obj.beta(neg);
            obj.clear_cols(pos_pairs);

            % Mark arglows
            obj.arglow(pos_pairs) = neg;

            % Record lowstars
            obj.lowstar(pos_pairs) = 0;
            obj.lowstar(neg) = pos_pairs;
        end

        % ===================
        % Persistence information vectors
        % ===================

        % -------------------
        % Create
        % -------------------

        function init(obj)
            % Wraps the creation of
            % low, lowstar, arglow, classes
            obj.create_low();
            obj.create_lowstar();
            obj.create_arglow();
            obj.create_classes();
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

        function create_lowstar(obj)
            if ~obj.has_lowstar
                % Create lowstar 
                obj.lowstar = -1*ones(1, obj.m);
                % Mark lowstar as created
                obj.has_lowstar = true;
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

        % -------------------
        % Mark
        % -------------------

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

        function set_unmarked(obj)
            % Mark negatives
            neg_ind = obj.classes == 0 & obj.low > 0;
            obj.arglow(neg_ind) = find(neg_ind);
            obj.mark_as_negative(neg_ind);
            % Mark positives
            pos = obj.low(neg_ind);
            obj.clear_cols(pos);
            % Mark essentials
            idx = obj.classes == 0;
            obj.mark_as_essential(idx);
        end

        % ===================
        % Auxiliary debugging functions
        % ===================

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

        function print_ph_status(obj)
            num_positives = nnz(obj.classes == 1);
            num_negatives = nnz(obj.classes == -1);
            num_essential = nnz(obj.classes == Inf);
            num_unknown = nnz(obj.classes == 0);
            fprintf('\nStats:\n');
            fprintf('\tpositives = %d\n', num_positives);
            fprintf('\tnegatives = %d\n', num_negatives);
            fprintf('\tessential = %d\n', num_essential);
            fprintf('\tunknown = %d\n', num_unknown);
            fprintf('\ttotal = %d\n', obj.m);
        end

        % ===================
        % Visualisation masks
        % ===================

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

    end
end


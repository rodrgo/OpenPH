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

        %% 
        % Query information
        %

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

        %% 
        % Modify matrix
        %
        
        function obj = as_dense(obj)
            obj.matrix = full(obj.matrix);
        end

        function obj = clear_cols(obj, idx)
            obj.matrix(:, idx) = 0;
            obj.low(idx) = 0;
            obj.mark_as_positive(idx);
        end

        function obj = reduce_col(obj, j)
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

        function obj = alpha_beta_reduce(obj)
            % Create relevant vectors
            obj = obj.create_alpha();
            obj = obj.create_beta();

            % Pick obvious lows from alpha/beta pairs
            neg = (obj.alpha == obj.beta) & (obj.beta > 0);
            obj.mark_as_negative(neg);

            % Clear and mark associated positive columns
            pos_pairs = obj.beta(neg);
            obj.clear_cols(pos_pairs);

            % Mark arglows
            obj.arglow(pos_pairs) = find(neg);
        end

        function obj = rho_clearing(obj)

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

        %% 
        % Persistence information vectors
        %

        % This function wraps the creation of
        % low, arglow, classes
        function obj = init(obj)
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

        function obj = create_classes(obj)
            if ~obj.has_classes
                obj.classes = zeros(1, obj.m);
                obj.has_classes = true;
            end
        end

        function obj = create_low(obj)
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

        function obj = create_arglow(obj)
            if ~obj.has_arglow
                % Create arglow
                obj.arglow = zeros(obj.m, 1);
                % Mark arglow as created
                obj.has_arglow = true;
            end
        end


        function obj = create_left(obj)
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

        function obj = create_alpha(obj)
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

        function obj = create_rho(obj)
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

        function obj = create_beta(obj)
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

        %% 
        % Persistence markers
        %

        function obj = mark_classes(obj)
            obj.mark_as_positive(obj.low == 0);
            obj.mark_as_negative(obj.low > 0);
            % Essentials are positive columns that are not
            % paired with any negative column.
            obj.mark_as_essential((obj.classes == 1)' & obj.arglow == 0);
        end

        % Positives can either be paired or essential
        function obj = mark_as_positive(obj, idx)
            obj.classes(idx) = +1;
        end

        function obj = mark_as_negative(obj, idx)
            obj.classes(idx) = -1;
        end

        function obj = mark_as_essential(obj, idx)
            obj.classes(idx) = Inf;
        end

        %% 
        % Visualisation masks
        %

        function mask = get_alpha_mask(obj)
            idx = find(obj.alpha);
            mask = sparse(obj.alpha(idx), idx, ones(size(idx)), obj.m, obj.m);
        end

        function mask = get_beta_mask(obj)
            idx = find(obj.beta);
            mask = sparse(obj.beta(idx), idx, ones(size(idx)), obj.m, obj.m);
        end

        function mask = get_left_mask(obj)
            idx = find(obj.left);
            mask = sparse(idx, obj.left(idx), ones(size(idx)), obj.m, obj.m);
        end

    end
end

%        function [lows, t] = reduce(obj, algorithm)
%
%            if strcmp(algorithm, 'std_sparse')
%                [lows, t] = obj.standard_reduction_sparse();
%            elseif strcmp(algorithm, 'std_sparse_opt') 
%                [lows, t] = obj.standard_reduction_sparse_opt();
%            elseif strcmp(algorithm, 'std_dense_opt') 
%                [lows, t] = obj.standard_reduction_dense_opt();
%            elseif strcmp(algorithm, 'twist_sparse') 
%                [lows, t] = obj.twist_reduction_sparse();
%            elseif strcmp(algorithm, 'twist_dense') 
%                [lows, t] = obj.twist_reduction_dense();
%            elseif strcmp(algorithm, 'rho_reduction') 
%                [lows, t] = obj.rho_reduction();
%            else
%                fprintf('Algorithm not identified\n');
%                assert(1 == 0);
%            end
%
%        end
%
%        function [lows, t] = twist_reduction_dense(obj)
%            tic;
%            R = full(obj.matrix);
%            obj = obj.create_low();
%            lows = obj.low;
%            complex_dim = obj.complex_dimension;
%            simplex_dim = obj.dimensions;
%            L = zeros(obj.m, 1);
%            for d = complex_dim:-1:1
%                for j = 1:obj.m
%                    if sum(simplex_dim(j)) == d
%						% L(lows(j)) == 0 iff lows(j) = low_star(j)
%                        while (lows(j) > 0) && (L(lows(j)) ~= 0)
%                            j0 = L(lows(j));
%                            R(:, j) = mod(R(:,j) + R(:,j0), 2);
%                            % Update lows
%                            if any(R(:, j))
%                                lows(j) = find(R(:,j), 1, 'last');
%                            else
%                                lows(j) = 0;
%                            end
%                        end
%                        if lows(j) > 0
%                            i = lows(j);
%                            L(i) = j;
%                            R(:, i) = 0;
%                        end
%                    end
%                end
%            end
%            t = toc;
%        end
%
%        function [lows, t] = twist_reduction_sparse(obj)
%            tic;
%            R = obj.matrix;
%            obj = obj.create_low();
%            lows = obj.low;
%            complex_dim = obj.complex_dimension;
%            simplex_dim = obj.dimensions;
%            L = zeros(obj.m, 1);
%            for d = complex_dim:-1:1
%                for j = 1:obj.m
%                    if sum(simplex_dim(j)) == d
%						% L(lows(j)) == 0 iff lows(j) = low_star(j)
%                        while (lows(j) > 0) && (L(lows(j)) ~= 0)
%                            j0 = L(lows(j));
%                            R(:, j) = mod(R(:,j) + R(:,j0), 2);
%                            % Update lows
%                            if any(R(:, j))
%                                lows(j) = find(R(:,j), 1, 'last');
%                            else
%                                lows(j) = 0;
%                            end
%                        end
%                        if lows(j) > 0
%                            i = lows(j);
%                            L(i) = j;
%                            R(:, i) = 0;
%                        end
%                    end
%                end
%            end
%            t = toc;
%        end
%
%        function [lows, t] = standard_reduction_sparse(obj)
%            tic;
%            R = obj.matrix;
%            obj = obj.create_low();
%            lows = obj.low;
%            for j = 1:obj.m
%                j0 = find(lows(1:(j-1)) ~= 0 & lows(1:(j-1)) == lows(j));
%                while ~isempty(j0)
%                    R(:, j) = mod(R(:, j) + R(:, j0), 2);
%                    if any(R(:, j))
%                        lows(j) = find(R(:, j), 1, 'last');
%                    else
%                        lows(j) = 0;
%                    end
%                    j0 = find(lows(1:(j-1)) ~= 0 & lows(1:(j-1)) == lows(j));
%                end
%            end
%            t = toc;
%        end
%
%        function [lows, t] = standard_reduction_sparse_opt(obj)
%            tic;
%            R = obj.matrix;
%            obj = obj.create_low();
%            lows = obj.low;
%            L = zeros(size(lows)); % lowstar(i) == L(i)
%            for j = 1:obj.m
%                R_j = R(:, j);
%                while (lows(j) > 0 && (L(lows(j)) ~= 0))
%                    j0 = L(lows(j));
%					R_j = mod(R_j + R(:, j0), 2);
%					if any(R_j)
%						lows(j) = find(R_j, 1, 'last');
%					else
%						lows(j) = 0;
%					end
%                end
%                R(:, j) = R_j;
%                if lows(j) > 0
%                    L(lows(j)) = j;
%                end
%            end
%            t = toc;
%        end
%
%        function [lows, t] = standard_reduction_dense_opt(obj)
%            tic;
%            R = full(obj.matrix);
%            obj = obj.create_low();
%            lows = obj.low;
%            L = zeros(size(lows));
%            for j = 1:obj.m
%                while (lows(j) > 0 && (L(lows(j)) ~= 0))
%                    j0 = L(lows(j));
%					R(:, j) = mod(R(:, j) + R(:, j0), 2);
%					if any(R(:, j))
%						lows(j) = find(R(:, j), 1, 'last');
%					else
%						lows(j) = 0;
%					end
%                end
%                if lows(j) > 0
%                    L(lows(j)) = j;
%                end
%            end
%            t = toc;
%        end
%        function [lows, t] = rho_reduction(obj)
%            tic;
%            R = obj.matrix;
%            obj = obj.create_left();
%            obj = obj.create_low();
%            obj = obj.create_alpha();
%            obj = obj.create_beta();
%            obj = obj.create_rho();
%            lows = obj.low;
%            % Look for real alpha/beta pairs
%            % True low stars
%            reduced = zeros(1, obj.m);
%            immediate_lows = obj.alpha == obj.beta';
%            reduced(immediate_lows) = 1;
%            zero_rho = find((obj.rho' == 0) & (1:obj.m > 1));
%            for j = zero_rho
%                    idx = (obj.beta' == 0) & (1:obj.m < j);
%                    lows(idx) = 0;
%                    reduced(idx) = 1;
%            end
%            not_reduced = find(reduced == 0);
%            % Reduce columns
%            L = lows; % lowstar(i) == L(i)
%            %idx = find(lows > 0);
%            %L(lows(idx)) = idx;
%            for j = 1:obj.m
%                R_j = R(:, j);
%                while (lows(j) > 0 && (L(lows(j)) ~= 0))
%                    j0 = L(lows(j));
%					R_j = mod(R_j + R(:, j0), 2);
%					if any(R_j)
%						lows(j) = find(R_j, 1, 'last');
%					else
%						lows(j) = 0;
%					end
%                end
%                R(:, j) = R_j;
%                if lows(j) > 0
%                    L(lows(j)) = j;
%                end
%            end
%            t = toc;
%        end

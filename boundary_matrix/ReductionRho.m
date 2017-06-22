classdef ReductionRho < BoundaryMatrix & handle
    properties

        rho;
        has_rho;

    end

    methods
        function obj = ReductionRho(stream, homology_mode)
            obj = obj@BoundaryMatrix(stream, homology_mode);
            obj.has_rho = false;
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

    end
end


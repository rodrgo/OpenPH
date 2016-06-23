classdef BoundaryMatrixReduction < BoundaryMatrix & handle
    methods
        function obj = BoundaryMatrixReduction(stream, mode)
            obj@BoundaryMatrix(stream, mode);
        end

        %%
        % reducer
        %
        function obj = reduce(obj, algorithm)
            if strcmp(algorithm, 'std_sparse')
                obj.std_red();
            elseif strcmp(algorithm, 'std_dense') 
                obj.as_dense();
                obj.std_red();
            elseif strcmp(algorithm, 'twist_sparse') 
                obj.twist_red();
            elseif strcmp(algorithm, 'twist_dense') 
                obj.as_dense();
                obj.twist_red();
            elseif strcmp(algorithm, 'rho_sparse') 
                obj.rho_red();
            else
                assert(false, 'Algorithm does not exist');
            end
        end

        %% 
        % std_red algorithm
        %
        function obj = std_red(obj)
            obj.create_low();
            obj.create_arglow();
            for j = 1:obj.m
                obj.reduce_col(j);
            end
        end

        %% 
        % twist_red algorithm
        %
        function obj = twist_red(obj)
            obj.create_low();
            complex_dim = obj.complex_dimension;
            for d = complex_dim:-1:1
                d_simplices = obj.get_simplices(d);
                for j = d_simplices
                    obj.reduce_col(j);
                    if obj.low(j) > 0
                        i = obj.low(j);
                        obj.clear_col(i);
                    end
                end
            end
        end

        %% 
        % rho_red algorithm
        %
        function obj = rho_red(obj)
            % Mark as lowstar all columns with alpha=beta
            obj.pick_obvious_lows();
            % Try clearning more via rho curve
            obj.create_rho();
            obj.rho_clearing();
            % Reduce remaining columns
            unreduced_cols = obj.get_unreduced_cols();
            for j = unreduced_cols
                obj.reduce_col(j);
            end
        end
    end
end

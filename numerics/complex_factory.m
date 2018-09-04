% COMPLEX_FACTORY.M
% example_factory wrapper

function [stream, str_cell] = complex_factory(complex_name, PARAMS)
    cn  = complex_name;
    md  = PARAMS.max_dim;
    mfv = PARAMS.max_filtr_val; 
    nd  = PARAMS.num_divs;
    np  = PARAMS.num_points;
    [stream, str_cell] = example_factory(cn, md, mfv, nd, np); 
end


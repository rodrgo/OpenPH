% name_change
%
function new_name = name_change(name)

    if strcmp(name, 'alpha_beta_parallel')
        new_name = 'pms';
    else
        new_name = name;
    end

end

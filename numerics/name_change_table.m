% name_change
%
function new_name = name_change_table(name)

    if strcmp(name, 'random_gaussian')
        new_name = 'Gaussian';
    elseif strcmp(name, 'random_figure_8')
        new_name = 'Figure-8';
    elseif strcmp(name, 'random_trefoil_knot')
        new_name = 'Trefoil-knot';
    elseif strcmp(name, 'sphere_product')
        new_name = 'Sphere-product';
    else
        new_name = name;
    end

end

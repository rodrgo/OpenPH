function [rows, cols, m] = pointcloud_to_boundary_matrix(filename, max_dimension, num_steps, max_filtration_value)

    import edu.stanford.math.plex4.*;

    % Load the input file
    point_cloud = load(filename);

    % If max_filtration_value is not provided, compute the
    % maximum distance between any two points
    if nargin < 4
        max_filtration_value = -1;
        for i = 1:size(point_cloud, 1)
            for ii = (i+1):size(point_cloud, 1) 
                dist = norm(point_cloud(i, :) - point_cloud(ii, :), 2);
                if dist > max_filtration_value
                    max_filtration_value = dist;
                end
            end
        end
    end

    stream = api.Plex4.createVietorisRipsStream(point_cloud, ...
        max_dimension,...
        max_filtration_value, ...
        num_steps);

    [rows, cols, m] = stream2cmo(stream);

end

function stream = pointcloud_to_vr_stream(filename, max_dimension, max_filtration_value, num_steps)
    load_javaplex;
    point_cloud = load(filename);
    filename = regexprep(filename,'.txt',''); 

    stream = api.Plex4.createVietorisRipsStream(point_cloud, ...
        max_dimension,...
        max_filtration_value, ...
        number_steps);
end


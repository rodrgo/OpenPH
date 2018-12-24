
init;
exclude = { 'Vicsek__particles_300_distance_1_noise_0.1_v0_0.03_box_5_timestep_150_of_300.txt', ...
            'Vicsek__particles_300_distance_1_noise_0.1_v0_0.03_box_5_timestep_300_of_300.txt', ...
            'Vicsek__particles_300_distance_1_noise_2_v0_0.03_box_7_timestep_300_of_600.txt', ...
            'Vicsek__particles_300_distance_1_noise_2_v0_0.03_box_7_timestep_600_of_600.txt', ...
            'celegans_weighted_undirected_reindexed_for_matlab.txt_maxdist_2.6429_SP_distmat.txt_point_cloud.txt'};
listing = dir('../datasets/pointclouds');

% Output here
% Only compute tests on complexes where Javaplex can create a stream 
fid = fopen('pointclouds_stream_success.log','w');

for i = 1:length(listing)
    if endsWith(listing(i).name, '.txt')
        fpath = [listing(i).folder '/' listing(i).name];
        display(listing(i).name);
        d = 5; % max_dimension 
        n = 5; % num_steps
        p = 3; % max_filtration_value
        try
            [r, c, m] = pointcloud_to_boundary_matrix(fpath, d, n, p);
            fprintf(fid, [listing(i).name '\n']);
            display(m);
        catch e
            display('Error');
            %e.message
            if(isa(e,'java.lang.OutOfMemoryError'))
                ex = e.ExceptionObject;
                assert(isjava(ex));
                ex.printStackTrace;
            end
        end
    end
end

fclose(fid);

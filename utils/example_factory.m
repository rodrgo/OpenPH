% EXAMPLE_FACTORY.M
% example_number = 1 creates Morozov cubic time example
% example_number > 1 creates Vietoris-Rips complexes from a point cloud
% sampled from a particular geometric object.
%
% api.Plex4.createVietorisRipsStream
% Creates Vietoris-Rips filtered simplicial complex.
% A simplex [v_0, ..., v_k] is in VR(r) if and only if the pairwise distances
% d(v_i, v_j) are less than r for all 0 <= i, j <= k. Thus a Vietoris-Rips
%
% max_dimension:        maximum dimension of the complex
% max_filtration_value: (max_distance) The maximum allowable distance (r).
% num_divisions:        Number of divisions in [0, max_filtration_value).
%
% e.g. (max_dimension, max_filtration_value, num_divisions) = (3, 4, 20)

function [stream, str_cell] = example_factory(example_name, max_dimension, max_filtration_value, num_divisions, num_points)

import edu.stanford.math.plex4.*;

is_vietoris_rips = 1 == 1;

if strcmp(example_name, 'morozov')
    % Morozov Cubic Time Example
    is_vietoris_rips = 1 == 0;
    order_complexity = 7;

    name_str = 'Morozov';
    geo_params_str = sprintf('order=%d', order_complexity);

    stream = examples.MorozovCubicTimeExample.getMorozovCubicTimeExample(order_complexity);
elseif strcmp(example_name, 'house')
    % House Example
    point_cloud = examples.PointCloudExamples.getHouseExample();

    name_str = 'House';
    geo_params_str = '';

    stream = api.Plex4.createVietorisRipsStream(point_cloud, max_dimension, max_filtration_value, num_divisions);

elseif strcmp(example_name, 'icosahedron')
    % Icosahedron Vertices
    point_cloud = examples.PointCloudExamples.getIcosahedronVertices();

    name_str = 'Icosahedron';
    geo_params_str = '';

    stream = api.Plex4.createVietorisRipsStream(point_cloud, max_dimension, max_filtration_value, num_divisions);

elseif strcmp(example_name, 'random_figure_8')
    % Random Figure-8 Points
    if nargin < 5
        num_points = 10;
    end
    point_cloud = examples.PointCloudExamples.getRandomFigure8Points(num_points);

    name_str = 'Random-Figure-8';
    geo_params_str = sprintf('|S|=%d', num_points);

    stream = api.Plex4.createVietorisRipsStream(point_cloud, max_dimension, max_filtration_value, num_divisions);

elseif strcmp(example_name, 'random_torus')
    % Random Torus Points
    if nargin < 5
        num_points = 10;
    end
    inner_radius = 2;
    outer_radius = 5;
    point_cloud = examples.PointCloudExamples.getRandomTorusPoints(num_points, inner_radius, outer_radius);

    name_str = 'Random-Torus';
    geo_params_str = sprintf('|S|=%d, r=%d, R=%d', num_points, inner_radius, outer_radius);

    stream = api.Plex4.createVietorisRipsStream(point_cloud, max_dimension, max_filtration_value, num_divisions);

elseif strcmp(example_name, 'sphere_product')
    % Random Sphere Product Points
    if nargin < 5
        num_points = 10;
    end
    sphere_dimension = 1;
    num_factors = 2;
    point_cloud = examples.PointCloudExamples.getRandomSphereProductPoints(num_points, sphere_dimension, num_factors);

    name_str = 'Random-Sphere-Product';
    geo_params_str = sprintf('|S|=%d, dim=%d, factors=%d', num_points, sphere_dimension, num_factors);

    stream = api.Plex4.createVietorisRipsStream(point_cloud, max_dimension, max_filtration_value, num_divisions);

elseif strcmp(example_name, 'random_trefoil_knot')
    % Random Trefoil Knot Points
    if nargin < 5
        num_points = 10;
    end
    point_cloud = examples.PointCloudExamples.getRandomTrefoilKnotPoints(num_points);

    name_str = 'Random-Trefoil-Knot';
    geo_params_str = sprintf('|S|=%d', num_points);

    stream = api.Plex4.createVietorisRipsStream(point_cloud, max_dimension, max_filtration_value, num_divisions);

elseif strcmp(example_name, 'random_gaussian')
    % Gaussian Points
    if nargin < 5
        num_points = 10;
    end
    dimension = 3;
    point_cloud = examples.PointCloudExamples.getGaussianPoints(num_points, dimension);

    name_str = 'Random-Gaussian';
    geo_params_str = sprintf('|S|=%d, dim=%d', num_points, dimension);    

    stream = api.Plex4.createVietorisRipsStream(point_cloud, max_dimension, max_filtration_value, num_divisions);
end

example_str = [name_str];
params_str = geo_params_str;
filtration_str = sprintf('m=%d', stream.getSize());

str_cell = {example_str, params_str, filtration_str};

end

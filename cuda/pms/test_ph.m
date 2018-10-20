
% =============
% Init
% =============

clc; clear all; close all;

% Set javaplex directory
JAVAPLEX_DIR = '/home/rodrigo/workspace/persistent_homology/javaplex';
MATLAB_JAVAPLEX_DIR = fullfile(JAVAPLEX_DIR, 'src/matlab/for_distribution');

javaaddpath([MATLAB_JAVAPLEX_DIR '/lib/javaplex.jar']);
import edu.stanford.math.plex4.*;

javaaddpath([MATLAB_JAVAPLEX_DIR '/lib/plex-viewer.jar']);
import edu.stanford.math.plex_viewer.*;

addpath([MATLAB_JAVAPLEX_DIR '/utility']);

% =============
% Import Javaplex
% =============

import edu.stanford.math.plex4.*;

% =============
% Get matrix from example
% =============

% House Example
point_cloud = examples.PointCloudExamples.getHouseExample();
stream = api.Plex4.createVietorisRipsStream(point_cloud, 2, 5, 5);

ccs  = streams.utility.StreamUtility.getCompressedBoundaryMatrix(stream);
m 	 = stream.getSize();
rows = double(cell2mat(cell(ccs.get(0).toArray())));
cols = double(cell2mat(cell(ccs.get(1).toArray())));
vals = ones(size(rows));

% Create matrix
R = sparse(rows, cols, vals); 
[r, c, v] = find(R);

col_width = 7;

% Solve with standard
[low_star, t] = std_test(R);

algorithms = {'standard_parallel', 'twist_parallel', 'ph_row_parallel', 'pms', ...
              'standard', 'twist'};

for i = 1:length(algorithms)
    alg = algorithms{i};
    % Test CUDA solvers
    tic;
    [low o2 o3 o4 o5 o6 o7 o8] = ph(alg, int32(r), int32(c), int32(m), int32(col_width), int32(low_star));
    display(sprintf('%s: OK  (%s seconds)', alg, toc));
    alg_ok = isequal(low_star, low);
    if ~alg_ok
        disp(sprintf('Mismatch in %s', alg));
    end
end


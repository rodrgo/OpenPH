
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
[rows, cols, vals] = find(R);

% Solve with standard
[low_test, t] = std_test(R);

% Test CUDA solvers
tic;
[low resRecord timeRecord] = ph('std', int32(rows), int32(cols), int32(vals), int32(m));
display(sprintf('STANDARD: The tests completed after %s seconds.',toc));
std_ok = isequal(low_test-1, low);
if ~isequal(low_test-1, low)
    display(low_test-1);
    display(low);
end

tic;
[low resRecord timeRecord] = ph('twist', int32(rows), int32(cols), int32(vals), int32(m));
display(sprintf('TWIST: The tests completed after %s seconds.',toc));
twist_ok = isequal(low_test-1, low);
if ~isequal(low_test-1, low)
    display(low_test-1);
    display(low);
end


tic;
[low resRecord timeRecord] = ph('pms', int32(rows), int32(cols), int32(vals), int32(m));
display(sprintf('PMS: The tests completed after %s seconds.',toc));
pms_ok = isequal(low_test-1, low);
if ~isequal(low_test-1, low)
    display(low_test-1);
    display(low);
end

if std_ok && twist_ok && pms_ok
    disp(sprintf('All good :)'));
else
    disp(sprintf('Mismatch in one of the algos :)'));
end


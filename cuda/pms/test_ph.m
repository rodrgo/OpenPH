
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

R = sparse(rows, cols, vals); 
[rows, cols, vals] = find(R);
rowsd = rows - 1;
colsd = cols - 1;

%disp(rowsd')
%disp(colsd')

rows = int32(rows);
cols = int32(cols);
vals = int32(vals);
m = int32(m);

%disp(rows)
%disp(cols)
%disp(vals)
%disp(m)

% vals is equal to all ones
tic;
[low resRecord timeRecord] = ph('std', rows, cols, vals, m);
%disp(low)
display(sprintf('The tests completed after %s seconds.',toc));

[low_test, t] = std_test(R);

disp(low);
disp(low_test-1);


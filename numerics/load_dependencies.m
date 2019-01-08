% This script loads javaplex components

clc; clear all; close all;

[dependency_dir, ~, ~] = fileparts(mfilename('fullpath'));
JAVAPLEX_DIR = fullfile(dependency_dir, 'dependencies', 'javaplex');
PHAT_DIR = fullfile(dependency_dir, 'dependencies', 'phat');

% -----------
% Javaplex
% -----------

MATLAB_JAVAPLEX_DIR = fullfile(JAVAPLEX_DIR, 'src/matlab/for_distribution');

javaaddpath([MATLAB_JAVAPLEX_DIR '/lib/javaplex.jar']);
import edu.stanford.math.plex4.*;

javaaddpath([MATLAB_JAVAPLEX_DIR '/lib/plex-viewer.jar']);
import edu.stanford.math.plex_viewer.*;

addpath([MATLAB_JAVAPLEX_DIR '/utility']);

% -----------
% Phat
% -----------

addpath(PHAT_DIR);


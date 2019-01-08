
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

% -----------
% Other paths
% -----------

addpath('./utils/');
addpath('../src/');
addpath('../src/cuda/pms');

% -----------
% Plot init
% -----------

FIGURE_DIR  = './figures/';

% Font size
fs          = [];
fs.title    = 20;
fs.legend   = 17;
fs.axis     = 20;
fs.ticks    = 20;


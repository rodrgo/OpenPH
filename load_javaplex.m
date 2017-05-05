% This script prepares the javaplex library for use

clc; clear all; close all;
%clear import;

MATLAB_JAVAPLEX_DIR = ['/home/mendozasmith/src/javaplex/src/matlab/for_distribution'];

javaaddpath([MATLAB_JAVAPLEX_DIR '/lib/javaplex.jar']);
import edu.stanford.math.plex4.*;

javaaddpath([MATLAB_JAVAPLEX_DIR '/lib/plex-viewer.jar']);
import edu.stanford.math.plex_viewer.*;

addpath([MATLAB_JAVAPLEX_DIR '/utility']);


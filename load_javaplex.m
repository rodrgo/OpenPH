% This script prepares the javaplex library for use

clc; clear all; close all;
clear import;

matlab_javaplex_dir = '../javaplex/src/matlab/for_distribution';

javaaddpath([matlab_javaplex_dir '/lib/javaplex.jar']);
import edu.stanford.math.plex4.*;

javaaddpath([matlab_javaplex_dir '/lib/plex-viewer.jar']);
import edu.stanford.math.plex_viewer.*;

addpath([matlab_javaplex_dir '/utility']);


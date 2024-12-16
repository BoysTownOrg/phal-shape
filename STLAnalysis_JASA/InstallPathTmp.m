% Run script to temporarily reset the MATLAB search path and add the
% required folders for the current MATLAB session only.
restoredefaultpath
addpath(genpath([pwd '\altmany-export_fig-3.47.0.0']))
addpath(genpath([pwd '\extrema']))
addpath(genpath([pwd '\gridtrimesh']))
addpath(genpath([pwd '\polygon2voxel']))
addpath(genpath([pwd '\utilities']))
clear all
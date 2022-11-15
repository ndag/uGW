% Determine where this file's folder is.
folder = fileparts(which('setup.m')); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

mex uGW-distance/construct_cost_one.cpp
mex uGW-distance/construct_ultra_cost_mat.cpp
mex dGW-distance/construct_cost_mat_dGW.cpp
mex dGW-distance/construct_cost_one_dGW.cpp

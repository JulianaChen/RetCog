%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% MATLAB CODE RETIREMENT AND COGNITION %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

%% Set up Parameters
run setup_params.m

%% Set up State Space
S = sspace(eparams,G); % optimal assets for poly, sspace
%S = sspace_lin(eparams,G); % linear assets, sspace3

%% Solution (loop all 18 types)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% MATLAB CODE RETIREMENT AND COGNITION %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

%% Set up Parameters
run setup_params.m

%% Set up State Space
S = sspace(G); % optimal assets for poly, sspace
%S = sspace(eparams,G); % with eparams

%% Solution (loop all 18 types)
function [numeric_version]=return_dynare_version(dynare_version)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% return_dynare_version.m
%
% This file constructs a numeric version of dynare from the string provided
% by Dynare.
% 
%Included in the package nonlinear_MA
%
%THIS VERSION: 1.0.8 Jun. 3, 2013
%
%Copyright: Hong Lan and Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%This software is provided as is with no guarantees of any kind
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynare_ver=dynare_version;
[dynare_ver_1,dynare_ver_remain] = strtok(dynare_ver, '.');
[dynare_ver_2,dynare_ver_remain2] = strtok(dynare_ver_remain, '.');
[dynare_ver_3,dynare_ver_remain3] = strtok(dynare_ver_remain2, '.');
numeric_version=str2num([dynare_ver_1 '.' dynare_ver_2 dynare_ver_3]);
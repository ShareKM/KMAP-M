% This is a demo file to test run_ParametricImage with a real total-body
% data
clear; clc;
%% setup function path
% PLOT_v1.0 package can be downloaded from 
% https://wanglab.faculty.ucdavis.edu/code
run('E:\Github_repo\PLOT_v1.0\setup');

addpath('../data/');
run('../setup.m');
% the image data can be downloaded from
% https://drive.google.com/drive/folders/1OWG0FMlhuZL3vnv_fswXSaqXo8_762_5?usp=drive_link
dynImgInp = 'HS01_test_data.mat';
foldOut = './output_paramImg';
fileBIF = 'Input_func_AA.mat';

Dopt.frameType  = '1H29';
Dopt.kernelFlag = 1;
Dopt.smoothFlag = 1;
Dopt.delayFlag  = 1;
if ~exist(foldOut, 'dir')
    mkdir(foldOut)
end
run_ParametricImage(dynImgInp, foldOut, fileBIF, Dopt)
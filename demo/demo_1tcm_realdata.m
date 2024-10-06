% This is a demo file to test 1TCM model with real TAC data
clear; clc

run('../setup.m');
addpath('../data/');
% load data
load('demo_1tcm_realdata.mat');
tac = WM;
scant = [cumsum([0; dt(1:end-1)]), cumsum(dt)];
t = mean(scant,2);

% optional setting
opt.ScanTime = scant;
opt.Decay    = 0;
opt.TimeStep = 1.0;
opt.PbrParam = [1 0 0];

opt.MaxIter    = 100;
opt.LowerBound = [0 0 0 0];
opt.UpperBound = [1 2 2 20];
opt.PrmSens    = [1 1 1 1];
kinit          = [0.1 0.1 0.1 0]';
wt = ones(size(tac));

% fit the TAC
[kfit, cfit] = kfit_1tcm(tac, Cp, scant, kinit, opt, wt);

figure
plot(t,tac, '-.o', t, cfit)
legend('Noisy TAC','Fitted TAC','Location','Best')

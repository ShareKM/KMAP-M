% This is a demo file to test 2TCM model with real TAC data
clear; clc

run('../../setup.m');
addpath('./data/');
% load data
load('demo_2tcm_realdata.mat');
tac = GM;
scant = [cumsum([0; dt(1:end-1)]), cumsum(dt)];
t = mean(scant,2);

% optional setting
opt.ScanTime = scant;
opt.Decay    = log(2)/109.8;
opt.TimeStep = 1.0;
opt.PbrParam = [1 0 0];

opt.MaxIter    = 100;
opt.LowerBound = [0 0 0 0 0 0];
opt.UpperBound = [1 2 2 2 0.5 20];
opt.PrmSens    = [1 1 1 1 1 1];
kinit          = [0.1 0.1 0.1 0.1 0.01 0]';
wt = ones(size(tac));

% fit the TAC
[kfit, cfit] = kfit_2tcm(tac, Cp, scant, kinit, opt, wt);

figure
plot(t,tac, '-.o', t, cfit)
legend('Noisy TAC','Fitted TAC','Location','Best')

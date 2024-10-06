% This is a demo file to test 2TCM model with simulation data
clear; clc

run('../setup.m');
% scan time setting
dt = [ones(6,1)*10;  ones(2,1)*30;  ones(6,1)*60; ones(5,1)*120; ones(4,1)*180; ones(6,1)*300;];
scant = [cumsum([0; dt(1:end-1)]), cumsum(dt)];
t = mean(scant,2);

% optional setting
opt.ScanTime = scant;
opt.Decay    = 0;
opt.TimeStep = 1.0;
opt.PbrParam = [1 0 0];

% input function
Cp = feng(t);
Cp(1) = 0;
% generating noise-free TAC
k0 = [0.05 0.81 0.38 0.10 0.009 10]'; % = [fv K1 k2 k3 k4 time_delay]
c = ktac_2tcm(k0, Cp, scant, opt);

% test least squares fitting
opt.MaxIter    = 100;
opt.LowerBound = [0 0 0 0 0 0];
opt.UpperBound = [1 2 2 2 0.5 20];
opt.PrmSens    = [1 1 1 1 1 1];
kinit          = [0.1 0.1 0.1 0.1 0.01 0]';
wt = ones(size(c));

% noisy TAC
tac = c + 3*diag(sqrt(c./dt))*randn(size(c));

% fit the noisy TAC
[kfit, cfit] = kfit_2tcm(tac, Cp, scant, kinit, opt, wt);

figure
plot(t,tac, '-.o', t, cfit, t, c)
legend('Noisy TAC','Fitted TAC','True TAC','Location','Best')

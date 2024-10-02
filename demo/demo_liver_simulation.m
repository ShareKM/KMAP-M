% This is a demo file to test 2T Liver model with simulation data
clear; clc

run('../setup.m');
% scan time setting
dt = [ones(30,1)*10;  ones(10,1)*60; ones(9,1)*300;];
scant = [cumsum([0; dt(1:end-1)]), cumsum(dt)];
t = mean(scant,2);

% optional setting
opt.ScanTime = scant;
opt.Decay    = log(2)/109.8;
opt.TimeStep = 1.0;
opt.PbrParam = [1 0 0];

% input function
Cp = feng(t);
Cp(1) = 0;
% generating noise-free TAC
k0 = [0.1, 0.9, 1.1, 0.01, 0.01, 1, 0.05 20]'; % = [fv K1 k2 k3 k4 ka fa time_delay]

c = ktac_liver(k0, Cp, scant, opt);

% test least squares fitting
opt.MaxIter    = 200;
opt.LowerBound = [0, 0, 0, 0, 0, 0, 0, 0];
opt.UpperBound = [1.0, 10, 10.0, 1.0, 0.1, 10, 1, 50];
opt.PrmSens    = [1 1 1 1 1 1 1 1];
kinit          = [0.01, 1, 1, 0.01, 0.01, 1, 0.01, 0]';
wt = ones(size(c));

% noisy TAC
tac = c + 3*diag(sqrt(c./dt))*randn(size(c));
% tac = max(1e-1, tac);

% fit the noisy TAC
[kfit, cfit] = kfit_liver(tac, Cp, scant, kinit, opt, wt);

figure
plot(t,tac, '-.o', t, cfit, t, c)
legend('Noisy TAC','Fitted TAC','True TAC','Location','Best')



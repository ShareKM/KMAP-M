clear; clc
codepath = '../../';
run([codepath,'DIRECT_v1.0/setup.m']);

% scan time setting
dt = [ones(4,1)*20;  ones(4,1)*40;  ones(4,1)*60; ones(4,1)*180; ones(8,1)*300;];
scant = [cumsum([0; dt(1:end-1)]), cumsum(dt)];
t = mean(scant,2);

% optional setting
opt.ScanTime = scant;
opt.Decay    = log(2)/109.8;
opt.TimeStep = 1.0;
opt.PbrParam = [1 0 0];

% input function
Cp = feng(t);

% generating noise-free TAC
k0 = [0.05 0.81 0.38 0.10 0.009]'; % = [fv K1 k2 k3 k4]
c = ktac_2t5p(k0, Cp, scant, opt);
 
% test least squares fitting
opt.MaxIter    = 100;
opt.LowerBound = [0 0 0 0 0];
opt.UpperBound = [1 2 2 2 0.5];
opt.PrmSens    = [1 1 1 1 1];
kinit          = [0.1 0.1 0.1 0.1 0.01]';

% noisy TAC
tac = c + 5*diag(sqrt(60./dt))*randn(size(c));
tac = max(1e-1, tac);

% fit the noisy TAC
wt = ones(size(c));
[kfit, cfit] = kfit_2t5p(tac, Cp, scant, kinit, opt, wt);
figure,plot(t,tac, '-.o', t, cfit(:,1), t, c)
legend('Noisy TAC','Fitted TAC','True TAC','Location','Best')

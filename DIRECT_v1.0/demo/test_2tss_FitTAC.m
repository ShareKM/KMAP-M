clear; clc
codepath = '/home/gbwang/work/';
run([codepath,'DIRECT_v1.0/setup.m']);

%% scan time setting
dt = [ones(40,1)*2;  ones(16,1)*10;  ones(8,1)*30; ones(3,1)*180;];
scant = [cumsum([0; dt(1:end-1)]), cumsum(dt)];
t = mean(scant,2);

% optional setting
opt.ScanTime = scant;
opt.Decay    = log(2)/109.8;
opt.TimeStep = 1.0;
opt.PbrParam = [1 0 0];

% input function
Cp = feng(t);
Cp([1 2]) = [0 10];

% generating noise-free TAC
k0 = [0.1 0.3 0.81 0.68 0.20 0.09]'; % = [vb ve K1 k2 k3 k4]
[c0,s0] = ktac_2tss(k0, Cp, scant, opt);

%% test least squares fitting
opt.MaxIter    = 100;
opt.LowerBound = [0 0 0 0 0 0];
opt.UpperBound = [1 1 2 2 2 0.5];
opt.PrmSens    = [1 0 1 1 1 1];
kinit          = [0.1 k0(2) 0.1 0.1 0.1 0.01]';

% noisy TAC
tac = c0 + 1*diag(sqrt(60./dt))*randn(size(c0));
tac = max(1e-1, tac);

% fit the noisy TAC
wt = ones(size(c0));
[kfit, cfit] = kfit_2tss(tac, Cp, scant, kinit, opt, wt);
figure,plot(t,tac, '-.o', t, cfit(:,1), t, c0)
legend('Noisy TAC','Fitted TAC','True TAC','Location','Best')
[k0 kfit]
Ki0 = mic2mac(k0(3:5))
Ki = mic2mac(kfit(3:5))

    
%% using conventinal model
opt1 = opt;
opt1.LowerBound = [0 0 0 0 0 ];
opt1.UpperBound = [1 2 2 2 0.5];
opt1.PrmSens    = [1 1 1 1 1];
kinit1          = [0.1 0.1 0.1 0.1 0.01]';

% noise-free
[kfit1, cfit1] = kfit_2t5p(c0, Cp, scant, kinit1, opt1, wt);
figure,plot(t,c0, '-.o', t, cfit1(:,1))
legend('Measured TAC','Fitted TAC','Location','Best')
Ki1 = mic2mac(kfit1(2:4))
[k0([1 3:end]) kfit1]

% noisy
[kfit1, cfit1] = kfit_2t5p(tac, Cp, scant, kinit1, opt1, wt);
figure,plot(t,tac, '-.o', t, cfit1(:,1), t, c0)
legend('Noisy TAC','Fitted TAC','True TAC','Location','Best')
Ki1 = mic2mac(kfit1(2:4))
[kfit([1 3:end]) kfit1]
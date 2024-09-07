clear; clc;

% scan time setting
dt = [ones(6,1)*5;  ones(9,1)*10;  ones(6,1)*30;];
dt = [ones(4,1)*20;  ones(4,1)*40;  ones(4,1)*60; ones(4,1)*180; ones(10,1)*300;];
scant = [cumsum([0; dt(1:end-1)]), cumsum(dt)];
t = mean(scant,2);

% parameters
fv = 0.1; k1 = 1.40; k2 = 0.20; 

% optional setting
opt.Decay = log(2)/20.4;
opt.ScanTime = scant;
opt.TimeStep = 1.0;
opt.KinType = '1t3p';
opt.LowerBound = [0 1e-3 1e-3];
opt.UpperBound = [10 2 2];
opt.ParSens = [0 1 1];
kinit = [fv 0.01 0.01]';
opt.MaxIter = 1000;
opt.PbrPar = [1.0 0 0];
opt.NoiseModel = 'gsn';
opt.Algorithm = 'lma';

% input function
cp = gamma3(t);
cwb = cp;

% test for 1-tissue
k0 = [fv k1 k2]';
c = tac_fun(k0, cp, opt);

% noisy TAC
tac = max(1e-2, c + 1000*diag(1./dt)*randn(size(c)) );

% nonlinear fitting
kinit = repmat(kinit(:),[1 size(tac,2)]);
[kfit, cfit] = fit_fun(kinit, cp, opt, tac, []);
figure,plot(t,c,'-.', t, tac(:,1), ':o', t, cfit(:,1))
legend('true','noisy','fitted');

% global solution
kk1 = linspace(opt.LowerBound(2), opt.UpperBound(2), 1000);
kk2 = linspace(opt.LowerBound(3), opt.UpperBound(3), 1000);
L = zeros(length(kk1),length(kk2));
for i = 1:length(kk1)
    kk = repmat([fv, kk1(i), 0]', [1 length(kk2)]);
    kk(3,:) = kk2;
    cc = tac_fun(kk, cp, opt);
    Li = sum(( cc - repmat(tac,[1 length(kk2)]) ).^2, 1);
    L(i,:) = Li;
end
[a, b] = min(L(:));
[idx1, idx2] = ind2sub([length(kk1), length(kk2)], b);
kopt = [kk1(idx1), kk2(idx2)]';
[kfit(2:3), kopt]
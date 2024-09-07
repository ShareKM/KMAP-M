clear; clc
run('/home/gbwang/work/kmap2020/setup.m')

% scan time setting
dt = [ones(4,1)*20;  ones(4,1)*40;  ones(4,1)*60; ones(4,1)*180; ones(8,1)*300;];
scant = [cumsum([0; dt(1:end-1)]), cumsum(dt)];
t = mean(scant,2);

% optional setting
opt.ScanTime = scant;
opt.Decay = log(2)/109.8;
opt.TimeStep = 1.0;
opt.PbrPar = [1 0 0];
opt.KinType = '2tss';

% input function
Cp = feng(t);

% noise-free TAC
k0 = [0.05 0.81 0.38 0.10 0.009]'; % = [fv k1 ki k3 k4]
c = tac_fun(k0, Cp, opt);

% test least squares fitting
opt.MaxIter = 100;
opt.LowerBound = [0 0 0 0 0];
opt.UpperBound = [1 2 2 2 0.5];
opt.ParSens = [1 1 1 1 1];
kinit = [0.1 0.1 0.1 0.1 0.01]';

if 1
    wt = ones(size(c));
    for i = 1:10
        tac(:,i) = c + 1*diag(sqrt(1./dt))*randn(size(c));
    end
    tac = max(1e-1, tac);
    kinit = repmat(kinit(:),[1 size(tac,2)]);
    
    tic;
    [kfit1, cfit1] = fit_fun(kinit, Cp, opt, tac, wt);
    toc;
    figure,plot(t,tac(:,1), '-.o', t, cfit1(:,1), t, c)
    [(mean(kfit1,2)-k0)./k0*100, std(kfit1,0,2)./k0*100]

    tt = 1:18;
    opt2 = opt;
    opt2.ScanTime = scant(tt,:);
    tic;
    [kfit2, cfit2] = fit_fun(kinit, Cp(tt), opt2, tac(tt,:), wt(tt));
    toc;
    figure,plot(t,tac(:,1), '-.o', t(tt), cfit2(tt,1))
    [(mean(kfit2,2)-k0)./k0*100, std(kfit2,0,2)./k0*100]
    
end

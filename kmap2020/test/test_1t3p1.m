clear; clc

% scan time setting
dt = [ones(4,1)*20;  ones(4,1)*40;  ones(4,1)*60; ones(4,1)*180; ones(4,1)*300;];
scant = [cumsum([0; dt(1:end-1)]), cumsum(dt)];
t = mean(scant,2);

% optional setting
opt.ScanTime = scant;
opt.Decay = log(2)/20.4;
opt.TimeStep = 1.0;
opt.PbrPar = [1 0 0];

% input function
Cp = feng(t);

% the reference tac
VDr = 1.0;
vbr = 0.0; k1r = 0.5; k2r = k1r/VDr;
opt.KinType = '1t3p';
cr = tac_fun([vbr k1r k2r]', Cp, opt);

% the target tac
vb = 0.0; k1 = .6; k2 = k1/VDr; BP = 1;
k2a = k2/(1+BP);
ct = tac_fun([vb k1 k2a]', Cp, opt);

% test srtm
R1 = k1/k1r;
k0 = [R1 (BP+1-R1)*k2a k2a]';
newopt = opt;
newopt.KinType = '1t3p1';
newopt.NoiseModel = 'psn';

tic;
Cr = cr.*exp(newopt.Decay*t/60);
[c,s] = tac_fun(k0, Cr, newopt);
toc;
%figure, plot(t,c,t,ct)

% test least squares fitting
newopt.NumOfPar = 3;
newopt.ParSens = [1 1 1];
newopt.LowerBound = [1e-4 1e-4 1e-4];
newopt.UpperBound = [5 5 5];
kinit = [1 .1 .5]';

if 1
    wt = ones(size(c));
    for i = 1:1e3
        tac(:,i) = c + 10*diag(sqrt(c./dt))*randn(size(c));
    end
    tac = max(1e-1, tac);
    kinit = repmat(kinit(:),[1 size(tac,2)]);
    
    newopt.Algorithm = 'bfm';
    newopt.MaxIter = 100;
    tic;
    [kfit1, cfit1] = fit_fun(kinit, Cr, newopt, tac, dt);
    toc;
    figure,plot(t,tac(:,1), '-.o', t, cfit1(:,1), t, c)
    [(mean(kfit1,2)-k0)./k0*100, std(kfit1,0,2)./k0*100]

    newopt.Algorithm = 'bfm2';
    newopt.MaxIter = 10;
    tic;
    [kfit2, cfit2] = fit_fun(kinit, Cr, newopt, tac, dt);
    toc;
    figure,plot(t,tac(:,1), '-.o', t, cfit2(:,1), t, c)
    [(mean(kfit2,2)-k0)./k0*100, std(kfit2,0,2)./k0*100]
    
    newopt.Algorithm = 'lma';
    newopt.MaxIter = 50;
    tic;
    [kfit3, cfit3] = fit_fun(kinit, Cr, newopt, tac, dt);
    toc;
    figure,plot(t,tac(:,1), '-.o', t, cfit3(:,1), t, c)
    [(mean(kfit3,2)-k0)./k0*100, std(kfit3,0,2)./k0*100]
    
end


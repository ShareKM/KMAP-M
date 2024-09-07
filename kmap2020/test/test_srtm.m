clear; clc

% scan time setting
dt = [ones(4,1)*20;  ones(4,1)*40;  ones(4,1)*60; ones(4,1)*180; ones(14,1)*300;];
scant = [cumsum([0; dt(1:end-1)]), cumsum(dt)];
t = mean(scant,2);

% optional setting
opt.ScanTime = scant;
opt.Decay = log(2)/20.4;
opt.TimeStep = 1.0;
opt.PbrPar = [1.25 0 0];

% input function
Cp = feng(t);
%Cp = Cp.*exp(-opt.Decay*t/60); % for consistance with COMKAT

% the reference tac
VDr = 1.0;
vbr = 0.0; k1r = 0.5; k2r = k1r/VDr;
opt.KinType = '1t3p';
cr = tac_fun([vbr k1r k2r]', Cp, opt);

% the target tac
vb = 0.0; k1 = 1; k2 = k1/VDr; BP = 1;
k2a = k2/(1+BP);
ct = tac_fun([vb k1 k2a]', Cp, opt);

% test srtm
R1 = k1/k1r;
k0 = [R1 k2 BP]';
newopt = opt;
newopt.KinType = 'srtm';
newopt.NoiseModel = 'gsn';
newopt.Algorithm = 'lma';

tic;
[c,s] = tac_fun(k0, cr, newopt);
toc;
%figure, plot(t,c,t,ct)

% test sensitivity functions
if 0
    [c1] = tac_fun([R1+1e-5 k2 BP]', cr, newopt);
    s1 = (c1-c)/1e-5; figure,plot([s(:,1),s1]);
    [c1] = tac_fun([R1 k2+1e-5 BP]', cr, newopt);
    s1 = (c1-c)/1e-5; figure,plot([s(:,2),s1]);
    [c1] = tac_fun([R1 k2 BP+1e-5]', cr, newopt);
    s1 = (c1-c)/1e-5; figure,plot([s(:,3),s1]);
end

% test least squares fitting
newopt.MaxIter = 200;
newopt.NumOfPar = 3;
newopt.ParSens = [1 1 1];
newopt.LowerBound = [1e-4 1e-4 1e-4];
newopt.UpperBound = [5 5 5];
kinit = [.01 0.01 .01]'*1e1;

if 1
    wt = ones(size(c));
    for i = 1:100
        tac(:,i) = c + 100*diag(1./dt)*randn(size(c));
    end
    tac = max(0, tac);
    kinit = repmat(kinit(:),[1 size(tac,2)]);
    tic;
    [kfit, cfit] = fit_fun(kinit, cr, newopt, tac);
    toc;
    figure,plot(t,tac(:,1), '-.o', t, cfit(:,1), t, c)
    [(mean(kfit,2)-k0)./k0*100, std(kfit,0,2)./k0*100]
end


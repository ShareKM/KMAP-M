clear; 

% scan time setting
dt = [ones(6,1)*5;  ones(9,1)*10;  ones(6,1)*30;];
dt = [ones(4,1)*20;  ones(4,1)*40;  ones(4,1)*60; ones(4,1)*180; ones(10,1)*300;];
scant = [cumsum([0; dt(1:end-1)]), cumsum(dt)];
t = mean(scant,2);

% parameters
fv = 1.15; k1 = 2.40; k2 = 1.42; 

% optional setting
opt.Decay = log(2)/20.4;
opt.ScanTime = scant;
opt.TimeStep = 1.0;
opt.KinType = '1t3p1';
opt.LowerBound = [0 1e-3 1e-3];
opt.UpperBound = [10 10 10];
opt.ParSens = [1 1 1];
opt.MaxIter = 20;
opt.PbrPar = [1.0 0 0];
opt.NoiseModel = 'gsn';
opt.Algorithm = 'lma';

% input function
cp = gamma3(t);
cwb = cp;

% test for 1-tissue
k0 = [fv k1 k2]';
[c,s] = tac_fun(k0, cp, opt);
%figure, plot(t, c)

if 1
    tic; c1 = tac_fun([fv+1e-5 k1 k2]', cp, opt);toc;
    s1 = (c1-c(:,1))/1e-5; figure,plot([s(:,1),s1]);
    tic; c1 = tac_fun([fv k1+1e-5 k2]', cp, opt);toc;
    s1 = (c1-c(:,1))/1e-5; figure,plot([s(:,2),s1]);
    tic; c1 = tac_fun([fv k1 k2+1e-5]', cp, opt);toc;
    s1 = (c1-c(:,1))/1e-5; figure,plot([s(:,3),s1]);
end

% test least squares fitting
if 1
    wt = ones(size(c));
    for i = 1:10
        tac(:,i) = max(1e-2, c + 100*diag(1./dt)*randn(size(c)) );
    end
    kinit = [0.01 0.01 0.01]'*1e0;
    kinit = repmat(kinit(:),[1 size(tac,2)]);
    tic;
    [kfit, cfit] = fit_fun(kinit, cp, opt, tac, dt);
    toc;
    figure,plot(t,c,'-.', t, tac(:,1), ':o', t, cfit(:,1))
    [k0 mean(kfit,2)]
end
%clear; 

% parameters
fv = 0; k1 = 1; k2 = 0.005;

[scant, cp, ct] = fitanal(k1, k2);
t = mean(scant,2);
dt=scant(:,2)-scant(:,1); dt=dt/60;
% optional setting
opt.Decay = log(2)/1e9;
opt.ScanTime = scant;
opt.TimeStep = 1;
opt.KinType = '1t3p1';
opt.LowerBound = [0 1e-3 1e-3];
opt.UpperBound = [10 10 10];
opt.ParSens = [0 1 1];
opt.MaxIter = 20;
opt.PbrPar = [1.0 0 0];
opt.NoiseModel = 'gsn';
opt.Algorithm = 'lma';

% input function
cwb = cp;

% test for 1-tissue
k0 = [fv k1 k2]';
[c] = tac_fun(k0, cp, opt);
figure, plot(t, ct, 'b:', t, c, 'r');
legend('Analytical truth', 'Numerical calculation');

if 0
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
    for i = 1:1
        tac(:,i) = max(0e-2, ct + 0*diag(1./dt)*randn(size(c)) );
    end
    kinit = [0.000 0.01 0.01]'*1e0;
    kinit = repmat(kinit(:),[1 size(tac,2)]);
    tic;
    [kfit, cfit] = fit_fun(kinit, cp, opt, tac, dt);
    toc;
    figure,plot(t,ct,'-.', t, tac(:,1), ':o', t, cfit(:,1))
    [k0 mean(kfit,2)]
end
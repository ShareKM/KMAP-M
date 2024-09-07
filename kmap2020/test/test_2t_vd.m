clear; clc

% scan time setting
dt = [ones(4,1)*20;  ones(4,1)*40;  ones(4,1)*60; ones(4,1)*180; ones(14,1)*300;];
scant = [cumsum([0; dt(1:end-1)]), cumsum(dt)];
t = mean(scant,2);

% parameters
k0 = [0.05 0.11 0.11/0.19 0.10 0.009]'; % = [fv k1 k1/k2 k3 k4]
k0 = [0.04 0.088 0.088/0.055 0.096 0.003]';
%k0 = [0.03 0.050 0.050/0.082 0.049 0.001]';
%k0 = [0.05 0.079 0.079/0.101 0.059 0.002]';

% optional setting
opt.Decay = log(2)/109.8;
opt.TimeStep = 1.0;
opt.LowerBound = [0.0 0 0 0 0  ];
opt.UpperBound = [0.5 2 5 2 0.5];
opt.MaxIter = 20;
opt.PrmSens = [1 1 1 1 1];
kinit = [0.05 0.01 0.01/0.01 0.01 0.001]';

% input function
cp = feng(t);
cwb = cp;
cp = cp.*exp(-opt.Decay*t/60); % for consistance with COMKAT

% test for 2-tissue
tic; 
[c, s] = ktac_2t5pv(k0, cp, scant, opt, cwb);
toc;
% figure, plot(t, c)

if 0
    for i = 1:length(k0)
        kk = k0; kk(i) = kk(i) + 1e-5;
        tic; [c1]=ktac_2t5pv(kk, cp, scant, opt, cwb);toc;
        s1 = (c1-c(:,1))/1e-5; figure,plot([s(:,i),s1]);
    end
end

% test least squares fitting
if 1
    wt = ones(size(c));
    for i = 1:10000
        tac(:,i) = c + 15*sqrt(c./dt).*randn(size(c));
    end
    tac = max(0, tac);
    
    kinit = repmat(kinit(:),[1 size(tac,2)]);
    tic;
    [kfit, cfit] = kfit_2t5pv(tac, cp, scant, kinit, opt, wt, cwb);
    toc;
    figure,plot(t,tac(:,1), '-.o', t, cfit(:,1))
    [(mean(kfit,2)-k0)./k0*100, std(kfit,0,2)./k0*100]
    
    ki0 = mic2mac2(k0(2:end),'ki');
    ki  = mic2mac2(kfit(2:end,:),'ki');
    [(mean(ki,2)-ki0)./ki0*100, std(ki,0,2)./ki0*100]
    
    
    % the regular model
    opt.UpperBound = [0.5 2 2 2 0.5];
    kinit = [0.05 0.01 0.01 0.01 0.01]';
    kinit = repmat(kinit(:),[1 size(tac,2)]);
    tic;
    [kfit, cfit] = kfit_2t5p(tac, cp, scant, kinit, opt, wt, cwb);
    toc;
    k0(3) = k0(2)/k0(3);
    [(mean(kfit,2)-k0)./k0*100, std(kfit,0,2)./k0*100]
    
    ki  = mic2mac(kfit(2:end,:),'ki');
    [(mean(ki,2)-ki0)./ki0*100, std(ki,0,2)./ki0*100]
    
end
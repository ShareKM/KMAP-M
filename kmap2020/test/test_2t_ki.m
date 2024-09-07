clear; clc;

% scan time setting
dt = [ones(4,1)*20;  ones(4,1)*40;  ones(4,1)*60; ones(4,1)*180; ones(14,1)*300;];
scant = [cumsum([0; dt(1:end-1)]), cumsum(dt)];
t = mean(scant,2);

% parameters
k0 = [0.05 0.11 0.038 0.10 0.009]'; % = [fv k1 ki k3 k4]

% optional setting
opt.Decay = log(2)/109.8;
opt.TimeStep = 1.0;
opt.LowerBound = [0.0 0 0 0 0  ];
opt.UpperBound = [0.3 2 0.1 2 0.5];
opt.MaxIter = 20;
opt.PrmSens = [1 1 1 1 1];
kinit = [0.05 0.01 0.01 0.01 0.01]';

% input function
cp = feng(t);
cwb = cp;
cp = cp.*exp(-opt.Decay*t/60); % for consistance with COMKAT


% test for 1-tissue
tic; 
[c, s] = ktac_2t5pk(k0, cp, scant, opt, cwb);
toc;
%figure, plot(t, c)

% test jac computation
if 0
    for i=1:length(k0)
        kk = k0; kk(i) = kk(i) + 1e-5;
        c1 = ktac_2t5pk(kk, cp, scant, opt, cwb);
        figure, plot([(c1-c)/1e-5, s(:,i)])
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
    [kfit, cfit] = kfit_2t5pk(tac, cp, scant, kinit, opt, wt, cwb);
    toc;
    figure,plot(t,tac(:,1), '-.o', t, cfit(:,1))
    disp('---- estimate Ki directly -----')
    [(mean(kfit,2)-k0)./k0*100, std(kfit,0,2)./k0*100]
    
    
    disp('---- estimate k2 --------------')
    opt.UpperBound = [0.3 2 2 2 0.5];
    tic;
    [kfit, cfit] = kfit_2t5p(tac, cp, scant, kinit, opt, wt, cwb);
    toc;
    figure,plot(t,tac(:,1), '-.o', t, cfit(:,1))
    k0(3) = (k0(2)/k0(3)-1)*k0(4);
    [(mean(kfit,2)-k0)./k0*100, std(kfit,0,2)./k0*100]
    
    ki0 = mic2mac(k0(2:end),'ki');
    ki  =  mic2mac(kfit(2:end,:),'ki');
    [(mean(ki,2)-ki0)./ki0*100, std(ki,0,2)./ki0*100]
    
end



clear; clc;

% scan time setting
dt = [ones(4,1)*20;  ones(4,1)*40;  ones(4,1)*60; ones(4,1)*180; ones(14,1)*300;];
scant = [cumsum([0; dt(1:end-1)]), cumsum(dt)];
t = mean(scant,2);

% parameters
k0 = [0.05 0.11 0.19 0.10 0.009]'; % = [fv k1 k2 k3 k4]

% optional setting
opt.Decay = log(2)/109.8;
opt.TimeStep = 1.0;
opt.LowerBound = [0.0 0 0 0 0  ];
opt.UpperBound = [0.3 2 2 2 0.5];
opt.MaxIter = 2;
opt.PrmSens = [1 1 1 1 1];
kin = [0.05 0.01 0.01 0.01 0.01]';

% input function
cp = feng(t);
cwb = cp;
cp = cp.*exp(-opt.Decay*t/60); % for consistance with COMKAT

% test for 2-tissue
tic; 
[c, s] = ktac_2t5p(k0, cp, scant, opt, cwb);
toc;

% test least squares fitting
wt = ones(size(c));
tac = zeros(length(c),500001);
for i = 1:size(tac,2)
    tac(:,i) = c + 10*sqrt(c./dt).*randn(size(c));
end
tac = max(0, tac);
kin = repmat(kin(:),[1 size(tac,2)]);

if 0
    disp('---- Single in M and Sequential in C ----')
    tic;
    [kfit, cfit] = kfit_2t5p(tac, cp, scant, kin, opt, wt, cwb);
    toc;   
end
if 0
    disp('---- Sequential in M and Single in C ----')
    kfit = zeros(size(kin));
    cfit = zeros(size(tac));
    tic;
    for i = 1:size(tac,2)
        [kfit_i, cfit_i] = kfit_2t5p(tac(:,i), cp, scant, kin(:,i), opt, wt, cwb);
        kfit(:,i) = kfit_i;
        cfit(:,i) = cfit_i; 
    end
toc;
end
if 0
    disp('---- Parallel in M and Single in C ----')
    tic;
    matlabpool
    parfor i = 1:size(tac,2)
        [kfit_i, cfit_i] = kfit_2t5p(tac(:,i), cp, scant, kin(:,i), opt, wt, cwb);
        kfit(:,i) = kfit_i;
        cfit(:,i) = cfit_i; 
    end
    matlabpool close
    toc;
end    
if 1
    disp('---- Parallel in M and Sequential in C ----')
    tic;
    matlabpool
    num_work = matlabpool('size');
    disp(sprintf('---- Total number of workers: %d', num_work));
    num_pix = size(tac,2);
    N = ceil(size(tac,2)/num_work);
    if num_pix<N*num_work
        nn = num_pix + ( 1:(N*num_work-num_pix) );
        tac(:,nn) = repmat(tac(:,1), [1 length(nn)]);
        kin(:,nn) = repmat(kin(:,1), [1 length(nn)]);
    else
        nn = [];
    end
    kfitN = zeros(N*size(kin,1), num_work);
    cfitN = zeros(N*size(tac,1), num_work);
    parfor i = 1:num_work
        ii = (1:N) + (i-1)*N;
        [kfit_i, cfit_i] = psnkfit_2t5p(tac(:,ii), cp, scant, kin(:,ii), opt, dt, cwb);
        kfitN(:,i) = kfit_i(:);
        cfitN(:,i) = cfit_i(:); 
    end
    matlabpool close
    kfit = reshape(kfitN(:), [size(kin,1) N*num_work]);
    cfit = reshape(cfitN(:), [size(tac,1) N*num_work]);
    tac(:,nn) = [];
    kin(:,nn) = [];
    kfit(:,nn) = [];
    cfit(:,nn) = [];
    toc;
end
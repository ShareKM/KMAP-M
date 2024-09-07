clear; clc

% scan time setting
dt = [ones(6,1)*5;  ones(9,1)*10;  ones(6,1)*30;];
scant = [cumsum([0; dt(1:end-1)]), cumsum(dt)];
t = mean(scant,2);

% optional setting
opt.Decay = log(2)/20.4;
opt.TimeStep = 1.0;
opt.LowerBound = [0 0 0];
opt.UpperBound = [0.5 2 2];
opt.MaxIter = 2;
opt.PrmSens = [1 1 1]';

% input function
Cp = gamma3(t);
Cwb = cumsum(dt/60.*Cp)*0.5;

% generate samples
fv = 0.05; k1 = 0.40; k2 = 0.42; 
c = ktac_1t3p([fv k1 k2]', Cp, scant, opt, Cwb);
for i = 1:100
    ct(:,i) = max(0, c + 20*diag(1./dt)*randn(size(c)) );
end
num_c = size(ct,2);
for i = 1:2
    ii =  num_c + i;
    cp = max(0, Cp+ 50*diag(1./dt)*randn(size(Cp)) );
    ct(:, ii) = (1-0.05)*cp + 0.05*Cwb;
    broi(ii) = 10;
end
broi = logical(broi(:));

% fitting
wt = ones(size(c));
kinit = [0.01 0.1 0.1]';
kinit = repmat(kinit(:),[1 size(ct,2)]);
Cp0 = mean(ct(:,broi),2);
tic;
Cpinit = Cp0;
for i = 1:10
    [k, cfit, Cpnew] = kpim_1t3p(ct, Cpinit, scant, kinit, opt, wt, Cwb, broi, 1, 20);
    Cpinit = Cpnew;
end
toc;
figure, plot(t,Cp,t,Cp0,t,Cpnew, t, mean(ct(:,broi),2))

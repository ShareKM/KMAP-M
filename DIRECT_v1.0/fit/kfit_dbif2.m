function [kfit, cfit] = kfit_dbif2_bf(ct, cp, scant, k0, opt, wt, cwb)
%--------------------------------------------------------------------------
% Fit time activity curves using the two-tissue compartmental model with 
% five kinetic parameters plus dispersion, fa and two time delays in the dual-blood
% input function. The fitting algorithm is the brut-force approach.
%
% The algorithm is developed by Yiran Wang
%
% Guobao Wang @ 4-14-2020
%
%--------------------------------------------------------------------------

% check nan
if any(isnan(ct(:))) | any(isnan(k0(:))) | any(isnan(cp(:)))
    error('NaN value in input');
end

% check size mismatch 
if size(ct,1)~=size(scant,1)
    error('Mismatched frame size.');
end

% make sure the unit of scan time is second
if size(scant,2)~=2
    error('Incorrect scan time input.')
end
if max(scant(:))<180
    %error('The scan time requires to be in seconds!');
end

% initial kinetic parameters  [va K1 k2 k3 k4 ka fa]'
if isempty(k0)
    k0 = [0.05 0.1 0.1 0.1 0.01 1 0];
end
if isvector(k0)
    k0 = repmat(k0(:), [1,size(ct,2)]);
end
kinit = zeros(7,size(k0,2));
if size(k0,1)==7
    kinit(1:size(k0,1),:) = k0;
else
    error('Unmatched size of kinetic parameter input.')
end

% the option paramters
if nargin<5
    opt = setkopt;
end
if isempty(opt.Decay)
    dk = 0;
else
    dk = opt.Decay;
end
if isempty(opt.LowerBound)
    lb = [0, 0, 0, 0, 0, 0, 0];
elseif length(opt.LowerBound)==7
    lb = opt.LowerBound;
else
    error('Size of lower bounds is 7');
end
if isempty(opt.UpperBound)
    ub = [1.0, 10, 5.0, 1.0, 0.1, 100, 1];
elseif length(opt.UpperBound)==7
    ub = opt.UpperBound;
else
    error('Size of upper bounds is 7');
end
if isempty(opt.PrmSens)
    ps = [1, 1, 1, 1, 1, 1, 1];
elseif length(opt.PrmSens)==7
    ps = opt.PrmSens;
else
    error('Size of active label is 7');
end
if isempty(opt.MaxIter)
    maxit = 20;
else
    maxit = opt.MaxIter;
end
if isempty(opt.TimeStep)
    td = 1.0;
else
    td = opt.TimeStep;
end

% frame weights
if nargin<6 | isempty(wt)
    wt = ones(size(ct,1),1);
elseif isvector(wt)
    wt = repmat(wt(:),[1,size(ct,2)]);
end

% whole blood volume
if nargin<7 | isempty(cwb)
    if length(opt.PbrParam)==3
        pbrp = opt.PbrParam;
    else
        pbrp = [1 0 0];
    end
    if length(cp)==size(scant,1)
        t = mean(scant,2);
    elseif length(cp)==scant(end,2)
        t = [1:length(cp)]';
    end    
    pbr = pbrp(1) * exp( -pbrp(2)*t/60) + pbrp(3);
    cwb = cp ./ pbr;
end

% input function
Ca = cp;
if length(cp)<scant(end,2)
    cp = finesample(scant, cp, td);
end
if length(cwb)<scant(end,2)
    cwb = finesample(scant, cwb, td);
end

% check nan
if any(isnan(ct(:))) | any(isnan(cp(:))) | any(isnan(cwb(:))) | any(isnan(kinit(:)))
    error('NaN value in input');
end

%% search range
ka_s = logspace(log10(lb(6)),log10(ub(6)),50);
fa_s = lb(7):0.02:ub(7);
td_s = 0:2:40;

% Portal vein input
B = zeros(length(cp), length(ka_s)*length(fa_s)*length(td_s)^2);
P = zeros(4,size(B,2));
n = 1;
for i = 1:length(ka_s)
    ka = ka_s(i);
    prm = [0 ka ka 0 0];
    Cpv = ktac_2t5p_mex(prm(:), scant, cp, cwb, dk, td);
    for j = 1:length(td_s)
        delay_PV = td_s(j);
        Cpv_delay = finesample_delay(scant, Cpv, [], [], delay_PV);
        for k = 1:length(fa_s)
           	fa = fa_s(k);
            for l = 1:length(td_s)
                delay_HA = td_s(l);
                Ca_delay = finesample_delay(scant, Ca, [], [], delay_HA);
                B(:,n) = (1-fa)*Cpv_delay + fa*Ca_delay;
                P(:,n) = [ka, delay_PV, fa, delay_HA]';
                n = n + 1;
            end
        end
    end
end

%% fitting
parfor n = 1:size(P,2)
    bt = B(:,n);
    [ks(:,n),cs(:,n)] = kfit_2t5p_mex(ct, wt, scant, bt, bt, dk, kinit(1:5), lb(1:5), ub(1:5), ps(1:5), maxit, td);
end
                
mse = sum(diag(wt)*(cs -repmat(ct,[1 size(cs,2)])).^2,1);
[tmp, I] = min(mse);
kfit = [ks(:,I); P(:,I)];
cfit = cs(:,I);

    
    

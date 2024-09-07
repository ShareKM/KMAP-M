function [kfit, cfit] = gsn_2t5p_lma(kin, cp, opt, ct, wt)
%--------------------------------------------------------------------------
% Fit time acitivity curves to nonliear kinetic model using the Gaussian 
% likelihood function (or eqivalently the weighted least squares)
%

% check the option terms
if nargin<3 | isempty(opt)
    setopt_fit;
end
if isempty(opt.Decay)
    dk = 0;
else
    dk = opt.Decay;
end
if isempty(opt.TimeStep)
    td = 1.0;
else
    td = opt.TimeStep;
end
if isempty(opt.ScanTime)
    error('no input for scan time')
else
    scant = opt.ScanTime;
end
if isempty(opt.MaxIter)
    maxit = 20;
else
    maxit = opt.MaxIter;
end

% blood input function
t = mean(scant,2);
pbr = p2blood_ratio(t, opt.PbrPar);
cp  = cp .* exp(-t/60*dk);
Cwb = cp ./ pbr;
B = interp_psf(scant, td, 'linear');
cin = B * cp;
cwb = B * Cwb;

% kinetic model
fitfun = str2func('kfit_2t5p_mex');
num_par = 5;

% Bound constraints
if isempty(opt.LowerBound)
    lb = zeros(num_par,1);
elseif length(opt.LowerBound)==num_par
    lb = opt.LowerBound;
else
    error(sprintf('Size of lower bounds is %d', num_par));
end
if isempty(opt.UpperBound)
    ub = ones(num_par,1);
elseif length(opt.UpperBound)==num_par
    ub = opt.UpperBound;
else
    error(sprintf('Size of upper bounds is %d', num_par));
end
if isempty(opt.ParSens)
    ps = ones(num_par,1);
elseif length(opt.ParSens)==num_par
    ps = opt.ParSens;
else
    error(sprintf('Size of active labels is %d', num_par));
end

% initial kinetic parameters        
if isempty(kin)
    kin = 0.01*ones(num_par,1);
end
if isvector(kin)
    kin = repmat(kin(:), [1,size(ct,2)]);
end
if size(kin,1)~=num_par
    error('mismatched size of kinetic parameters');
end

% frame weights
if nargin<5 | isempty(wt)
    wt = ones(size(ct,1),1);
elseif isvector(wt)
    wt = repmat(wt(:),[1,size(ct,2)]);
end
if size(wt,1)~=size(ct,1)
    error('mismatched sizes of TAC and weighting factors');
end

% check NaN
if any(isnan(kin(:))) | any(isnan(cp)) | any(isnan(ct(:))) | any(isnan(wt(:)))
    error('NaN value in input');
end

% fit time activity curves
[kfit,cfit] = fitfun(ct, wt, scant, cin, cwb, dk, kin, lb, ub, ps, maxit, td);

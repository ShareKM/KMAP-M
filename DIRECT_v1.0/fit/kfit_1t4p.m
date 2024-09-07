function [k, cfit] = kfit_1t4p(ct, cp, scant, k0, opt, wt, cwb)
%--------------------------------------------------------------------------
% Fit time activity curves using the one-tissue compartmental model with 
% four parameters. Gaussian noise model is used.
%
% Guobao Wang @ 6-10-2014
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

% initial kinetic parameters  [va vb K1 k2]'
if isempty(k0)
    k0 = [0.01 0.01 0.01 0.01];
end
if isvector(k0)
    k0 = repmat(k0(:), [1,size(ct,2)]);
end
kinit = zeros(4,size(k0,2));
if size(k0,1)==4
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
    lb = [0, 0, 0, 0];
elseif length(opt.LowerBound)==4
    lb = opt.LowerBound;
else
    error('Size of lower bounds is 4');
end
if isempty(opt.UpperBound)
    ub = [1.0, 1.0, 10, 1.0];
elseif length(opt.UpperBound)==4
    ub = opt.UpperBound;
else
    error('Size of upper bounds is 4');
end
if isempty(opt.PrmSens)
    ps = [1, 1, 1, 1];
elseif length(opt.PrmSens)==4
    ps = opt.PrmSens;
else
    error('Size of active label is 4');
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
    cwb = cp;
end

% input function
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

% fit time activity curves
[k,cfit] = kfit_1t4p_mex_omp(ct, wt, scant, cp, cwb, dk, kinit, lb, ub, ...
                         ps, maxit, td);




   
    
    

function [k, cfit] = kfit_2t5p(ct, cp, scant, k0, opt, wt, cwb)
%--------------------------------------------------------------------------
% Fit time activity curves using the two-tissue compartmental model with 
% five parameters. Gaussian noise model is used.
%
% Guobao Wang @ 12-10-2009
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
% if max(scant(:))<180
%     error('The scan time requires to be in seconds!');
% end

% initial kinetic parameters  [vb K1 k2 k3 k4]'
if isempty(k0)
    k0 = [0.05 0.1 0.1 0.1 0.01];
end
if isvector(k0)
    k0 = repmat(k0(:), [1,size(ct,2)]);
end
kinit = zeros(5,size(k0,2));
if size(k0,1)==5
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
    lb = [0, 0, 0, 0, 0];
elseif length(opt.LowerBound)==5
    lb = opt.LowerBound;
else
    error('Size of lower bounds is 5');
end
if isempty(opt.UpperBound)
    ub = [1.0, 2.0, 2.0, 1.0, 1.0];
elseif length(opt.UpperBound)==5
    ub = opt.UpperBound;
else
    error('Size of upper bounds is 5');
end
if isempty(opt.PrmSens)
    ps = [1, 1, 1, 1, 1];
elseif length(opt.PrmSens)==5
    ps = opt.PrmSens;
else
    error('Size of active label is 5');
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
end
if isvector(wt)
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
ceps = .1%mean(ct(:))*1e-1;
for it = 1:round(maxit/20)
    if it==1
        wt_i = wt;
    else
        wt_i = wt./sqrt(wt.*(ct-cfit).^2+ceps);
        kinit = k;
    end
    [k,cfit] = kfit_2t5p_mex_omp(ct, wt_i, scant, cp, cwb, dk, kinit, lb, ub, ...
                         ps, 20, td);
    
end


    
    

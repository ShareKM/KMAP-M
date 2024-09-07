function [kfit, cfit] = psnfit_srtm(kin, Cr, opt, ct, wt)
%--------------------------------------------------------------------------
% Fit time acitivity curves to the simplified reference tissue model using 
% the Poisson likelihood function
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

% interpolation matrix
B = interp_psf(scant, td, 'linear');
cr = B * Cr;

% kinetic model
switch opt.KinType
    case 'srtm'
        tacfun = str2func('ktac_srtm_mex');
        fitfun = str2func('kfit_srtm_mex');
        num_par = 3;
    otherwise
        error('unknown kinetic model')
end

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
    error(sprintf('Size of upper bovbunds is %d', num_par));
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

% adapt for c code
lb = [0; lb(:)];
ub = [0; ub(:)];
ps = [0; ps(:)];
kin = [zeros(1,size(kin,2)); kin];

% frame weights
if nargin<5 | isempty(wt)
    wt = ( scant(:,2) - scant(:,1) ) / 60;
end
if isvector(wt)
    wt = repmat(wt(:),[1,size(ct,2)]);
end
if size(wt,1)~=size(ct,1)
    error('mismatched sizes of TAC and weighting factors');
end

% check NaN
if any(isnan(kin(:))) | any(isnan(cr)) | any(isnan(ct(:))) | any(isnan(wt(:)))
    error('NaN value in input');
end

% fit time activity curves using Poisson likelihood -----------------------

% initialize
num_frm = size(ct,1);
R1 = kin(2,:);
cfit = tacfun(kin, scant, cr, cr, dk, td);
cc_fit = cfit - Cr*R1;

% fit time activity curves    
subit = 2;
ceps = mean(ct(:))*1e-9;
for it = 1:maxit
    
    % estimate linear coefficients
    for i = 1:5
        R1 = R1 ./ (Cr'*wt) .* ( Cr' * ( wt .* ct ./ cfit ) );
        R1 = max(min(ub(2), R1),lb(2));
        idx = R1 > ( kin(3,:)+1 - 1e-4 );
        %R1(idx) = kin(3,idx)+1 - 1e-4;
        bt = Cr * R1;
        cfit = cc_fit + bt;
    end
    
    % estimate the TACs (transformed into least squares)
    nt = poisson_curv(ct, [], cc_fit, bt, cfit, 'oc');
    cfit(cfit==0) = ceps;
    cc = cc_fit + ( ct./cfit - 1 ) ./ nt;
    cc = max(0, cc);
    c  = cc + bt;
    
    % least squares fitting
    kin(2,:) = R1;
    ps(2) = 0;
    [kin, cfit] = fitfun(c, wt.*nt, scant, cr, cr, dk, kin, lb, ub, ps, subit, td);
    
    % update the TAC
    cc_fit = cfit - bt;
    
end
kfit = kin(2:end,:);


%--------------------------------------------------------------------------
function ni = poisson_curv(yi, ai, li, ri, yb, ctype)
%--------------------------------------------------------------------------
% Transfer the Poisson likelihood function into a paraboloidal function. 
% This code is adapted from Jeff Fessler's image reconstruction toolbox. 
% See his MIC98 conference paper for more details.
%

% check input
if isempty(ai)
    ai = ones(size(yi));
end
if isempty(yb) | nargin<5
    yb = ai .* li + ri;
end
if isempty(ctype) | nargin<6
    ctype = 'oc';
end

% precomputed curvature
if strcmp(ctype, 'pc')     
	li = (yi - ri) ./ ai;
	li = max(li, 0);
	yb = ai .* li + ri;
	if any(yb(:)<=0 & yi(:)>0)
        error 'model mismatch'
    end
	ni = zeros(size(li));
	idx = li > 0;
	ni(idx) = ai(idx).^2 ./ max(yi(idx),1);
	idx = li<=0 & ri>0;
	ni(idx) = (ai(idx) ./ ri(idx)).^2 .* yi(idx);
    
% optimal curvature
elseif strcmp(ctype, 'oc')  
    
    % small nonzero value for yi=0
    ni = 0.001 * (ai ./ ri).^2;     

    % curvature at l=0 and yi>0
    ni_max = yi .* (ai ./ ri).^2;	
    tmp = log(yb ./ ri) - ai.*li ./ yb;
    iy = yi ~= 0;
    il = li < 0.1 * ri ./ ai;
    idx = iy & il;
    ni(idx) = ni_max(idx);

    % curvature at l>0 and yi>0
    idx = iy & ~il;
    ni(idx) = 2 * yi(idx) ./ li(idx).^2 .* tmp(idx);

    if any(ni <= 0)
        idx = ni(:)<=0;
        [yi(idx) li(idx) ri(idx)]
        error 'nonpositive ni'
    end
    if any((ni > ni_max) & yi)
        plot([ni_max(:) ni(:) ni(:)>ni_max(:)])
        error 'large ni'
    end
else
    error('unknown curvature type');
end

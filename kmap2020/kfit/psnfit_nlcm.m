function [kfit, cfit] = psnfit_nlcm(kin, cp, opt, ct, wt)
%--------------------------------------------------------------------------
% Fit time acitivity curves to nonlinear kinetic model using the Poisson 
% likelihood function
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
if opt.ParSens(1)==1
    flagvb = 1;
else
    flagvb = 0;
end

% blood input function
pbr = p2blood_ratio(mean(scant,2), opt.PbrPar);
Cwb = cp ./ pbr;
B = interp_psf(scant, td, 'linear');
bld = B * cp;
cwb = B * Cwb;

% kinetic model
switch opt.KinType
    case '1t3p'
        tacfun = str2func('ktac_1t3p_mex');
        fitfun = str2func('kfit_1t3p_mex');
        num_par = 3;
    case '2t5p'
        tacfun = str2func('ktac_2t5p_mex');
        fitfun = str2func('kfit_2t5p_mex');
        num_par = 5;
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
    wt = scant(:,2) - scant(:,1);
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

% fit time activity curves using Poisson likelihood -----------------------

% initialize
num_frm = size(ct,1);
vb = kin(1,:);
ri = Cwb * vb;
ai = 1 - vb;
ps(1) = 0;
kin(1,:) = 0.0;
cc_fit = tacfun(kin, scant, bld, cwb, dk, td);
cfit = repmat(ai,[num_frm 1]) .* cc_fit + ri;

% fit time activity curves    
subit = 2;
ceps = mean(ct(:))*1e-9;
for it = 1:maxit
    
    % estimate linear coefficients K1(1-vb) and vb
    for i = 1:5
        if flagvb
            vb = vb ./ (Cwb'*wt) .* ( Cwb' * ( wt .* ct ./ cfit ) );
            vb = min(ub(1), vb);
            ri = Cwb * vb;
        end
        ai = ai ./ sum(cc_fit.*wt,1) .* sum(cc_fit.*wt.*ct./cfit,1);
        cfit = repmat(ai,[num_frm 1]) .* cc_fit + ri;
    end
    
    % transformed into least squares
    ai = repmat(ai,[num_frm 1]);
    nt = poisson_curv(ct, ai, cc_fit, ri, cfit, 'pc');
    cfit(cfit==0) = ceps;
    cc = cc_fit + ai./nt .* ( ct./cfit - 1 );
    cc = max(0, cc);
    
    % least squares fitting with vb = 0
    [kfit,cc_fit] = fitfun(cc, wt.*nt, scant, bld, cwb, dk, kin, lb, ub, ps, subit, td);
    kin = kfit; 
    
    % update the total activity
    cfit = ai .* cc_fit + ri;
    ai = ai(1,:);
    
end

% pass vb to k vector
kfit(1,:) = vb;
kfit(2,:) = kfit(2,:) .* ai ./ (1-vb);


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

    if any(ni <= 0), error 'nonpositive ni', end
    if any((ni > ni_max) & yi)
        plot([ni_max(:) ni(:) ni(:)>ni_max(:)])
        error 'large ni'
    end
else
    error('unknown curvature type');
end

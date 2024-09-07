function [kfit, cfit] = psn_1t3p1_lma1(kin, cp, opt, ct, wt)
%--------------------------------------------------------------------------
% Fit time acitivity curves to nonlinear kinetic model using the Poisson 
% likelihood function and the Levenberg-Marquardt algorithm
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
t = mean(scant,2);
pbr = p2blood_ratio(t, opt.PbrPar);
cp  = cp .* exp(-t/60*dk);
Cwb = cp ./ pbr;
B = interp_psf(scant, td, 'linear');
cin = B * cp;
cwb = B * Cwb;

% kinetic model
tacfun = str2func('ktac_1t3p1_mex');
fitfun = str2func('kfit_1t3p1_mex');
num_par = 3;

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
    wt = ( scant(:,2) - scant(:,1) ) / 60;
end
if isvector(wt)
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
vb(vb<=0) = 1e-3;
ri = Cwb * vb;
cfit = tacfun(kin, scant, cin, cwb, dk, td);
cc_fit = cfit - ri;

% fit time activity curves    
subit = 2;
ceps = mean(ct(:))*1e-9;
for it = 1:maxit
    
    % estimate linear coefficients vb
    if flagvb
        va = 1;
        for i = 1:5        
            vb = vb ./ (Cwb'*wt) .* ( Cwb' * ( wt .* ct ./ cfit ) );
            va = va ./ sum(cc_fit.*wt, 1) .* sum(cc_fit.*wt.*ct./cfit, 1);
            ri = Cwb * vb;
            Va = repmat(va,[size(cc_fit,1) 1]);
            cfit = Va .* cc_fit + ri;
        end
        kin(1,:) = vb;
        kin(2,:) = kin(2,:).*va;
        cc_fit = Va .* cc_fit;
    end    
    
    % transformed into least squares
    nt = poisson_curv(ct, [], cc_fit, ri, cfit, 'oc');
    cfit(cfit==0) = ceps;
    cc = cfit + ( ct./cfit - 1 ) ./ nt;
    cc = max(ceps, cc);
    
    % least squares fitting
    ps(1) = 0;
    [kfit,cfit] = fitfun(cc, wt.*nt, scant, cin, cwb, dk, kin, lb, ub, ps, subit, td);
    kin = kfit;
    
    % update
    cc_fit = cfit - ri;
    
end


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

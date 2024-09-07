function [k, xt, ks, L, cpucost] = otsp(yt, a, P, xt, rt, maxit, beta, imgsize, k, cp, opt, step)

%--------------------------------------------------------------------------
% check inputs for reconstruction

% total frame number
num_frame = size(yt,2);

% check inputs
if isempty(a)
    a = ones(size(P,1),1);
end
if isempty(rt)
    rt = repmat(mean(yt,1)*1e-8,[size(yt,1),1]);
end
if nargin<6 | isempty(beta)
    beta = 0;
end
if nargin<7
    imgsize = [];
end
if nargin<12
    step = [];
end
if beta>0 & isempty(imgsize)
    error('image size is underdetermined');
end
if nargin<9
    kflag = 0;
else
    kflag = 1;
end

% sensitivity map
at = repmat(a,[1,num_frame]);
p1 = repmat(sum(P,2),[1,num_frame]);
pj = P'*a;
mask = pj>0;
pj(~mask) = mean(pj)*1e-9;
wt = repmat(pj,[1 num_frame]);

%--------------------------------------------------------------------------
% check inputs for fitting
num_pixel = size(P,2);
if kflag
    if isempty(k)
        k = [0.05 0.01 0.01 0.01 0.01]';
    end
    if size(k,2)<num_pixel
        k = repmat(k(:,1), [1 num_pixel]);
    end
    k(:,~mask) = 0;

    % time duration in minute
    scant = opt.ScanTime;
    dt = ( scant(:,2) - scant(:,1) ) / 60; 

    % initialization
    if isempty(xt)
        cc = ktac_2t5p(k, cp, scant, opt);
    else
        ct = diag(1./dt) * xt';
        cc = ct;
        pp = opt; pp.MaxIter = 20;
        [k(:,mask), cc(:,mask)] = kfit_2t5p(ct(:,mask), cp, scant, k(:,mask), pp, dt);  
    end
    cc(:,~mask) = 0;
    xt = cc' * diag(dt);
else
    xt = ones(num_pixel, num_frame) * 0.1;
end

ks = [];
if nargout>2 & not(isempty(step))
    ks = zeros(size(k,1),size(k,2),floor(maxit/step));
end

%--------------------------------------------------------------------------
% iterative loop
if nargin>4
    cpucost = zeros(maxit,2);
end
for it = 1:maxit   
    
    % save data
    if not(isempty(step)) & nargout>2 & rem(it,step)==0 & kflag
        ks(:,:,it/step) = k;
    end  
    if nargin>4
        cpucost(it,1) = cputime;
    end
    
    % sp reconstruction
    yb = at.*(P*xt) + rt;
    yy = yt ./ ( yb + mean(yt(:))*1e-9 );
    yy( yb==0 & yt==0 ) = 1;    
    gt = P' * (at.*(yy-1));
    nt = poisson_curvature(yt, at, [], rt, yb, 'oc');   
    wt = P'*(nt.*p1);
    xt_sp = xt + gt./wt;
    idx = wt==0;
    xt_sp(idx) = 0;    
    
    % smoothed estimate 
    if beta>0
        [xt_reg, wt_reg] = smoother(xt, imgsize);
        xt_reg(~mask,:) = 0;
        wt_reg(~mask,:) = 0;
    else
        xt_reg = 0;
        wt_reg = 0;
    end
    
    % objective function value
    if nargout > 3 
        idx = yb>0;
        tmp = yt(idx) .* log(yb(idx)) - yb(idx);
        L(it) = sum(tmp(:)); 
        if beta>0
            L(it) = L(it) - sum( beta .* sum(wt_reg.*(xt-xt_reg).*(xt/2), 1) );
        end
    end
    
    % fitting
    if nargin>4
        cpucost(it,2) = cputime;
    end
    if kflag
        ct = diag(1./dt) * ( (wt.*xt_sp + beta*wt_reg.*xt_reg)./(wt+beta*wt_reg) )';
        wt = diag(dt.^2) * ( wt + beta*wt_reg )';
        ct(wt==0) = 0; %ct = max(0,ct);
        [k(:,mask), cc(:,mask)] = kfit_2t5p(ct(:,mask), cp, scant, k(:,mask), opt, wt(:,mask));
        xt = cc' * diag(dt);
    else
        wt = wt + beta*wt_reg;
        xt = (wx.*xt_sp + beta*wt_reg.*xt_reg)./wt;
        xt(wt==0) = 0;
        xt(xt<0) = 0;
    end
    
end


%--------------------------------------------------------------------------
function [xt_flt, wgt] = smoother(xt, imgsize)
%--------------------------------------------------------------------------
% construct the 2D or 3D filter template and filter the images xt
%

% define the filter
if length(imgsize)==2
    a = 1/sqrt(2);
    weight = [a, 1, a;
              1, 0, 1;
              a, 1, a];
elseif length(imgsize)==3
    a = 1/sqrt(2);
    b = 1/sqrt(3);
    weight(:,:,1) = [b, a, b;
                     a, 1, a;
                     b, a, b];
    weight(:,:,2) = [a, 1, a;
                     1, 0, 1;
                     a, 1, a];
    weight(:,:,3) = weight(:,:,1);
else
    error('4 or higher-D not included yet')
end
template = weight;
wgt = sum(weight(:));
template(weight==0) = wgt;
template = template / (2*wgt);

% invariant flitering
xt_flt = zeros(size(xt));
for m = 1:size(xt,2)
    xf = convn(reshape(xt(:,m),imgsize), template, 'same');
    xt_flt(:,m) = xf(:);
end
wgt = wgt * ones(size(xt_flt));

%--------------------------------------------------------------------------
function ni = poisson_curvature(yi, ai, li, ri, yb, ctype)
%--------------------------------------------------------------------------
% Transfer the Poisson likelihood function into a paraboloidal function. 
% This code is adapted from Jeff Fessler's image reconstruction toolbox. 
% See his MIC98 conference paper for more details.
%

% check input
if isempty(ai)
    ai = ones(size(yi));
end
if isempty(li)
    li = (yb-ri)./ai;
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

function [kp, out, Phi] = eml_ot_em(yt, nt, G, Gopt, rt, maxit, beta, R, cp, kopt, kp0)
%--------------------------------------------------------------------------
% The optimization transfer algorithm with the EM surrogate for penalized
% likelihood direct reconstruction of parameteric images from dynamic PET
% data. Reference: 
%
%   G. B. Wang, J. Qi, An optimization transfer algorithm for nonlinear 
%   parametric image reconstruction from dynamic PET data, IEEE Tran. on 
%   Medical Imaging, 31(10): 1977-1988, 2012.
%
%--------------------------------------------------------------------------
% INPUT:
%   yt      dynamic sinogram, each column of which is a sinogram.
%   nt      multiplicative factor (normalization, attenuation)
%   G       system matrix
%   Gopt    option set for G, can be empty if G is a matlab sparse matrix
%   rt      addtive factor (randoms and scatters)
%   maxit   maximum iteration number
%   beta    regularization parameter
%   R       regularization operator, can be set to empty if beta=0
%   cp      blood or reference input function
%   kopt    option set for kinetic model
%   kp0     initial estimate of kinetic paraemters
%
%   yt, nt, rt have the same size of [Ni Nm] where Ni is the total number
%   of sinogram bins and Nm is the total number of time frames. kp0 is of
%   size [Nj Np] where Nj is the total number of pixels in image and Np is 
%   the total number of kinetic parameters in the kinetic model.
%
% OUTPUT
%   kp      Estimate of kinetic parametric images
%   out     output
%   Phi     objective function
%
%--------------------------------------------------------------------------
% Programmer: Guobao Wang @ UC Davis, Qi Lab, 10-05-2007
% Last Modified: 06-28-2013
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% check inputs for reconstruction

% total frame number
numfrm = size(yt,2);
numpix = prod(Gopt.imgsiz);

% preprocess sinograms
[yt, rt, nt] = sino_preprocess(yt, rt, nt, numfrm);

% set Gopt
Gopt = setGopt(nt, G, Gopt);

% maximum iteration number
if isempty(maxit)
    maxit = 50;
end

% regularization
if isempty(beta) | nargin<8
    beta = 0;
end
R = setRopt(R);
if beta>0
    C1 = sum(abs(R.C),2);
end

% kinetic model


% initialization


x    = max(mean(x0(:))*1e-9,x0(:)); x(~Gopt.mask) = 0;
yb   = ni.*proj_forw(G, Gopt, x) + ri;
cx   = R.C*x;
yeps = mean(yi(:))*1e-9;
wx   = Gopt.sens;

% sensitivity map
at = repmat(a,[1,num_frame]);
pj = P'*a;
fov = pj>0;
pj(~fov) = mean(pj)*1e-9;
wt = repmat(pj,[1 num_frame]);
mask = mask(:) & fov;


%--------------------------------------------------------------------------
% check inputs for fitting

if isempty(k)
    k = [0.05 0.01 0.01 0.01 0.01]';
end
num_pixel = size(P,2);
if size(k,2)<num_pixel
    k = repmat(k(:,1), [1 num_pixel]);
end
k(:,~mask) = 0;

% time duration in minute
scant = opt.ScanTime;
dt = ( scant(:,2) - scant(:,1) ) / 60; 

% plasma-to-blood ratio function
if isempty(opt)
    opt = setkopt;       
end
if length(opt.PbrParam)==3
    pbrp = opt.PbrParam;
else    
    pbrp = [1 0 0];
end
pbr = pbrp(1) * exp( -pbrp(2) * mean(scant,2)/60 ) + pbrp(3);
bpr = 1./pbr;
        
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

ks = [];
if nargout>2 & not(isempty(step))
    ks = zeros(size(k,1),size(k,2),floor(maxit/step));
end
if nargin>4
    cpucost = zeros(maxit,2);
end
%--------------------------------------------------------------------------
% iterative loop

for it = 1:maxit
       
    % save data
    if not(isempty(step)) & nargout>2 & rem(it,step)==0
        ks(:,:,it/step) = k;
    end
    if nargin>4
        cpucost(it,1) = cputime;
    end
    
    % em reconstruction
    yb = at.*(P*(xt+repmat(mk(:)',[size(P,2),1]))) + rt;
    yy = yt ./ ( yb + mean(yt(:))*1e-9 );
    yy( yb==0 & yt==0 ) = 1;    
    xt_em = (xt+repmat(mk(:)',[size(P,2),1])) ./ wt .* ( P' * (at.*yy) );
    xt_em(~mask,:) = 0;
    
    % smoothed estimate 
    if beta>0
        [xt_reg, wt_reg] = smoother(xt, imgsize);
        xt_reg(~mask,:) = 0;
        wt_reg(~mask,:) = 0;
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
    ct = diag(1./dt) * xt_em';
    if beta>0
        bet = beta./wt.*wt_reg*(diag(dt.^2));
        [k(:,mask), cc(:,mask)] = hfit_2t5p(ct(:,mask), cp, scant, k(:,mask), ...
            opt, dt, [], bet(mask,:)', diag(1./dt)*xt_reg(mask,:)', diag(1./dt)*mk(:));
    else
        [k(:,mask), cc(:,mask)] = hfit_2t5p(ct(:,mask), cp, scant, k(:,mask), ...
            opt, dt, [], [], [], diag(1./dt)*mk(:));
    end
    
    % update
    xt = cc' * diag(dt);
    
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

function [varargout] = kfit_srtm_bfm(cr, tac, wt, scant, lam, dk)
%--------------------------------------------------------------------------
% fitsrtm_bfm implements simplified reference tissue modeling (SRTM) for 
% estimating receptor binding parameters using the basis function algorithm.
%
% INPUT
%   cr      tac in the reference region
%   tac     time acitivity curve. It can be a vector or a matrix of which
%           each column is a tac.
%   wt      temporal weighting factors, can be empty
%   scant   [ts te], the start and end time for each frame. if a vector, it
%           is t = (ts+te)/2;
%   lam     vector of rate constants used for the BFM algorithm. lam can be
%           a three-element vector which gives [mink maxk knum].[mink maxk] 
%           is also acceptable which has a defaut knum = 100. 
%   dk      tracer decay constant
%
% OUTPUT
%   R1      the ratio of transport rate in the target region, R1 = K1/K1_r
%   BP      binding potential
%   tacFit  fitted time activity curves
%
% USAGE
%   [R1, BP, k2] = kfit_srtm_bfm(tac, cr, wt, scant, lam, dk)
%   [R1, BP, k2, tacFit] = kfit_srtm_bfm(tac, cr, wt, scant, lam, dk)
%   [R1, BP, k2, tacFit, k2r] = kfit_srtm_bfm(tac, cr, wt, scant, lam, dk)
%
% PROGRAMER
%   Guobao Wang @ UC Davis
%   Sept 10, 2008
% 
% Last Updated: 
%   Nov 4, 2009
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% check input
%
if min(size(tac))==1
    tac = tac(:);
end
if isempty(wt)
    wt = ones(size(tac,1),1);
end
W = diag(sqrt(wt));
if size(tac,1)~=length(wt(:))
    error('unmatched size')
end
if size(scant,2)~=2
    error('unmatched size of scant')
end
if length(lam(:))<=3
    if isempty(lam)
        lam = [0.01 1.5 100];
    end
    if length(lam(:))==2
        lam = [lam(:)' 100];
    end
    lam = logspace(log10(lam(1)), log10(lam(2)), lam(3));
end
NumOfBasis = length(lam);
NumOfPixel = size(tac,2);
NumOfFrame = size(tac,1);


%--------------------------------------------------------------------------
% least squares estimation of microparameters
%

% prepare the basis functions
B = GenerateBasisFunctions(scant, cr, lam, dk);

% initialization
R_1 = zeros(size(tac,2),1);
k_2a = zeros(size(tac,2),1);
k_2r = zeros(size(tac,2),1);
M = zeros(NumOfBasis*2,size(tac,1));
for i = 1:NumOfBasis
    A = [cr B(:,i)];
    [Q,R] = qr(W*A);
    M(2*i-1:2*i,:) = R\Q';
end

% The following implementation is fast for parametric imaging
rss = zeros(NumOfBasis,NumOfPixel);
tht = zeros(2,NumOfPixel,NumOfBasis);
for i = 1:NumOfBasis
    tht_i = M(2*i-1:2*i,:) * W * tac;
    tht(:,:,i) = tht_i;
    ct = [cr B(:,i)] * tht_i;
    r = tac - ct;   
    rss(i,:) = sum((W*r).^2, 1);
end
[tmp,i_min] = min(rss);       
for j = 1:NumOfPixel
    tht_j = tht(:,j,i_min(j));
    tht_new(:,j) = tht_j;
end

% calculate the binding parameters
R_1 = tht_new(1,:)';
k_2a = lam(i_min)';
idr = R_1 > 1e-3; 
k_2r(idr) = tht_new(2,idr)' ./ R_1(idr) + k_2a(idr);
k_2 = tht_new(2,:)' + R_1.*lam(i_min)';

%--------------------------------------------------------------------------
% output
%
R_1(R_1<0) = 0;
varargout{1} = R_1;
if nargout > 1
    B_P = R_1 .* k_2r ./ k_2a;
    B_P(B_P<0) = 0;
    varargout{2} = B_P;
end
if nargout > 2
    varargout{3} = k_2;
end
if nargout > 3
    tacFit = zeros(NumOfFrame,NumOfPixel);
    for m = 1:NumOfFrame
        tacFit(m,:) = tht_new(1,:)*cr(m) + tht_new(2,:).*B(m,i_min);
    end
end
if nargout > 3
    varargout{4} = tacFit;
end
if nargout > 4
    varargout{5} = k_2r;
end
    
    
%--------------------------------------------------------------------------
function B = GenerateBasisFunctions(scant,Ca,k,dk)
%--------------------------------------------------------------------------
NumOfBasis = length(k);
NumOfFrame = size(scant,1);

% interpolation
res = 1/60;   % time resolution: 1/60 min = 1 sec
scant = scant / res;
tm = mean(scant,2);
td = 1:scant(end,2);
Ca = interp1(tm,Ca,td,'spline','extrap') * res; 

B = zeros(NumOfFrame,NumOfBasis);
for i = 1:NumOfBasis
    vec = exp_conv(Ca(:),1,(k(i)+dk)*res,td(:));
    for j = 1:NumOfFrame
        tt = im2single([(scant(j,1)+1):scant(j,2)]);
        dt = scant(j,2)-scant(j,1);
        B(j,i) = sum(vec(tt))/dt;
    end
end
    
%--------------------------------------------------------------------------
function convolution = exp_conv(inp,k1,k2,t)
% convolute the exponential function with the input function. This is 
% originated from Dr Gunn's RPM code.
%--------------------------------------------------------------------------
convolution = zeros(1,length(t));
ek2=exp(-k2);
prev=0;
switch ek2
    case 1
        convolution=k1*cumsum(inp);   
    otherwise
        for i=1:length(t)
            prev=prev.*ek2+k1.*inp(i).*(1-ek2)./k2;
            convolution(i)=prev;
        end
end

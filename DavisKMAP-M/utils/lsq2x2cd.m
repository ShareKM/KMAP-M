function [x1, x2, r] = lsq2x2cd(y, a1, a2, w, x1, x2, maxit, vlower, bupper)
%--------------------------------------------------------------------------
% Coordinate descent algorithm for least squares estimation of 2-parameter 
% vector x under the model y = A*x
% NOTE: x1>=0

% the weighting factor
if nargin < 4 | isempty(w)
    w = ones(size(y));
end
if isvector(w)
    w = repmat(w(:), [1, size(y,2)]);
end
if nargin<5 | isempty(x1)
    x1 = zeros(size(y,2),1);
end
if nargin<6 | isempty(x2)
    x2 = zeros(size(y,2),1);
end
if nargin<7 | isempty(maxit)
    maxit = 200;
end
if nargin<8 | isempty(vlower)
    vlower = 0;
end
if nargin<9 | isempty(bupper)
    bupper = 1e3;
end

% run coordinate descent
r  = y - a1.*repmat(x1',[size(a1,1),1]) - a2.*repmat(x2',[size(a2,1),1]);
b1 = sum(a1.^2,1)'; w1 = 1./b1; w1(b1==0)=0;
b2 = sum(a2.^2,1)'; w2 = 1./b2; w2(b2==0)=0;
for it = 1:maxit
    xi = x1;
    x1 = x1 + sum(a1.*r,1)'.*w1;
    x1 = max(vlower,x1);
    r  = r - a1.*repmat(x1-xi,[1 size(a1,1)])';
    xi = x2;
    x2 = x2 + sum(a2.*r,1)'.*w2;
    idx     = abs(x2)>bupper;
    x2(idx) = sign(x2(idx))*bupper;
    r  = r - a2.*repmat(x2-xi,[1 size(a2,1)])';
end

function [x1, x2, r] = lsq2x2nonneg(y, a1, a2, w, x1, x2, maxit)
%--------------------------------------------------------------------------
% Coordinate descent algorithm for least squares estimation of 2-parameter 
% vector x under the model y = A*x
% NOTE: x1>=0

% the weighting factor
if isvector(a1) & length(a1)==size(y,1)
    a1 = repmat(a1,[1 size(y,2)]);
end
if isvector(a2) & length(a2)==size(y,1)
    a2 = repmat(a2,[1 size(y,2)]);
end
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
    maxit = 100;
end

% first run, no nonnegativity on v
[v, b, r] = lsq2x2(y, a1, a2); 
idx    = v>0 & abs(b)<1e3;
x1(idx)= v(idx);
x2(idx)= b(idx);

% second run for nonnegativity on v
idx = ~idx;
if any(idx)
    [x1(idx), x2(idx), r(:,idx)] = lsq2x2cd(y(:,idx), a1(:,idx), a2(:,idx), w(:,idx), x1(idx), x2(idx), maxit);
end

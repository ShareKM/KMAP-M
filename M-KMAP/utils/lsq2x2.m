function [x1, x2, r] = lsq2x2(y, a1, a2, w)
%--------------------------------------------------------------------------
% Least squares estimation of 2-parameter vector x under the model y = A*x
%

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

% calculate the determinant
a11 = sum(a1.*w.*a1, 1);
a12 = sum(a1.*w.*a2, 1);
a22 = sum(a2.*w.*a2, 1);
adet = a11.*a22 - a12.*a12; 
g1 = sum(a1.*w.*y, 1);
g2 = sum(a2.*w.*y, 1);

% the least-square solution
x1 = zeros(size(y,2),1);
x2 = zeros(size(y,2),1);
idd = adet~=0;
x1(idd) = ( a22(idd).*g1(idd) - a12(idd).*g2(idd)) ./ adet(idd);
x2(idd) = (-a12(idd).*g1(idd) + a11(idd).*g2(idd)) ./ adet(idd);

% residual
if nargin>2
    r = y - a1.*repmat(x1',[size(a1,1) 1]) - a2.*repmat(x2',[size(a2,1) 1]);
end

function newK = make_sym(K, maxit)
%--------------------------------------------------------------------------
% An iterative algorithm for enforing a matrix K to become symmetric and 
% normalized. The sum of each column or row is equal to 1.
% 
%   P. Milanfar, Symmetrizing smoothing filters, SIAM Journal of Imaging
%   Science, 6(1), 263-284, 2013
%
% Thanks Dr. Jian Zhou for bringing up this paper to my attention.
% 
% G. Wang @ UC Davis, 12/01/2013
%
if nargin<2
    maxit = 100;
end

N = size(K,1);
r = ones(N,1);
for it = 1:maxit
    c = 1./ (K'*r);
    r = 1./ (K *c);
end
C = spdiags(c, 0, N, N);
R = spdiags(r, 0, N, N);
newK = R * K * C;

newK = ( newK' + newK )/2;


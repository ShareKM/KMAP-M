function [x, out] = tripatlak(cp, ct, scant, maxit, x0)
%--------------------------------------------------------------------------
% patlakplot implements Patlak plot for estimating influx rate parameter of
% irreversible tracer.
%
% INPUT
%   cp      input function
%   ct      time acitivity curve. It can be a vector or a matrix of which
%           each column is a tac.
%   scant   [ts te], the start and end time for each frame
%
% OUTPUT
%   x       parameter estimates
%   out     output information
%
% PROGRAMER
%   Guobao Wang @ UC Davis
%   May 30, 2013
% 
% Last Updated: 
%   May 31, 2013
%
%--------------------------------------------------------------------------

%% check input
[num_frm, num_pix] = size(ct);
if nargin<4 | isempty(maxit)
    maxit = 100;
end
if nargin<5 | isempty(a0)
    x0 = zeros(3,size(ct,2));
end

%% prepare data for the tri-variate Patlak plot
scant = scant / 60;
dt = scant(:,2) - scant(:,1);
sp = cumsum(diag(dt)*cp, 1);
delta = mean(cp(:))*1e-9;
if delta<=0
    delta = 1e-9;
end
cp(cp==0) = delta;
A(:,1) = sp;
A(:,2) = cp;
A(:,3) = 1;

% estimation
w = A'*sum(A,2);
x = x0;
ax = A*x;
for it = 1:maxit
    r = ct - ax;
    out.objf(it) = sum(r(:).^2)/2;
    g = -A'*r;
    xnew = x - diag(1./w)*g;
    dx = xnew - x;
    adx = A*dx;
    a = sum(r.*adx,1)./sum(adx.^2,1); a(isnan(a))=1;   
    x = x + dx.*repmat(a,[size(x,1) 1]);
    ax = ax + adx.*repmat(a,[size(ax,1) 1]);
end
x = max(0,x);
    

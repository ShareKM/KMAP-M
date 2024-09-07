function [k, b, out] = sbpatlak(cp, ct, scant, s0, sbflag, maxit)
%--------------------------------------------------------------------------
% patlakplot implements Patlak plot for estimating influx rate parameter of
% irreversible tracer.
%
% INPUT
%   cp      input function
%   ct      time acitivity curve. It can be a vector or a matrix of which
%           each column is a tac.
%   scant   [ts te], the start and end time for each frame
%   s0      integral of cp from 0 to t*
%   sbflag  indicator if semiblind estimation is used
%
% OUTPUT
%   k       Patlak slope
%   b       Patlak intercept
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
if nargin<4 | isempty(s0)
    s0 = 0;
end
if nargin<5 | isempty(sbflag)
    sbflag = 1;
end
if sbflag==0 & s0==0
    warning('It seems s0=0 is incorrect');
end
if nargin<6 | isempty(maxit)
    maxit = 10;
end
if sbflag==0
    maxit = 1;
end

%% prepare data for the Patlak plot
scant = scant / 60;
dt = scant(:,2) - scant(:,1);
sp = cumsum(diag(dt)*cp, 1);
delta = mean(cp(:))*1e-9;
if delta<=0
    delta = 1e-9;
end
cp(cp==0) = delta;
a1 = sp + s0;
a2 = cp;
[k, b, r] = lsq2x2(ct, a1, a2);
out.objf = zeros(maxit,1);
out.sest = zeros(maxit,1);

%% iterative loop
for it = 1:maxit
    
    % patlak plot
    [k, b, r] = lsq2x2nonneg(ct, a1, a2, [], k, b, 2);
    
    % objective function
    out.objf(it) = sum(r(:).^2)/2;
    out.sest(it) = s0;
    
    % estimation of s0
    if sbflag
        K = repmat(k',[size(r,1) 1]);
        ds = (r(:)'*K(:))/(K(:)'*K(:));
        s0 = s0 + ds;
        a1 = sp + s0;
    end
    
end


    

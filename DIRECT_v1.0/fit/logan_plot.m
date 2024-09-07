function [v, b, x, y, yt] = loganplot(cp, scant, dk, ct, tx)
%--------------------------------------------------------------------------
% loganplot implements the Logan plot for estimating the volume of 
% distribution of reversible tracers.
%
% INPUT
%   cp      input function
%   ct      time acitivity curve. It can be a vector or a matrix of which
%           each column is a tac.
%   scant   [ts te], the start and end time for each frame
%   tx      start time point
%   dk      decay coefficient
%
% OUTPUT
%   v       Logan slope
%   b       Logan intercept
%   x       Logan x coordinate
%   y       Logan y coordinate
%
% PROGRAMER
%   Guobao Wang @ UC Davis
%   Sept 10, 2008
% 
% Last Updated: 
%   04-30-2012
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% check input
%
[num_frm, num_pix] = size(ct);
if isvector(cp)
    cp = repmat(cp(:), [1, num_pix]);
end      
if isempty(dk) 
   dk = 0;
end
if isempty(tx) | nargin < 5
   tx = 1;
end
if tx>num_frm
    error('incorrect start time point')
end

% correct for decay
ck = decaycoef(scant, dk);
cp = diag(ck) * cp; 
ct = diag(ck) * ct; 

% prepare data for the Logan plot
dt = scant(:,2) - scant(:,1);
sp = cumsum(diag(dt)*cp, 1);
st = cumsum(diag(dt)*ct, 1);
delta = mean(ct(:))*1e-9;
if delta<=0
    delta = 1e-9;
end
ct(ct==0) = delta;
a1 = sp ./ ct;
a2 = ones(size(a1));
yt = st ./ ct;

% estimate the slope and intercept
tt = tx:num_frm;    % effective time frames  
[v, b] = leastsquare2x2(yt(tt,:), a1(tt,:), a2(tt,:));

% output the x and y coordinates of the plot
if nargout>3
    x = a1;
    y = x .* repmat(v, [num_frm, 1]) + repmat(b, [num_frm, 1]);
end

%--------------------------------------------------------------------------
function [x1, x2] = leastsquare2x2(y, a1, a2, w)
%--------------------------------------------------------------------------
% Least squares estimation of 2-parameter vector x under the model y = A*x
%

% the weighting factor
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
x1 = zeros(1, size(y,2));
x2 = zeros(1, size(y,2));
idd = adet~=0;
x1(idd) = ( a22(idd).*g1(idd) - a12(idd).*g2(idd)) ./ adet(idd);
x2(idd) = (-a12(idd).*g1(idd) + a11(idd).*g2(idd)) ./ adet(idd);

    
%--------------------------------------------------------------------------    
function cc = decaycoef(scant, dk, c)
%--------------------------------------------------------------------------
% decay correction on time activity curves
%
ts = scant(:,1);
dt = scant(:,2)-scant(:,1);
if dk>0
    c1 = exp( - dk * ts );
    c2 = exp( - dk * ( ts + dt ) );
    cc = dk * dt ./ ( c1 - c2 ); 
else
    cc = ones(size(scant,1),1);
end
if nargin>2
    cc = c .* repmat(cc(:),[1 size(c,2)]);
end

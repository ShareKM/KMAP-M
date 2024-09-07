function [k, b, x, y, yt] = firstpplot(cp, ct, scant, tx, alg)
%--------------------------------------------------------------------------
% implement the first-pass plot for estimating influx rate parameter of
% irreversible tracer.
%
% INPUT
%   cp      input function
%   ct      time acitivity curve. It can be a vector or a matrix of which
%           each column is a tac.
%   scant   [ts te], the start and end time for each frame
%   tx      start time point
%
% OUTPUT
%   k       Patlak slope
%   b       Patlak intercept
%   x       Patlak x coordinate
%   y       Patlak y coordinate
%
% PROGRAMER
%   Guobao Wang @ UC Davis
%   May 10, 2013
% 
% Last Updated: 
%   may 23, 2013
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% check input
%
[num_frm, num_pix] = size(ct);
if isvector(cp)
    cp = repmat(cp(:), [1, num_pix]);
end      
if isempty(tx) | nargin < 4 
   tx = 1;
end
if tx>num_frm
    error('incorrect start time point')
end
if nargin<5
    alg = 0;
end

% prepare data for the Patlak plot
scant = scant / 60;
dt = scant(:,2) - scant(:,1);
sp = cumsum(diag(dt)*cp, 1);
delta = mean(cp(:))*1e-9;
if delta<=0
    delta = 1e-9;
end
cp(cp==0) = delta;
a1 = sp;
if alg
    a2 = cp;
else
    a2 = ones(size(a1));
end
yt = ct; 

% estimate the slope and intercept
tt = tx:num_frm;    % effective time frames
[k, b] = leastsquare2x2(yt(tt,:), a1(tt,:), a2(tt,:));

% output the x and y coordinates of the plot
x = a1;
y = x .* repmat(k, [num_frm, 1]) + repmat(b, [num_frm, 1]);


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

    

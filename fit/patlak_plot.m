function [Ki, b, out] = patlak_plot(cp, ct, scant, tx, atype)
%--------------------------------------------------------------------------
% patlak_plot implements the Patlak graphical plot for estimating the 
% influx rate parameter Ki and intercept b for an irreversible tracer.
%
% INPUT
%   cp      input function
%   ct      time acitivity curve. It can be a vector or a matrix of which
%           each column is a tac.
%   scant   [ts te], the start and end time for each frame. Unit: seconds
%   tx      start time point
%   atype   type of algorithm
%
% OUTPUT
%   Ki      Patlak slope
%   b       Patlak intercept
%   out     additional output data
%
% PROGRAMER: Guobao Wang @ UC Davis, Sept 10, 2008
% Last Updated: 08-26-2013
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% check input
%
[num_frm, num_pix] = size(ct);    
if nargin < 4 
   tx = [];
end
if tx>num_frm
    error('incorrect start time point')
end
if nargin<5 | isempty(atype)
    atype = 0;
end

% scan time
scant = scant / 60;
dt = scant(:,2) - scant(:,1);
t  = mean(scant,2);

% start time
if isempty(tx)
    for m = 1:size(scant,1)
        if scant(m,1)>=30
            tx = m;
            break;
        end
    end
end

tt = tx:num_frm; % effective time frames

% prepare data for the Patlak plot
sp = cumsum(diag(dt)*cp, 1);
delta = mean(cp(:))*1e-2;
if delta<=0
    delta = 1e-9;
end
cp(cp<delta) = delta;


% algorithms
switch atype
    
    case 0

        % coordinate x and y
        xt = repmat(sp./cp,[1 size(ct,2)]);
        yt = ct ./ repmat(cp,[1 size(ct,2)]); 

        % estimate the slope and intercept
        xt_mean = repmat(mean(xt(tt,:),1),[length(tt), 1]);
        yt_mean = repmat(mean(yt(tt,:),1),[length(tt), 1]);
        Ki = sum((xt(tt,:)-xt_mean).*(yt(tt,:)-yt_mean),1)./sum((xt(tt,:)-xt_mean).^2);
        b  = mean(yt(tt,:),1) - Ki.*mean(xt(tt,:),1);

        % output
        if nargout>2
            out.xt = xt;
            out.yt = yt;
            out.yf = xt .* repmat(Ki(:)', [num_frm 1]) + repmat(b(:)', [num_frm 1]);
        end

    case 1
        
        a1 = repmat(sp,[1 size(ct,2)]);
        a2 = repmat(cp,[1 size(ct,2)]);

        % estimate the slope and intercept
        [Ki, b] = lsq2x2(ct(tt,:), a1(tt,:), a2(tt,:));

        % output
        if nargout>2
            out.At = [sp, cp];
            out.ct = ct; % original TACs
            out.cf = a1 .* repmat(Ki(:)',[num_frm 1]) ...
                + a2 .* repmat(b(:)',[num_frm 1]); % fitted TACs
            out.xt = repmat(sp./cp,[1 size(ct,2)]);
            out.yt = ct ./ repmat(cp,[1 size(ct,2)]); 
            out.yf = out.xt .* repmat(Ki(:)', [num_frm 1]) + repmat(b(:)', [num_frm 1]);
        end

end

%--------------------------------------------------------------------------
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

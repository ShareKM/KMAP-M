function [k, x, y, yt, A] = patlak_plot(cp, ct, scant, tx, dk, ptype, Cwb)
%--------------------------------------------------------------------------
% patlakplot implements Patlak plot for estimating influx rate parameter of
% irreversible tracer.
%
% INPUT
%   cp      input function
%   ct      time acitivity curve. It can be a vector or a matrix of which
%           each column is a tac.
%   scant   [ts te], the start and end time for each frame
%   tx      start time point
%   dk      tracer decay constant
%   ptype   type of Patlak model
%
% OUTPUT
%   k       [Patlak slope, Patlak intercept, blood volume]
%   x       Patlak x coordinate
%   y       Patlak y coordinate
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
if nargin<5 | isempty(dk)
    dk = 0;
end
if nargin<6 | isempty(ptype)
    ptype = 1;
end
if nargin<7 | isempty(Cwb)
    Cwb = cp;
end

% scan time and decay correction
scant = scant / 60;
dt = scant(:,2) - scant(:,1);
t  = mean(scant,2);
cc = decaycoef(scant, dk);
cp = diag(cc)*cp;
Cwb = diag(cc)*Cwb;

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

switch ptype
    
    case 1
        
        a1 = repmat(sp./cp,[1 size(ct,2)]);
        a2 = ones(size(a1));
        ct = diag(cc)*ct;
        yt = ct ./ repmat(cp,[1 size(ct,2)]); 

        % estimate the slope and intercept
        [Ki, b] = lsq2x2(yt(tt,:), a1(tt,:), a2(tt,:));

        % output the x and y coordinates of the plot
        k = [Ki, b];
        if nargout>1
            x = a1;
            y = x .* repmat(Ki(:)', [num_frm 1]) + repmat(b(:)', [num_frm 1]);
        end
        if nargout>4
            A = [a1(tt) a2(tt)];
        end
        
    case 2
        
        A(:,1) = sp./cc;
        A(:,2) = cp./cc;
        A(:,3) = Cwb./cc; % whole blood
        A = A(tt,:);
        B = pinv(A'*A)*A';
        k = B * ct(tt,:);
        x = t(tt); y = A*k;  yt = ct(tt,:);
        k = k';

end

% plasma to whole blood ratio
function pb = pbr(t)

pb = 0.386*exp(-0.191*t) + 1.165;


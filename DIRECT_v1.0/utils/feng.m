function cp = feng(t, p, dk, tdelay)
%--------------------------------------------------------------------------
% generates the blood input time activity curve using Feng's model
%

% check input
if nargin<2 | isempty(p)
    p = [851.1 20.8 21.9 4.1 0.01,0.12];
end
if nargin<3 | isempty(dk)
    dk = 0;
end
if nargin<4 | isempty(tdelay)
    tdelay = 0;
end
    
% from seconds to minutes
t = max(0, t-tdelay) / 60;

% time activity of blood
cp = ( p(1)*t + p(2) + p(3) ) .* exp(-p(4)*t) + ...
     p(2) * exp(-p(5)*t) + ...
     p(3) * exp(-p(6)*t);

% decayed 
cp = cp.*exp(-dk*t);

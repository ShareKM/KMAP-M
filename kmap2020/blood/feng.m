function cp = feng(t, p)
%--------------------------------------------------------------------------
% generates the blood input time activity curve using Feng's model
%

% check input
if nargin<2 | isempty(p)
    p = [851.1 20.8 21.9 4.1 0.01,0.12];
end
    
% from seconds to minutes
t = t / 60;
    
% time activity of blood
cp = ( p(1)*t + p(2) + p(3) ) .* exp(-p(4)*t) + ...
     p(2) * exp(-p(5)*t) + ...
     p(3) * exp(-p(6)*t);

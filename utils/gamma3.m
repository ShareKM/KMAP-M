function cp = gamma3(t, p)
% Generate blood input function using the gamma function

% check input
if nargin<2 | isempty(p)
    p = [0.055 10.0 0.7 6 1 1 1.8 22 180];
end

% time activity of blood
cp = p(1) * t.^p(4) .* exp( -t / p(7) ) + ...
     p(2) * t.^p(5) .* exp( -t / p(8) ) + ...
   	 p(3) * t.^p(6) .* exp( -t / p(9) );
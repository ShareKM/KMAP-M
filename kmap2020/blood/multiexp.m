function cp = multiexp(t, p)
% Use four exponential functions to generate the input function
% 

% check input
if isempty(p) | nargin<2
    p = [650 146 105 21 6.7 0.25 0.03 0.0001];
end

% from seconds to minutes
t = t / 60;

% time activity of blood
cp = p(1) * exp( -p(5) * t ) + ...
     p(2) * exp( -p(6) * t ) + ...
     p(3) * exp( -p(7) * t ) + ...
     p(4) * exp( -p(8) * t );


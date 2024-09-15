function [p, u] = gamma3_fit(cp, t, p0, tdelay)
% Fit a blood input function using the gamma function
%
% gbwagn@ucdavis.edu Aug-24-2017
%

% check input
if nargin<3 | isempty(p)
    p = [0.055 10.0 0.7 6 1 1 1.8 22 180];
end

t = max(0,t - tdelay);

% Gamma3 function
F = @(p,t) p(1) * t.^p(4) .* exp( -t / p(7) ) + ...
     p(2) * t.^p(5) .* exp( -t / p(8) ) + ...
   	 p(3) * t.^p(6) .* exp( -t / p(9) );
 
LB = zeros(size(p0)); 
p  = lsqcurvefit(F,p0,t,cp,LB);
u  = F(p,t);
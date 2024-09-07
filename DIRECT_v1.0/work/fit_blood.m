function Cp = fit_blood(t,u,t0, flag)
% fit noisy blood input function
%
% gbwang@ucdavis.edu, 02/11/2015
%
if nargin<3
    t0 = 0;
end
if nargin<4
    flag = 0;
end

[maxu, I] = max(u);
Cp = u;

% 1-st phase
if I>3
    tt = abs( t(1:I-1) - t(I-1) );
    uu = u(1:I-1);

    F  = @(p, tt) (p(1)*tt+p(2)).*exp(-p(3)*tt);
    
    p0 = [100 10 1];
    LB = zeros(size(p0));
    p  = lsqcurvefit(F,p0,tt,uu,LB);

    Cp(1:I-1) = F(p,tt);

end


% 2-nd phase
ii = I+2:length(t);
tt = t(ii) - t0;
mu = max(u(ii));
uu = u(ii)/mu;

F  = @(p, tt) p(1)*exp(-p(2)*tt) + p(3)*exp(-p(4)*tt) + p(5)*tt + p(6);
p0 = [1 .1 0 0 0 0];
LB = zeros(size(p0)); 
%UB = ones(size(p0));
p  = lsqcurvefit(F,p0,tt,uu,LB);

Cp(ii) = F(p,tt)*mu;

if flag
    figure, plot(t,u,'ro--',t,Cp,'b-','LineWidth',2);
end

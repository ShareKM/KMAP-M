function [C_end, Cp] = pred_blood(t,u,tend)
% predict blood input function for a late time point
%
% gbwang@ucdavis.edu, 02/11/2015
%

[maxu, I] = max(u);
Cp = u;

ii = I+3:length(t);
tt = t(ii);
mu = max(u(ii));
uu = u(ii)/mu;

F  = @(p, tt) p(1)*exp(-p(2)*tt) + p(3)*exp(-p(4)*tt) + p(5)*tt + p(6);
p0 = [1 .1 0 0 0 0];
LB = zeros(size(p0)); 
%UB = ones(size(p0));
p  = lsqcurvefit(F,p0,tt,uu,LB);

Cp(ii) = F(p,tt)*mu;
C_end = F(p,tend)*mu;

% figure, plot(t,u,'ro',t(ii),Cp(ii),'b-','LineWidth',2);
% hold on; plot(tend,C_end,'b+-','LineWidth',2)


function [v, b, res] = loganest(cr, scant, dk, xt, T, maxit)
%--------------------------------------------------------------------------
% function [v, b] = loganest(cr, scant, dk, xt, T)
% loganest implements the Logan plot using the transformed matrix form for 
% estimating the volume of distribution of reversible tracers.
%
% INPUT
%   cr      input function
%   scant   [ts te], the start and end time for each frame, (in min)
%   dk      decay coefficient
%   xt      dynamic images. each column is a tac (not time-averaged).
%   T       start time point
%   maxit   maximum iteration number
% OUTPUT
%   v       Logan slope
%   b       Logan intercept
%
% PROGRAMER
%   Guobao Wang @ UC Davis
%   Apr 30, 2012
% 
% Last Updated: 
%   04-30-2012
%
%--------------------------------------------------------------------------

if nargin<5 | isempty(T)
    T = 1;
end
if nargin<6 | isempty(maxit)
    maxit = 100;
end

% setting for Logan estimation
[At, Bt, Sr] = loganset(cr, scant, dk, T);
xt    = [sum(xt(:,1:T),2) xt(:,T+1:end)];

% enforce v>=0
[v, b, res] = lsq2x2nonneg(Bt*xt', repmat(Sr,[1 size(xt,1)]), At*xt', [], [], [], maxit); 


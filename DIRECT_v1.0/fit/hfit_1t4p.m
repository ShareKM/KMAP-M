function [k, cfit] = hfit_1t4p(ct, cp, scant, k0, opt, wt, cwb, wg, cg)
%--------------------------------------------------------------------------
% Fit time activity curves using the one-tissue compartmental model with 
% four parameters. Hybrid Poisson noise model is used, see the reference 
% below for more details:
%
%   Wang GB, Qi J, An optimization transfer algorithm for nonlinear 
%   parametric image reconstruction from dynamic PET data, IEEE Trans. on
%   Medical Imaging, 31(10): 1977-1988, 2012.
%
% INPUT
%   ct      size [nm np], time activity curves (TACs). nm is the total 
%           number of time frames and np is the total number of pixels or 
%           ROIs.
%   cp      size [nm 1], plasma input function. 
%   scant   size [nm 2], start and end time of frames in seconds.
%   k0      size can be [nk 1] or [nk np], initial esitmate of kinetic 
%           parameters. For 1T4P model, nk=4: [fa, fb, K1 k2]'. 
%   opt     option set for kinetic model, see setKopt.m
%   wt      size [nm 1], weighting factors for fitting
%   cwb     size [nm 1], whole blood function
%   wg      weighting factors in the weighted least squares data term
%   cg      TACs in the weighted least squares data term
%
% OUPUT
%   k       size [nk np], estimate of kinetic parameters
%   cfit    size [nm np], fitted time activity curves.
%
%--------------------------------------------------------------------------
% Guobao Wang @ UC Davis, July 8, 2014
% Last Modified: July 8, 2014
%--------------------------------------------------------------------------

% check nan
if any(isnan(ct(:))) | any(isnan(k0(:))) | any(isnan(cp(:)))
    error('NaN value in input');
end

% check size mismatch 
if size(ct,1)~=size(scant,1)
    error('Mismatched frame size.');
end

% make sure the unit of scan time is second
if size(scant,2)~=2
    error('Incorrect scan time input.')
end
if max(scant(:))<180
    error('The scan time requires to be in seconds!');
end

% initial kinetic parameters  [fa fb K1 k2]'
if isempty(k0)
    k0 = [0.01 0.01 0.01 0.01];
end
if isvector(k0)
    k0 = repmat(k0(:), [1,size(ct,2)]);
end
kinit = zeros(4,size(k0,2));
if size(k0,1)==4
    kinit(1:size(k0,1),:) = k0;
else
    error('Unmatched size of kinetic parameter input.')
end

% the option paramters
if nargin<5
    opt = setKopt;
end
if isempty(opt.Decay)
    dk = 0;
else
    dk = opt.Decay;
end
if isempty(opt.LowerBound)
    lb = [0, 0, 0, 0];
elseif length(opt.LowerBound)==4
    lb = opt.LowerBound;
else
    error('Size of lower bounds is 4');
end
if isempty(opt.UpperBound)
    ub = [1.0, 2.0, 2.0, 1.0];
elseif length(opt.UpperBound)==4
    ub = opt.UpperBound;
else
    error('Size of upper bounds is 4');
end
if isempty(opt.PrmSens)
    ps = [1, 1, 1, 1];
elseif length(opt.PrmSens)==4
    ps = opt.PrmSens;
else
    error('Size of active label is 4');
end
if isempty(opt.MaxIter)
    maxit = 20;
else
    maxit = opt.MaxIter;
end
if isempty(opt.TimeStep)
    td = 1.0;
else
    td = opt.TimeStep;
end

% frame weights
if nargin<6 | isempty(wt)
    wt = diag(scant(:,2)-scant(:,1))*ones(size(ct,1),1)/60;
elseif isvector(wt)
    wt = repmat(wt(:),[1,size(ct,2)]);
end

% whole blood volume
if nargin<7 | isempty(cwb)
    cwb = cp;
end

% input function
if length(cp)<scant(end,2)
    cp = finesample(scant, cp, td);
end
if length(cwb)<scant(end,2)
    cwb = finesample(scant, cwb, td);
end

% regularization term
if nargin<8 | isempty(wg)
    wg = zeros(size(wt));
end
if size(wg,1)~=size(wt,1)
    error('unmatched sizes');
end
if nargin<9 | isempty(cg)
    cg = zeros(size(ct));
end
if size(cg,1)~=size(ct,1) | size(cg,2)~=size(ct,2)
    error('unmatched sizes');
end

% check nan
if any(isnan(ct(:))) | any(isnan(cp(:))) | any(isnan(cwb(:))) | any(isnan(kinit(:)))
    error('NaN value in input');
end

% fit time activity curves
[k,cfit] = hfit_1t4p_mex_omp(ct, wt, scant, cp, cwb, dk, kinit, lb, ub, ...
                         ps, maxit, td, cg, wg);






   
    
    

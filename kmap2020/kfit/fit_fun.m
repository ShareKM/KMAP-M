function [kfit, cfit] = fit_fun(kin, cp, opt, ct, wt, cs, ws)
%--------------------------------------------------------------------------
%

% check inputs
if isfield(opt,'KinType')
    ktype = opt.KinType;
else
    ktype = '2t5p';
end
if isfield(opt,'NoiseModel')
    noise = opt.NoiseModel;
else
    noise = 'gsn';
end
if isfield(opt, 'Algorithm')
    algor = opt.Algorithm;
else
    algor = 'lma';
end
if nargin<5
    wt = [];
end
if nargin<6
    cs = [];
end
if nargin<7
    ws = [];
end

% command function
fitfun = str2func(sprintf('%s_%s_%s', noise, ktype, algor));
 
ct(isnan(ct(:))) = 0;
 
% implementation
if strcmp(noise,'png')
    [kfit, cfit] = fitfun(kin, cp, opt, ct, wt, cs, ws);
else
    [kfit, cfit] = fitfun(kin, cp, opt, ct, wt);
end
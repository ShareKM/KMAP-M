function [ct, st] = ktac_2t5p(k, cp, scant, opt, cwb)
%--------------------------------------------------------------------------
% Generate time-activity curves (TACs) using the two-tissue compartmental model.
%
% Inputs:
% - k: Kinetic parameters. This can be a vector or matrix containing the
%      parameters [vb, K1, k2, k3, k4]'.
% - cp: Plasma concentration data.
% - scant: Scan time data, expected as a two-column matrix where each row 
%          represents [start_time, end_time] in seconds.
% - opt: Optional settings for the kinetic modeling, including decay constant, 
%        time step, and blood volume correction parameters.
% - cwb: (Optional) Whole blood concentration data. If not provided, it will 
%        be calculated from `cp` and parameters in `opt`.
%
% Outputs:
% - ct: Computed time-activity curve (TAC) for the given kinetic parameters.
% - st: Sensitivity matrix, which provides the derivatives of the TAC with 
%       respect to each kinetic parameter. This is only computed if required.
%
% This function serves as a wrapper for the compiled MEX function `ktac_2t5p_mex`,
% which performs the actual computation of TACs based on the input kinetic parameters.
%
% Guobao Wang @ 12-10-2009
%
%--------------------------------------------------------------------------

% Adapt the input vector to [vb K1 k2 k3 k4]'
prm = zeros(5, size(k, 2));
if size(k, 1) == 3     % One-tissue compartment model
    prm(1:3, :) = k;
elseif size(k, 1) == 5 % Two-tissue compartment model
    prm(1:size(k, 1), :) = k;
else
    error('Unmatched size of kinetic parameter input.');
end

% Validate the scan time input
if size(scant, 2) ~= 2
    error('Incorrect scan time input.');
end

% (Optional) Uncomment this block if you need to ensure the scan time is in seconds
% if max(scant(:)) < 180
%     error('The scan time requires to be in seconds!');
% end

% Handle optional parameters
if nargin < 4
    opt = setkopt;
end
if isempty(opt.Decay)
    dk = 0;
else
    dk = opt.Decay;
end
if isempty(opt.TimeStep)
    td = 1.0;
else
    td = opt.TimeStep;
end

% Calculate whole blood volume if not provided
if nargin < 5 || isempty(cwb)
    if length(opt.PbrParam) == 3
        pbrp = opt.PbrParam;
    else
        pbrp = [1 0 0];
    end
    if length(cp) == size(scant, 1)
        t = mean(scant, 2);
    elseif length(cp) == scant(end, 2)
        t = (1:length(cp))';
    end    
    pbr = pbrp(1) * exp(-pbrp(2) * t / 60) + pbrp(3);
    cwb = cp ./ pbr;
end

% Resample input function if necessary
if length(cp) < scant(end, 2) / td
    cp = finesample(scant, cp, td);
end
if length(cwb) < scant(end, 2) / td
    cwb = finesample(scant, cwb, td);
end

% Check for NaN values in the inputs
if any(isnan(prm(:))) || any(isnan(cp(:))) || any(isnan(cwb)) 
    error('NaN value in input');
end

% Generate TACs and sensitivity matrix using the compiled MEX function
if nargout == 2 && size(prm, 2) == 1
    [ct, st] = ktac_2t5p_mex(prm, scant, cp, cwb, dk, td);
else
    ct = ktac_2t5p_mex(prm, scant, cp, cwb, dk, td);
    st = [];
end

end

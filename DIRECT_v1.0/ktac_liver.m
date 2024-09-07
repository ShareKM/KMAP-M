function [ct, st] = ktac_liver(k, cp, scant, opt, cwb)
%--------------------------------------------------------------------------
% Generate time activity curves (TACs) from kinetic parameters based on the
% two-tissue compartmental model for liver data.
%
% Inputs:
%   - k: Kinetic parameters (e.g., [vb, K1, k2, k3, k4, ka, fa]')
%   - cp: Plasma input function (concentration over time)
%   - scant: Scan time matrix [start_time, end_time] (in seconds)
%   - opt: Options structure (optional) containing additional settings
%   - cwb: Whole blood concentration (optional)
%
% Outputs:
%   - ct: Generated time activity curve
%   - st: Sensitivity matrix (optional)
%
% This function wraps the compiled MEX file 'ktac_liver_mex', which performs
% the actual computation. The function prepares the input parameters,
% handles optional arguments, and checks for input validity before
% calling the MEX file.
%
% Guobao Wang @ 12-10-2016
%--------------------------------------------------------------------------

% Adapt the input vector 'k' to a standard form of [vb, K1, k2, k3, k4, ka, fa]'
prm = zeros(7, size(k, 2));

% Determine the kinetic model based on the number of parameters
if size(k, 1) == 3  % One-tissue compartment model
    prm(1:3, :) = k;
elseif size(k, 1) == 5  % Two-tissue compartment model
    prm(1:5, :) = k;
elseif size(k, 1) == 7  % Two-tissue compartment model with additional parameters
    prm(1:7, :) = k;
else
    error('Unmatched size of kinetic parameter input.');
end

% Ensure that the scan time is provided in seconds
if size(scant, 2) ~= 2
    error('Incorrect scan time input. Scan time should be a [start_time, end_time] matrix.');
end

% Optional check for scan time being less than 180 seconds (commented out)
if max(scant(:)) < 180
    % error('The scan time requires to be in seconds!');
end

% Handle optional parameters and set defaults if not provided
if nargin < 4  % If options structure 'opt' is not provided
    opt = setkopt;  % Use default options
end

% Set decay constant (dk) from options or default to 0
if isempty(opt.Decay)
    dk = 0;
else
    dk = opt.Decay;
end

% Set time step (td) from options or default to 1 second
if isempty(opt.TimeStep)
    td = 1.0;
else
    td = opt.TimeStep;
end

% Calculate whole blood volume (cwb) if not provided
if nargin < 5 || isempty(cwb)
    % Set default parameters for blood partition ratio
    if length(opt.PbrParam) == 3
        pbrp = opt.PbrParam;
    else
        pbrp = [1 0 0];
    end
    
    % Generate time vector 't' based on the length of cp or scant
    if length(cp) == size(scant, 1)
        t = mean(scant, 2);  % Use the mean of scant as time points
    elseif length(cp) == scant(end, 2)
        t = (1:length(cp))';  % Use indices as time points
    end
    
    % Calculate whole blood concentration (cwb) from plasma input (cp) and pbrp
    pbr = pbrp(1) * exp(-pbrp(2) * t / 60) + pbrp(3);
    cwb = cp ./ pbr;
end

% Resample input functions (cp and cwb) if needed to match the time step
if length(cp) < scant(end, 2) / td
    cp = finesample(scant, cp, td);
end
if length(cwb) < scant(end, 2) / td
    cwb = finesample(scant, cwb, td);
end

% Check for any NaN values in inputs and throw an error if found
if any(isnan(prm(:))) || any(isnan(cp(:))) || any(isnan(cwb))
    error('NaN value found in input parameters.');
end

% Generate time activity curves and sensitivity matrix (if requested)
if nargout == 2 && size(prm, 2) == 1
    % Call the MEX function and return both ct and st
    [ct, st] = ktac_liver_mex(prm, scant, cp, cwb, dk, td);
else
    % Call the MEX function and return only ct
    ct = ktac_liver_mex(prm, scant, cp, cwb, dk, td);
    st = [];  % Sensitivity matrix not requested
end

end

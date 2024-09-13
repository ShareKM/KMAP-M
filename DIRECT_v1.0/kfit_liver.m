function [k, cfit, ctis] = kfit_liver(ct, cp, scant, k0, opt, wt, cwb)
%--------------------------------------------------------------------------
% Fit time-activity curves (TACs) using the two-tissue compartmental model 
% with seven parameters, including dispersion and fractional arterial blood (fa).
% The Gaussian noise model is assumed for the fitting process.
%
% Inputs:
%   - ct: Measured time-activity curve data.
%   - cp: Plasma input function.
%   - scant: Scan time data, in seconds.
%   - k0: Initial kinetic parameters [va, K1, k2, k3, k4, ka, fa, time_delay].
%   - opt: Options structure for fitting, including bounds and settings.
%   - wt: Frame weights (optional).
%   - cwb: Whole blood data (optional).
%
% Outputs:
%   - k: Fitted kinetic parameters.
%   - cfit: Fitted TAC curve.
%   - ctis: Tissue response function (optional).
%
% This function serves as a wrapper around the compiled MEX file 
% 'kfit_liver_mex_omp.mexa64', which performs the actual optimization using
% the Levenberg-Marquardt algorithm in a multi-threaded environment.
%
% Guobao Wang @ 6-3-2016
%--------------------------------------------------------------------------

% Check for NaN values in input data
if any(isnan(ct(:))) || any(isnan(k0(:))) || any(isnan(cp(:)))
    error('NaN value in input');
end

% Check for size mismatch between the time-activity curve (TAC) and scan times
if size(ct, 1) ~= size(scant, 1)
    error('Mismatched frame size.');
end

% Ensure the scan time input has two columns (start and end times in seconds)
if size(scant, 2) ~= 2
    error('Incorrect scan time input.');
end

% Initial kinetic parameters [va K1 k2 k3 k4 ka fa time_delay]
if isempty(k0)
    k0 = [0.05, 0.1, 0.1, 0.1, 0.01, 1, 0, 0]; % Default initial parameters
end
if isvector(k0)
    k0 = repmat(k0(:), [1, size(ct, 2)]); % Replicate initial parameters for all curves
end
kinit = zeros(8, size(k0, 2));
if size(k0, 1) == 8
    kinit(1:size(k0, 1), :) = k0; % Initialize kinit with the provided k0
else
    error('Unmatched size of kinetic parameter input.');
end

% Handle optional parameters
if nargin < 5
    opt = setkopt; % Set default options if not provided
end

% Set decay constant (dk)
if isempty(opt.Decay)
    dk = 0;
else
    dk = opt.Decay;
end

% Set lower bounds for the parameters
if isempty(opt.LowerBound)
    lb = [0, 0, 0, 0, 0, 0, 0, 0]; % Default lower bounds
elseif length(opt.LowerBound) == 8
    lb = opt.LowerBound;
else
    error('Size of lower bounds should be 8.');
end

% Set upper bounds for the parameters
if isempty(opt.UpperBound)
    ub = [1.0, 10, 5.0, 1.0, 0.1, 100, 1, 20]; % Default upper bounds
elseif length(opt.UpperBound) == 8
    ub = opt.UpperBound;
else
    error('Size of upper bounds should be 8.');
end

% Set parameter sensitivity (which parameters to fit)
if isempty(opt.PrmSens)
    ps = [1, 1, 1, 1, 1, 1, 1, 1]; % Default to fitting all parameters
elseif length(opt.PrmSens) == 8
    ps = opt.PrmSens;
else
    error('Size of active label should be 8.');
end

% Set the maximum number of iterations for the optimization
if isempty(opt.MaxIter)
    maxit = 20; % Default maximum iterations
else
    maxit = opt.MaxIter;
end

% Set the time step for the scan duration
if isempty(opt.TimeStep)
    td = 1.0; % Default time step
else
    td = opt.TimeStep;
end

% Handle frame weights
if nargin < 6 || isempty(wt)
    wt = ones(size(ct, 1), 1); % Default weights
elseif isvector(wt)
    wt = repmat(wt(:), [1, size(ct, 2)]); % Replicate weights for all curves
end

% Handle whole blood volume calculation
if nargin < 7 || isempty(cwb)
    if length(opt.PbrParam) == 3
        pbrp = opt.PbrParam; % Get parameters for blood partitioning ratio
    else
        pbrp = [1, 0, 0]; % Default blood partitioning ratio parameters
    end
    if length(cp) == size(scant, 1)
        t = mean(scant, 2); % Use mean scan time if cp matches scant size
    elseif length(cp) == scant(end, 2)
        t = (1:length(cp))'; % Create time vector if cp size doesn't match scant
    end    
    pbr = pbrp(1) * exp(-pbrp(2) * t / 60) + pbrp(3); % Compute blood partitioning ratio
    cwb = cp ./ pbr; % Calculate whole blood concentration
end

% Fine-sample input function to match the time resolution
if length(cp) < scant(end, 2) / td
    cp = finesample(scant, cp, td); % Fine-sample cp
end
if length(cwb) < scant(end, 2) / td
    cwb = finesample(scant, cwb, td); % Fine-sample cwb
end

% Check for NaN values again after processing
if any(isnan(ct(:))) || any(isnan(cp(:))) || any(isnan(cwb(:))) || any(isnan(kinit(:)))
    error('NaN value in input');
end

% Fit the time-activity curves using the compiled MEX file
maxNumCompThreads(20); % Set maximum number of computation threads to 20
[k, cfit] = kfit_liver_mex_omp(ct, wt, scant, cp, cwb, dk, kinit, lb, ub, ...
                         ps, maxit, td);

% Compute impulse response function (IRF) if requested
if nargout > 2
    u0 = zeros(size(cp)); % Create a zero input function
    u0(1) = 1; % Set the first element to 1
    ctis = ktac_liver(k, u0, scant, opt, cwb); % Calculate tissue response
end
end

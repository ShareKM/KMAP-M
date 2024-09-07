function [ct, st] = ktac_1t4p(k, cp, scant, opt, cwb)
%--------------------------------------------------------------------------
% Generate time activity curves from kinetic parameters based on the
% two-tissue compartmental model.
%
% Guobao Wang @ 12-10-2009
%
%--------------------------------------------------------------------------

% adapt the input vector to [va vb K1 k2]'
prm = zeros(4, size(k,2));
if size(k,1)==4     % one-tissue
    prm(1:4,:) = k;
else
    error('Unmatched size of kinetic parameter input.');
end

% make sure the unit of scan time is second
if size(scant,2)~=2
    error('Incorrect scan time input.')
end
if max(scant(:))<180
    error('The scan time requires to be in seconds!');
end

% the option paramters
if nargin<4
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

% whole blood volume
if nargin<5 | isempty(cwb)
    cwb = cp;
end

% input function
if length(cp)<scant(end,2)
    cp = finesample(scant, cp, td);
end
if length(cwb)<scant(end,2)
    cwb = finesample(scant, cwb, td);
end

% check nan
if any(isnan(prm(:))) | any(isnan(cp(:))) | any(isnan(cwb)) 
    error('NaN value in input');
end

% generate time activity curves and sensitivity matrix if needed
if nargout==2 & size(prm,2)==1
    [ct, st] = ktac_1t4p_mex(prm, scant, cp, cwb, dk, td);
else
    ct = ktac_1t4p_mex(prm, scant, cp, cwb, dk, td);
    st = [];
end



   
    
    

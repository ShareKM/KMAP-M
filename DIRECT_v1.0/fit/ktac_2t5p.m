function [ct, st] = ktac_2t5p(k, cp, scant, opt, cwb)
%--------------------------------------------------------------------------
% Generate time activity curves from kinetic parameters based on the
% two-tissue compartmental model.
%
% Guobao Wang @ 12-10-2009
%
%--------------------------------------------------------------------------

% adapt the input vector to [vb K1 k2 k3 k4]'
prm = zeros(5, size(k,2));
if size(k,1)==3     % one-tissue
    prm(1:3,:) = k;
elseif size(k,1)==5 % two-tissue
    prm(1:size(k,1),:) = k;
else
    error('Unmatched size of kinetic parameter input.');
end

% make sure the unit of scan time is second
if size(scant,2)~=2
    error('Incorrect scan time input.')
end
if max(scant(:))<180
    %error('The scan time requires to be in seconds!');
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
    if length(opt.PbrParam)==3
        pbrp = opt.PbrParam;
    else
        pbrp = [1 0 0];
    end
    if length(cp)==size(scant,1)
        t = mean(scant,2);
    elseif length(cp)==scant(end,2)
        t = [1:length(cp)]';
    end    
    pbr = pbrp(1) * exp( -pbrp(2)*t/60) + pbrp(3);
    cwb = cp ./ pbr;
end

% input function
if length(cp)<scant(end,2)/td
    cp = finesample(scant, cp, td);
end
if length(cwb)<scant(end,2)/td
    cwb = finesample(scant, cwb, td);
end

% check nan
if any(isnan(prm(:))) | any(isnan(cp(:))) | any(isnan(cwb)) 
    error('NaN value in input');
end

% generate time activity curves and sensitivity matrix if needed
if nargout==2 & size(prm,2)==1
    [ct, st] = ktac_2t5p_mex(prm, scant, cp, cwb, dk, td);
else
    ct = ktac_2t5p_mex(prm, scant, cp, cwb, dk, td);
    st = [];
end


   
    
    

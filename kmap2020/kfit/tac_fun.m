function [ct, st] = tac_fun(kin, cp, opt)
%--------------------------------------------------------------------------
% Generate time acitivity curves from kinetic parameters
%

% check the option terms
if nargin<3 | isempty(opt)
    setopt_fit;
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
if isempty(opt.ScanTime)
    error('no input for scan time')
else
    scant = opt.ScanTime;
end

% blood input function
% [B, tt] = interp_psf(scant, td, 'linear');
% if length(cp)<=scant(end,2)/td
%     pbr = p2blood_ratio(mean(scant,2), opt.PbrPar);
%     cp  = cp .* exp(-mean(scant,2)/60*dk);
%     Cwb = cp ./ pbr;
%     cin = B * cp;       
%     cwb = B * Cwb;
% else
%     error('incorrect input function');
% end    
tm = mean(scant,2);
cin = interp1(mean(scant,2), cp, 1:scant(end,2), 'linear','extrap');
cwb = cin;

% kinetic model
switch opt.KinType
    case '1t3p'
        tacfun = str2func('ktac_1t3p_mex');
        num_par = 3;
    case '1t3p1'
        tacfun = str2func('ktac_1t3p1_mex');
        num_par = 3;
    case '2t5p'
        tacfun = str2func('ktac_2t5p_mex');
        num_par = 5;    
    case 'srtm'
        tacfun = str2func('ktac_srtm_mex');
        num_par = 4;
        kin = [zeros(1,size(kin,2)); kin];
    otherwise
        error('unknown kinetic model')
end

% check the size of kinetic parameters
if size(kin,1)~=num_par
    error('mismatched size of kinetic parameters');
end

% check NaN
if any(isnan(kin(:))) | any(isnan(cp))
    error('NaN value in input');
end

% generate time activity curves and sensitivity matrix if needed
if nargout==2 & size(kin,2)==1
    [ct, st] = tacfun(kin, scant, cin, cwb, dk, td);
    if strcmp(opt.KinType,'srtm')
        st(:,1) = [];
    end
else
    ct = tacfun(kin, scant, cin, cwb, dk, td);
    st = [];
end   


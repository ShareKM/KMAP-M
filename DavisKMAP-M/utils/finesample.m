function [Ca,tt] = finesample(scant, blood, res, itype)
%--------------------------------------------------------------------------
% resample the blood curve into finer samples to adapt to convolution
% computation.
%

if nargin<3
    res = 1;   % time resolution: 1 sec
end
if nargin<4
    itype = 'linear';
end
ts = scant(:,1) / res;
te = scant(:,2) / res;
tm = ( ts + te ) / 2;
tt = 0.5:1:te(end);
ts(1) = max(1,ts(1));
switch itype
    case 'linear'
        Ca = max(0,interp1(tm,blood,tt,itype,'extrap') ); 
    case 'const'
        for i = 1:size(scant,1)
            Ca(ts(i):te(i)) = blood(i);
        end
    otherwise
        Ca = max(0,interp1(tm,blood,tt,itype,'extrap') ); 
end
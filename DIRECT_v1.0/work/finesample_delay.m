function [Ca_delay,tt] = finesample_delay(scant, blood, res, itype, tdelay, td)
%--------------------------------------------------------------------------
% resample the blood curve into finer samples to adapt to convolution
% computation.
%

if nargin<3 | isempty(res)
    res = 1;   % time resolution: 1 sec
end
if nargin<4 | isempty(itype)
    itype = 'const';
end
if nargin<5
    tdelay = 0;
end
ts = scant(:,1) / res;
te = scant(:,2) / res;
tm = ( ts + te ) / 2;
tt = 0.5:1:te(end);
ts(1) = max(1,ts(1));

% for i = 1:size(scant,1)
%     Ca(ts(i):te(i)) = blood(i);
% end
Ca_delay = max(0,interp1(tm,blood,tt-tdelay/res,'linear','extrap') ); 

% Ca_delay = zeros(size(Ca));
% Ca_delay(tdelay+1:end) = Ca(1:end-tdelay);

Ca_delay = Ca_delay(:);


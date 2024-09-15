function [blood,tt] = fine_sample(scant, blood, res, itype)
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
tt = 1:te(end);
blood = max(0,interp1(tm,blood,tt,itype,'extrap') ); 

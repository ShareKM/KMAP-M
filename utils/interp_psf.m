function [B, tt] = interp_psf(scant, res, itype)
%  Determine the interpolant matrix for interpolating the blood input
%  function with scan grid for numerical computing grid.
%  Guobao Wang @ Mar 20,2010
%

if nargin<2
    res = 1;   % time resolution: 1 sec
end
if nargin<3
    itype = 'linear';
end

% computing grid
ts = floor(scant(:,1) / res);
te = floor(scant(:,2) / res);
tm = ( ts + te ) / 2;
tt = [1:te(end)]';

% total number of scan and computing points
num_frame = size(scant,1);
num_time = length(tt);

% the matrix for time integration
I = zeros(num_frame, num_time);
for m = 1:num_frame
    I(m,[(ts(m)+1):te(m)]) = 1/(te(m)-ts(m));
end

% the response matrix
B = ( inv(I*I') * I )';


    

    

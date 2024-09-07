function [B, cin] = basis_func(k, Cp, scant, dk, td, Tc, E, K1, fv)
%--------------------------------------------------------------------------

if nargin<6
    Tc = 0;
end
if nargin<7
    E = 1;
end
if nargin<8
    K1 = 1;
end
if nargin<9
    fv = 0;
end

% sizes
NumOfBasis = length(k);
NumOfFrame = size(scant,1);

% interpolation
if length(Cp)<scant(end,2)
    [BB,tt] = interp_psf(scant, td, 'linear');
    cin = BB * Cp;
else
    tt = 1:scant(end,2);
    cin = Cp(:);
end

B = zeros(NumOfFrame,NumOfBasis*length(Tc)*length(E));
vec = zeros(length(cin),1);
for i = 1:NumOfBasis
    for n = 1:length(Tc)
        for m = 1:length(E)
            in = sub2ind([length(E) length(Tc) NumOfBasis], m, n, i);
            vec = exp_conv(cin, K1, (k(i)+dk), tt(:), td, Tc(n), E(m));
            vec = vec(:) + fv * cin;
            for j = 1:NumOfFrame
                jj = round([(scant(j,1)+1):scant(j,2)]); % im2single([(scant(j,1)+1):scant(j,2)]);
                dt = scant(j,2)-scant(j,1);
                B(j,in) = sum(vec(jj))/dt;
            end
        end
    end
end

%--------------------------------------------------------------------------
function c = exp_conv(Cp, k1, k2, t, td, Tc, E)
% convolute the exponential function with the input function. 
%--------------------------------------------------------------------------
k1 = k1*td/60;
k2 = k2*td/60;
ek2 = exp(-k2);
c = zeros(1,length(t));
if ek2==1
    c = k1*cumsum(Cp);
else
    if Tc>0
        c(1:Tc) = k1*cumsum(Cp(1:Tc));
    end
    tmp = (1-ek2)/k2;
    c0 = 0;
    for m = 1:length(t)-Tc
        c0 = c0*ek2 + E*k1*tmp*Cp(m);
        C(m+Tc) = c0;
    end
    if Tc>0
        for m = Tc+1:length(t)
            C0 = k1*sum(Cp(m-Tc+1:m));
            c(m) = C0 + C(m);
        end
    else
        c = C;
    end
    
end
c = c(:);

%--------------------------------------------------------------------------
function [B, tt] = interp_psf(scant, res, itype)
%  Determine the interpolant matrix for interpolating the blood input
%  function with scan grid for numerical computing grid.
%  Guobao Wang @ Mar 20,2010
%
%--------------------------------------------------------------------------

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
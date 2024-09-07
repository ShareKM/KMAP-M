function B = basis_func(k, Cp, scant, dk, td)
%--------------------------------------------------------------------------

% sizes
NumOfBasis = length(k);
NumOfFrame = size(scant,1);

% interpolation
[BB,tt] = interp_psf(scant, td, 'linear');
cin = BB * Cp;

B = zeros(NumOfFrame,NumOfBasis);
for i = 1:NumOfBasis
    vec = exp_conv(cin, 1, (k(i)+dk), tt(:), td);
    for j = 1:NumOfFrame
        jj = im2single([(scant(j,1)+1):scant(j,2)]);
        dt = scant(j,2)-scant(j,1);
        B(j,i) = sum(vec(jj))/dt;
    end
end

%--------------------------------------------------------------------------
function c = exp_conv(Cp, k1, k2, t, td)
% convolute the exponential function with the input function. 
%--------------------------------------------------------------------------
k1 = k1*td/60;
k2 = k2*td/60;
ek2 = exp(-k2);
prev = 0;
c = zeros(1,length(t));
if ek2==1
    c = k1*cumsum(Cp);
else
    tmp = (1-ek2)/k2;
    for m = 1:length(t)
        prev = prev*ek2 + k1*tmp*Cp(m);
        c(m) = prev;
    end
end

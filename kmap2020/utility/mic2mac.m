function K = mic2mac(k, ktype)
% transform the microparameters into macroparameters
%

if nargin<2
    ktype = 'ki';
end

switch ktype
    case 'ki'
        tmp = k(2,:)+k(3,:);
        tmp(tmp==0) = mean(tmp)*1e-9;
        K = k(1,:).*k(3,:) ./ tmp;
    case 'vd'
        K = k(1,:)./k(2,:) .* ( 1 + k(3,:)./k(4,:) );
    case 'bp'
        K = k(3,:)./k(4,:);
end
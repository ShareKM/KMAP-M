%--------------------------------------------------------------------------
% calculate the macroparameters from microparameters
%
function [Ki, VD, VD1] = mic2mac(k)
Ki= zeros(1,size(k,2),size(k,3));
idk = (k(2,:)+k(3,:))~=0;
Ki(idk) = k(1,idk).*k(3,idk)./(k(2,idk)+k(3,idk));    % KI
%Ki = Ki(:);

% idx = k(2,:)==0;
% Ki(idx) = 0;

if nargout>1
    VD = zeros(size(Ki));
    VD(idk) = k(1,idk).*k(2,idk)./(k(2,idk)+k(3,idk)+1e-3).^2;
end
if nargout>2
    k234 = k(2,:)+ k(3,:) + k(4,:);
    a1 = 1/2*k234 - 1/2*sqrt(k234.^2-4*k(2,:).*k(4,:));
    a2 = 1/2*k234 + 1/2*sqrt(k234.^2-4*k(2,:).*k(4,:));
    da = a2 - a1;
    VD1 = k(1,:)./da .* ( (k(4,:)-a1)./a1 + (a2-k(4,:))./a2 );
end

    
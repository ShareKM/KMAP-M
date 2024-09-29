function pbr = p2blood_ratio(t, pbr_par)
%--------------------------------------------------------------------------
% calculate the ratio function of plasma-to-blood
%
if nargin<2 | isempty(pbr_par)
    pbr_par = [1 0 0];
end

if length(pbr_par)==1
    pbr = pbr_par * ones(length(t),1);
    
elseif length(pbr_par)==3
    pbr = pbr_par(1) * exp( -pbr_par(2) * t/60 ) + pbr_par(3);
    
elseif length(pbr_par)==length(t)
    pbr = pbr_par(:);
    
else
    error('the input parameters are not correct')
end
    
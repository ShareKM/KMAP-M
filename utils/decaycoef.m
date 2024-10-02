%--------------------------------------------------------------------------    
function cc = decaycoef(scant, dk, c)
%--------------------------------------------------------------------------
% decay correction on time activity curves
%

% the unit of t should be minute

ts = scant(:,1);
dt = scant(:,2)-scant(:,1);
if dk>0
    c1 = exp( - dk * ts );
    c2 = exp( - dk * ( ts + dt ) );
    cc = dk * dt ./ ( c1 - c2 ); 
else
    cc = ones(size(scant,1),1);
end
if nargin>2
    cc = c .* repmat(cc(:),[1 size(c,2)]);
end
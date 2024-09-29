function [At, Bt, Sr] = loganset(cr, scant, dk, T)

if nargin<4 | isempty(T)
    T = 1;
end

% input function
xr    = diag(scant(:,2)-scant(:,1)) * cr(:);

% only using later frames
scant = [scant(1,1) scant(T,2); scant(T+1:end,:)];
xr    = [sum(xr(1:T)); xr(T+1:end)];

% matrices
dt = ( scant(:,2) - scant(:,1) ); 
Bt = zeros(length(dt)-1,length(dt));
for m = 1:size(Bt,1) 
    Bt(m,1:m+1) = 1; 
end
At = [zeros([size(Bt,1) 1]) diag(1./dt(2:end))];

% decay effect
Dk = diag(decaycoef(scant, dk));
Bt = Bt * Dk;
At = At * Dk;
Sr = Bt * xr;

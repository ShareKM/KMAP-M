function [mask,th] = get_mask(x, th, yi, ni, G, Gopt)
%--------------------------------------------------------------------------
% get the object mask using thresholding
%
% Guobao Wang, 09-20-2013
%

% if no image input, reconstruct an image from projection data
if isempty(x)
    
    % check inputs
    imgsiz = Gopt.imgsiz;
    numpix = prod(imgsiz);
    [yi, ri, ni] = sino_preprocess(yi, ri, ni);

    % set Gopt
    Gopt = setGopt(ni, G, Gopt);

    % backprojection
    x = proj_back(G, Gopt, ni.*(yi-ri))./Gopt.sens;
    x(~Gopt.mask) = 0;

end

% hist to find a threshold
if nargin<2 | isempty(th)
    [N,X] = hist(x(:),500);
    Z = [N(1:end-1)' N(2:end)']; 
    idx = find(Z(:,2)>Z(:,1));
    th = X(idx(1));
end

% mask
mask = x>th; 
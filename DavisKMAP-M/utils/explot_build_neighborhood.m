function [N, W] = explot_build_neighborhood(X, mask, nbrtyp, nbrpar, kertyp, kerpar, normflag, midflag, inplaneflag)
%
% Get the index and weight of neighboring voxes in a neighborhood for
% voxels in a mask
%
% gbwang@ucdavis.edu (12-23-2019)
%

%% check input variables
if nargin<1
    X = [];
end
if nargin<2 | isempty(mask)
    error('A mask must be defined');
end
imgsiz = size(mask);
imgdim = length(find(imgsiz>1));
if imgdim==1
    imgsiz = [imgsiz(imgsiz(:)>1) 1 1];
elseif imgdim==2
    imgsiz = [imgsiz(imgsiz(:)>1) 1];
end
if nargin<3 | isempty(nbrtyp)
    nbrtyp = 'cube';
end
if nargin<4 | isempty(nbrpar)
    nbrpar = 3;
end
if length(nbrpar)==1 
    if imgdim==3
        nbrpar = [nbrpar nbrpar nbrpar];
    elseif imgdim==2
        nbrpar = [nbrpar nbrpar 1]; 
    end
end
if nargin<5 | isempty(kertyp)
    kertyp = 'invdist';
end
if nargin<6 | isempty(kerpar)
    kerpar = 1;
end
if nargin<7 | isempty(normflag)
    normflag = 1;
end
if nargin<8 | isempty(midflag)
    midflag = 1;
end
if nargin<9 | isempty(inplaneflag)
    inplaneflag = 0;
end

J = find(mask);
numvox = length(J);

%% neighborhood type
switch nbrtyp
            
    case 'clique'
        
        s   = [1  0  1  1  0  1  0 -1  0  1 -1 -1  1;
               0  1  1 -1  0  0  1  0 -1  1  1 -1 -1;
               0  0  0  0  1  1  1  1  1  1  1  1  1];
        d   = sqrt(sum(abs(s).^2,1));
        
        % image dimension
        if imgdim==1
            idx = 1;
        elseif imgdim==2
            idx = 1:4;
        elseif imgdim==3
            idx = 1:13;
        end
        s = s(:,idx); d = d(idx);
        
        % neighborhood order
        if nbrpar==0.5  % for TV
            idx = d==1;
        elseif nbrpar==1
            idx = d==1;
        elseif nbrpar==2
            idx = d>0;
        end
        s = s(:,idx); d = d(idx);
        
        % neighboring voxel index
        N = zeros(numvox, size(s,2));
        W = zeros(size(N));
        for k = 1:size(s,2)
            N(:,k) = setBoundary3(J, s(:,k), imgsiz);
            W(:,k) = 1/d(k);
        end
        if nbrpar>=1
            for k = 1:size(s,2)
                N(:,k+size(s,2)) = setBoundary3(J, -s(:,k), imgsiz);
                W(:,k+size(s,2)) = 1/d(k);
            end
        end
        
        if not(isempty(X))
            W = W.*calc_wgt(X, J, N, kertyp, kerpar);
        end
        
    case 'cube'
        
        % sizes of neighborhood window
        wlen = 2*floor(nbrpar/2); % length of neighborhood window
        xidx = -wlen(1)/2:wlen(1)/2; 
        yidx = -wlen(2)/2:wlen(2)/2; 
        
        if imgsiz(3)==1 
            zidx = 0; 
        else
            zidx = -wlen(3)/2:wlen(3)/2;            
        end
        
        % image grid
        [I1, I2, I3] = ndgrid(1:imgsiz(1),1:imgsiz(2),1:imgsiz(3));
        I1 = I1(J);
        I2 = I2(J);
        I3 = I3(J);
        
        % weight
        if imgdim==2
            h = fspecial('gaussian', nbrpar(1), nbrpar(1)/(4*sqrt(2*log(2)))); 
        elseif imgdim==3
            h = fspecial3('gaussian', nbrpar); 
        end
        
        % index and distance
        N = zeros(numvox, length(h(:)));
        W = N;
        l = 1; n = 1;
        for x = xidx
            Xnew = setBoundary1(I1 + x, imgsiz(1));
            for y = yidx
                Ynew = setBoundary1(I2 + y, imgsiz(2));
                for z = zidx
                    Znew = setBoundary1(I3 + z, imgsiz(3));
                    if inplaneflag & (abs(x)>0&abs(y)>0&abs(z)>0)
                        disp(sprintf('skip %d',n));
                        n = n + 1;
                        continue;
                    end
                    N(:,l) = Xnew + (Ynew-1).*imgsiz(1) + (Znew-1)*imgsiz(1)*imgsiz(2);
                    W(:,l) = h(n);
                    l = l + 1;
                    n = n + 1;
                end
            end
        end
        
        if l<=size(N,2)
            N(:,l:end) = [];
            W(:,l:end) = [];
        end
        if ~midflag
            i = ceil(size(N,2)/2);
            N(:,i) = [];
            W(:,i) = [];
        end
        if not(isempty(X))
            W = W.*calc_wgt(X, J, N, kertyp, kerpar);
        end
        
    case 'knn'
        k = nbrpar(1);
        [N, D] = knnsearch(X, X, 'dist', 'seuclidean', 'k', k);
        if ~midflag
            N = N(:,2:end);
        end
        W = calc_wgt(X, J, N, kertyp, kerpar);
end

% normalized to 1 and output
if normflag
    W = W./repmat(sum(W,2),[1 size(W,2)]);
end


%% sub functions

%--------------------------------------------------------------------------
function J = setBoundary3(J, d, imgsiz)
%--------------------------------------------------------------------------
[x,y,z] = ind2sub(imgsiz,J);
x = setBoundary1(x+d(1), imgsiz(1)); 
y = setBoundary1(y+d(2), imgsiz(2)); 
z = setBoundary1(z+d(3), imgsiz(3));
J = sub2ind(imgsiz,x,y,z);

%--------------------------------------------------------------------------
function x = setBoundary1(x, N)
%--------------------------------------------------------------------------
if N==1
    x = x(:).^0;
else
    idx = x(:)>N;
    if any(idx)
        x(idx) = N - (x(idx)-N); 
    end
    idx = x(:)<1;
    x(idx) = 1 + (1-x(idx)); 
    x = x(:);
end

%--------------------------------------------------------------------------
function W = calc_wgt(X, J, N, kertyp, kerpar)
%--------------------------------------------------------------------------
D  = zeros(size(N));
for i = 1:size(N,2)
    D(:,i) = sqrt(mean(((X(J,:)-X(N(:,i),:))*diag(1./std(X,0,1))).^2,2));
end

switch kertyp
    case 'invdist'
        W = 1./D;
        
    case 'radial' % radial Gaussian
        W = exp(-D.^2./(2*kerpar.^2));
        
    case 'poly' % polynomial
        for i = 1:size(N,2)
            W(:,i) = ( mean(X(J,:).*X(N(:,i),:),2) + kerpar ).^2;
        end
        
    case 'acos'
        n = kerpar;
        for i = 1:size(N,2)
            xnorm = sqrt(sum(X(J,:).^2,2));
            ynorm = sqrt(sum(X(N(:,i),:).^2,2));
            theta = acos(sum(X(J,:).*X(N(:,i),:),2)./(xnorm.*ynorm));
            theta = real(theta);
            theta(isnan(theta)) = 0;
            if n==0
                Jn = pi - theta;
            elseif n==1
                Jn = sin(theta) + (pi-theta).*cos(theta);
            elseif n==2
                Jn = 3*sin(theta).*cos(theta)+(pi-theta).*(1+2*cos(theta).*cos(theta));
            end
            W(:,i) = 1/pi*xnorm.^n.*ynorm.^n.*Jn;
        end
    otherwise
        error('unknown kernel type')
end 

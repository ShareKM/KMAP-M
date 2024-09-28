function K = explot_build_kernel(U, imgsiz, fold)
% build a kernel matrix from composite image priors (CIP)
%
% gbwang@ucdavis.edu, 12-23-2019
%

if nargin<2 | isempty(imgsiz)
    imgsiz = [256 256 830];
end
if nargin<3 | isempty(fold)
    fold = '..';
end

% Normalization of CIP
U = U * diag(1./std(U,1));

% build an object mask
x = mean(U(:,2:end),2);
[mask,th] = get_mask(x, mean(x)/10);
mask = reshape(mask,imgsiz);
save(sprintf('%s/cip.mask.mat',fold),'mask');

% building the neighborhood and weights
J = find(mask(:));
num_parts = 4;
numk = ceil(length(J)/num_parts);
nmax = 50;
N = [];
W = [];
tic;
for k = 1:num_parts
    k
    mask_k = logical(zeros(imgsiz));
    kk = (1:numk)+(k-1)*numk;
    kk(kk>length(J)) = [];
    mask_k(J(kk)) = 1;
    [N_k, W_k] = explot_build_neighborhood(U, mask_k, 'cube', [9 9 9], 'radial', 1, 1);
    N_k = N_k'; 
    W_k = W_k';
    [temp,I] = sort(W_k,1,'descend');
    I = I + repmat([0:size(W_k,2)-1]*size(W_k,1),[size(W_k,1),1]);
    W_k = W_k(I); 
    N_k = N_k(I);
    W_k = W_k(1:nmax,:)'; 
    N_k = N_k(1:nmax,:)';
    W_k = W_k ./ repmat(sum(W_k,2),[1 size(W_k,2)]);
    N = [N; N_k];
    W = [W; W_k];
    toc;
end

% build the kernel matrix
numvox = prod(imgsiz);
J = find(1-mask(:));
w = ones(length(J),1);
K = sparse(J, J, w, numvox, numvox);

J = find(mask(:));
for n = 1:size(N,2)
    K_n = sparse(J, N(:,n), W(:,n), numvox, numvox);
    K = K + K_n;
end

% make symmetric
K = (K'+K)/2;
K = make_sym(K);
     
disp('Saving the kernel matrix')
save(sprintf('%s/cip.ker.mat',fold),'K','-v7.3');






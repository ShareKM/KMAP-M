function [kfit, cfit] = gsn_srtm_bfm(kin, Cr, opt, ct, wt)
%--------------------------------------------------------------------------
% Estimate kinetic parameters of the simplified reference tissue model
% using the Gaussian likelihood and basis functions algorithm
%
% Guobao Wang @ UC Davis, 10-05-2010
%

% check the option terms
if nargin<3 | isempty(opt)
    setopt_fit;
end
if isempty(opt.Decay)
    dk = 0;
else
    dk = opt.Decay;
end
if isempty(opt.TimeStep)
    td = 1.0;
else
    td = opt.TimeStep;
end
if isempty(opt.ScanTime)
    error('no input for scan time')
else
    scant = opt.ScanTime;
end
if isempty(opt.MaxIter)
    maxit = 20;
else
    maxit = opt.MaxIter;
end

% reference input function
t = mean(scant,2);
pbr = p2blood_ratio(t, opt.PbrPar);
Cr  = Cr .* exp(-t/60*dk);

% Bound constraints
num_par = 3;
if isempty(opt.LowerBound)
    lb = zeros(num_par,1);
elseif length(opt.LowerBound)==num_par
    lb = opt.LowerBound;
else
    error(sprintf('Size of lower bounds is %d', num_par));
end
if isempty(opt.UpperBound)
    ub = ones(num_par,1);
elseif length(opt.UpperBound)==num_par
    ub = opt.UpperBound;
else
    error(sprintf('Size of upper bounds is %d', num_par));
end
if isempty(opt.ParSens)
    ps = ones(num_par,1);
elseif length(opt.ParSens)==num_par
    ps = opt.ParSens;
else
    error(sprintf('Size of active labels is %d', num_par));
end

% frame weights
if nargin<5 | isempty(wt)
    wt = ones(size(ct,1),1);
end
% if isvector(wt)
%     wt = repmat(wt(:),[1,size(ct,2)]);
% end
if size(wt,1)~=size(ct,1)
    error('mismatched sizes of TAC and weighting factors');
end

% total number of TACs and time frames
num_voxel = size(ct,2);
num_frame = size(ct,1);

% get the basis functions
if isfield(opt,'LamVec')
    lam = opt.LamVec;
else
    lam = [0.01, 1.5]; 
end
if isempty(lam) | lam(2)<lam(1)
    error('incorrect lambda vector');
end
if length(lam(:))<=3
    if isempty(lam)
        lam = [0.01 1.0 100];
    end
    if length(lam(:))==2
        lam = [lam(:)' 100];
    end
    lam = logspace(log10(lam(1)), log10(lam(2)), lam(3));
end
B = basis_func(lam, Cr, scant, dk, td);

% least squares estimation
W = diag(sqrt(wt));
R1 = zeros(num_voxel,1);
k2a = zeros(num_voxel,1);
M = zeros(size(B,2)*2,num_frame);
for i = 1:size(B,2)
    A = [Cr B(:,i)];
    [Q,R] = qr(W*A);
    M(2*i-1:2*i,:) = R\Q';
end
rss = zeros(size(B,2),num_voxel);
tht = zeros(2,num_voxel,size(B,2));
for i = 1:size(B,2)
    tht_i = M(2*i-1:2*i,:) * W * ct;
    tht(:,:,i) = tht_i;
    cc = [Cr B(:,i)] * tht_i;
    r = ct - cc;   
    rss(i,:) = sum((W*r).^2, 1);
end
[tmp,i_min] = min(rss);       
for j = 1:num_voxel
    tht_j = tht(:,j,i_min(j));
    tht_new(:,j) = tht_j;
end

% output
R_1 = tht_new(1,:);
k_2a = lam(i_min);
idr = R_1 > 1e-3; 
k_2r = zeros(size(k_2a));
k_2r(idr) = tht_new(2,idr) ./ R_1(idr) + k_2a(idr);
k_2 = tht_new(2,:) + R_1.*lam(i_min);
B_P = R_1.*k_2r./k_2a-1;
kfit = [R_1; k_2; B_P];
cfit = zeros(num_frame,num_voxel);
for m = 1:num_frame
    cfit(m,:) = tht_new(1,:)*Cr(m) + tht_new(2,:).*B(m,i_min);
end


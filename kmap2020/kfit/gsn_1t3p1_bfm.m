function [kfit, cfit] = gsn_1t3p1_bfm(kin, Cr, opt, ct, wt)
%--------------------------------------------------------------------------
% Estimate kinetic parameters of the one-tissue model using the Gaussian
% likelihood and basis functions algorithm
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
if isvector(wt)
    wt = repmat(wt(:),[1,size(ct,2)]);
end
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
[tht, idx, cfit] = lsbfm_2x2(Cr, B, ct, wt);

% output
k2a = lam(idx);
kfit = [tht; k2a];

%--------------------------------------------------------------------------
function [tht, idx, tac_fit] = lsbfm_2x2(cr, B, tac, wt)
%--------------------------------------------------------------------------
% least squares basis function algorithm
%
num_frm = size(B,1);
num_par = size(B,2);
num_pix = size(tac,2);

a11 = repmat((cr.^2)'*wt, [num_par,1]);
a12 = ( repmat(cr,[1,num_par]) .* B )' * wt;
a22 = ( B.^2 )' * wt;
x1  = repmat(cr'*(wt.*tac), [num_par,1]);
x2  = B' * (wt.*tac);

theta1 = zeros(num_par,num_pix);
theta2 = theta1;
det_a = a11.*a22-a12.^2;
idd = det_a~=0;
theta1(idd) = ( a22(idd).*x1(idd) - a12(idd).*x2(idd)) ./ det_a(idd);
theta2(idd) = (-a12(idd).*x1(idd) + a11(idd).*x2(idd)) ./ det_a(idd);

rss = zeros(num_par, num_pix);
for i = 1:num_par
    rss(i,:) = sum(wt.*(tac-[cr B(:,i)]*[theta1(i,:); theta2(i,:)]).^2, 1);
end
[tmp,idx] = min(rss,[],1);

tht = zeros(2,num_pix);
S = zeros(num_frm,num_pix);
for i = 1:num_par
    idi = idx == i;
    if ~isempty(idi)
        tht(:,idi) = [theta1(i,idi);theta2(i,idi)];
        S(:,idi) = repmat(B(:,i),[1,length(find(idi))]);
    end
end
tac_fit = cr * tht(1,:) + S .* repmat(tht(2,:),[num_frm,1]);

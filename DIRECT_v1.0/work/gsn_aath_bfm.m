function [kfit, cfit] = gsn_aath_bfm(Tcmax, Cr, opt, ct, wt)
%--------------------------------------------------------------------------
% Estimate kinetic parameters of the AATH model using the Gaussian
% likelihood and basis functions algorithm
%
% Guobao Wang @ UC Davis, 04-12-2017
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
%pbr = p2blood_ratio(t, opt.PbrParam);
Cr  = Cr .* exp(-t/60*dk);

% Bound constraints
num_par = 4;

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

% blood component 
Sp = basis_func(0, Cr, scant, dk, td);

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

Tc = 0:td:Tcmax;
for i = 1:length(Tc)
    
    % bases
    c = Sp;
    c(scant(:,1)>Tc(i)) = 0;
    c(c==0) = 1e-3;
    B = basis_func(lam, Cr, scant, dk, td, Tc(i));

    % least squares estimation
    [tht, idx, cs, rss] = lsbfm_2x2(c, B, ct, wt);
    k2 = lam(idx);
    ks = [tht; k2];
    
    % output
    if i==1 
        kfit = ks;
        cfit = cs;
        rss_min = rss;
    else
        ridx = rss<rss_min;
        kfit(:,ridx) = ks(:,ridx);
        cfit(:,ridx) = cs(:,ridx);
        rss_min(ridx) = rss(ridx);
    end
    
end



%--------------------------------------------------------------------------
function [tht, idx, tac_fit, rss_min] = lsbfm_2x2(cr, B, tac, wt)
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
[rss_min,idx] = min(rss,[],1);

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

function [kfit, cfit] = gsn_1t4p_bfm(kin, Cr, opt, ct, wt, C_RV)
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
if isfield(opt,'ScanTimeBlood')
    scantBlood = opt.ScanTimeBlood;
else
    scantBlood = scant;
end

% input function
t = mean(scantBlood,2);
pbr = p2blood_ratio(t, opt.PbrPar);
Cr  = Cr .* exp(-t/60*dk);

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

% spill-over effect
if nargin<6 | isempty(C_RV)
    C_RV = zeros(size(ct,1),1);
end

% total number of TACs and time frames
num_voxel = size(ct,2);
num_frame = size(ct,1);

% get the basis functions
if isfield(opt,'LamVec')
    lam = opt.LamVec;
else
    lam = logspace(log10(0.01), log10(1.5), 100); 
end
if isempty(lam) | lam(2)<lam(1)
    error('incorrect lambda vector');
end

if isfield(opt,'TcVal')
    TcVal = opt.TcVal;
else
    TcVal = 0; 
end
if isfield(opt,'ExtFrac')
    E = opt.ExtFrac;
else
    E = 0:0.05:1; 
end
if isfield(opt,'fb')
    fb = opt.fb;
else
    fb = 0; 
end

B = basis_func(lam, Cr, scant, dk, td, TcVal, E);

B0 = basis_func(0, Cr, scant, dk, td, 0, 1, 0, 1);

% least squares estimation
[tht, idx, cfit] = lsbfm_3x3(B0, B, ct, wt, C_RV);

% output
[i,j,k] = ind2sub([length(E) length(TcVal) length(lam)],idx);
k2a = lam(k);
Tc  = TcVal(j);
Ef  = E(i);
kfit = [tht; k2a; Tc; Ef];

%--------------------------------------------------------------------------
function [tht, idx, tac_fit] = lsbfm_3x3(cr, B, tac, wt, cb)
%--------------------------------------------------------------------------
% least squares basis function algorithm
%
num_frm = size(B,1);
num_par = size(B,2);
num_pix = size(tac,2);

a = repmat((cr.*cr)'*wt, [num_par,1]);
b = ( repmat(cr,[1,num_par]) .* B )' * wt;
c = repmat((cr.*cb)'*wt, [num_par,1]);
d = b;
e = ( B.^2 )' * wt;
f = ( repmat(cb,[1,num_par]) .* B )' * wt;
g = c;
h = f;
i = repmat((cb.*cb)'*wt, [num_par,1]);

A = e.*i -f.*h;
BB = -(d.*i-f.*g);
C = d.*h-e.*g;
E = (a.*i-c.*g);
F = -(a.*h-b.*g);
I = a.*e-b.*d;
D = BB;
G = C;
H = F;

det_a = a.*A + b.*BB +c.*C;

x1  = repmat(cr'*(wt.*tac), [num_par,1]);
x2  = B' * (wt.*tac);
x3  = repmat(cb'*(wt.*tac), [num_par,1]);

theta1 = zeros(num_par,num_pix);
theta2 = theta1;
theta3 = theta1;
idd = det_a~=0;
theta1(idd) = ( A(idd).*x1(idd) + BB(idd).*x2(idd) + C(idd).*x3(idd)) ./ det_a(idd);
theta2(idd) = ( D(idd).*x1(idd) + E(idd).*x2(idd) + F(idd).*x3(idd)) ./ det_a(idd);
theta3(idd) = ( G(idd).*x1(idd) + H(idd).*x2(idd) + I(idd).*x3(idd)) ./ det_a(idd);

theta1(theta1<0) = 0;

rss = zeros(num_par, num_pix);
for i = 1:num_par
    rss(i,:) = sum(wt.*(tac-[cr B(:,i) cb]*[theta1(i,:); theta2(i,:); theta3(i,:)]).^2, 1);
end
[tmp,idx] = min(rss,[],1);

tht = zeros(3,num_pix);
S = zeros(num_frm,num_pix);
for i = 1:num_par
    idi = idx == i;
    if ~isempty(idi)
        tht(:,idi) = [theta1(i,idi);theta2(i,idi);theta3(i,idi)];
        S(:,idi) = repmat(B(:,i),[1,length(find(idi))]);
    end
end
tac_fit = cr * tht(1,:) + S .* repmat(tht(2,:),[num_frm,1]) + cb * tht(3,:);


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
    
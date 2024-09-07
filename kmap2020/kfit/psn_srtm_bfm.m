function [kfit, cfit] = psn_srtm_bfm(kin, Cr, opt, ct, wt)
%--------------------------------------------------------------------------
% Estimate kinetic parameters of the simplified reference tissue model
% using the Poisson likelihood and basis functions algorithm
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

% total number of TACs and time frames
num_voxel = size(ct,2);
num_frame = size(ct,1);

% get the basis functions
if isfield(opt,'LamVec')
    lam = opt.LamVec;
else
    lam = [0.01 1.5]; 
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
B = basis_func(lam, Cr, scant, dk, td) * diag(lam);

% time integration
dt = ( scant(:,2) - scant(:,1) )/60;
Cr = Cr .* dt;
B  = diag(dt) * B;
ct = diag(dt) * ct;

% initial kinetics and TACs
if isempty(kin)
    R1 = 1.0*ones(1, num_voxel);
    Ra = 1.01*ones(1, num_voxel) - R1;
    bidx = ones(1, num_voxel);
else
    R1 = kin(1,:);
    Ra = R1;
    [tmp, bidx] = histc(kin(3,:), lam-diff([0 lam])/10);
end
cc = Cr*R1 + repmat(Ra,[num_frame,1]).*B(:,bidx);
logB = log(B*diag(1./sum(B,1)));

% iterative loop
for it = 1:maxit

    % first component
    R1 = R1 ./ sum(Cr) .* ( Cr' * (ct ./ cc) );
    
    % second component
    ca = repmat(Ra,[num_frame,1]).*B(:,bidx) ./ cc .* ct;
    [tmp, bidx] = max( logB' * ca, [], 1); 
    Ra = sum(ca,1)./sum(B(:,bidx),1);
    
    % total concentration
    cc = Cr*R1 + repmat(Ra,[num_frame,1]).*B(:,bidx);

end

% output
BP = R1 + Ra - 1;
k2 = lam(bidx) .* (BP+1);
kfit = [R1; k2; BP];
cfit = diag(1./dt) * cc;


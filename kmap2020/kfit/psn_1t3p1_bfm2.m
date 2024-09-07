function [kfit, cfit] = psn_1t3p1_bfm(kin, cp, opt, ct, wt)
%--------------------------------------------------------------------------
% Estimate kinetic parameters of the one-tissue kinetic model using the
% Poisson likelihood and basis functions algorithm
%
% Guobao Wang @ UC Davis, 10-06-2010
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
cp  = cp .* exp(-t/60*dk);
pbr = p2blood_ratio(t, opt.PbrPar);
cwb = cp ./ pbr;

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
B = basis_func(lam, cp, scant, dk, td) * diag(lam);

% time integration
dt = ( scant(:,2) - scant(:,1) )/60;
cwb = cwb .* dt;
B  = diag(dt) * B;
ct = diag(dt) * ct;

% initial kinetics and TACs
if isempty(kin)
    vb = 1.0*ones(1, num_voxel);
    k1 = 1.01*ones(1, num_voxel) - vb;
    bidx = ones(1, num_voxel);
else
    lamm = lam - diff([0 lam])/10;
    lamm(1) = -inf; lamm(end) = inf;
    [tmp, bidx] = histc(kin(3,:), lamm);
    vb = kin(1,:);
    k1 = kin(2,:)./lam(bidx);
end
cc = cwb*vb + repmat(k1,[num_frame,1]).*B(:,bidx);
logB = log(B);

% iterative loop
for it = 1:maxit

    % linear coefficients
    for i = 1:20
        vb = vb ./ sum(cwb) .* ( cwb' * (ct ./ cc) );
        k1 = k1 ./ sum(B(:,bidx)) .* sum(B(:,bidx)./cc.*ct);
        cc = cwb*vb + repmat(k1,[num_frame,1]).*B(:,bidx);
    end
    
    % nonlinear coefficient
    ca = repmat(k1,[num_frame,1]).*B(:,bidx) ./ cc .* ct;
    [tmp, bidx] = max( logB' * ca - sum(B,1)'*k1, [], 1); 
    cc = cwb*vb + repmat(k1,[num_frame,1]).*B(:,bidx);
   
end

% output
kfit = [vb; k1.*lam(bidx); lam(bidx)];
cfit = diag(1./dt) * cc;


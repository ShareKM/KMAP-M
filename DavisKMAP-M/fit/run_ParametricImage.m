function run_ParametricImage(dynImgInp, foldOut, fileBIF, Dopt)
% Generate parametric image from dynamic PET data with model selection
% Inputs:
% - dynImgInp: path to input dynamic images
% - foldOut: path for output parametric images
% - fileBIF: path to input function
% - Dopt: Options structure for parametric image generation
%% check inputs
if nargin<1 | isempty(dynImgInp)
    error('No input dynamic PET found.');
end
if nargin<2 | isempty(foldOut)
    error('Please provide a path to store output parametric images.');
end
if nargin<3 | isempty(fileBIF)
    fileBIF  = foldOut;
end
if nargin<4 | isempty(Dopt)
    Dopt = [];
end
if not(isfield(Dopt,'frameType'))
    Dopt.frameType = '1H29';
end
if not(isfield(Dopt,'kernelFlag'))
    Dopt.kerFlag = 1;
end
if not(isfield(Dopt,'delayFlag'))
    Dopt.delayFlag = 1;
end
if not(isfield(Dopt,'smoothFlag'))
    Dopt.smoothFlag = 1;
end
if not(isfield(Dopt,'scale4SUV'))
    Dopt.scale4SUV = 1;
end
if not(isfield(Dopt,'weightType'))
    Dopt.weightType = 0;
end


%% load the dynamic images
load(dynImgInp);   
% scan duration
dt = imghdr.dt(:);
scant = [cumsum([0; dt(1:end-1)]), cumsum(dt)];
scant(end)=round(scant(end));
% image
imgsiz = imghdr.imgsiz;
numvox = prod(imgsiz);

X = reshape(X, [numvox, length(X(:))/numvox]);

%% Composite Image Prior (CIP) for kernel matrix

switch Dopt.frameType
    case '1H29'
        M = {1:11, 12:20, 21:25, 26:29}; % composite frames
end

maxv = 3e4;
zz = 1:imgsiz(3);
if Dopt.kernelFlag
    for i = 1:length(M)
        if M{i}(end)<=size(X,2)
            disp(sprintf('CIP %d, %d-%g min.', i, scant(M{i}(1),1)/60, scant(M{i}(end),2)/60))
        else
            disp('WARNING: incomplete dynamic frames');
            break;
        end
        
        mm = M{i};
        x = (X(:,mm)*dt(mm))/sum(dt(mm));
        U(:,i) = x;

        % display
        x = reshape(x, imgsiz);
        mip = squeeze(max(x(:,:,zz),[],2));
        mip = flipdim(mip,1);    
        figure
        imagesc(imrotate(mip,-90),[0 maxv]);
        axis off;
        colormap(flipud(gray));
        daspect([4, 4 1]);
    end
end

%% build the kernel matrix

if Dopt.kernelFlag
    K = explot_build_kernel(U, imgsiz, foldOut);
else
    load(sprintf('%s/cip.ker.mat',foldOut));
end

%% smoothing dynamic images

% kernel smoothing
if Dopt.smoothFlag
    X = K * X; 
end



%% load the blood input data

% scan duration
t = mean(scant,2);
load(fileBIF);
Ca  = ascending_aorta(:,1);

%% TACs

% SUV scale for match with the input function
C = X' * Dopt.scale4SUV; 
clear X;

% define a mask
x = reshape(mean(C,1),imgsiz);
x = gaussfilt(x,imgsiz,2);
mask = logical(zeros(imgsiz));
mask(:,:,:) = 1;
mask = mask & x>mean(x(:))*1e-1;
figure,imagesc(flipud(imrotate(squeeze(max(mask,[],2)),90)),[0 1]);
daspect([4 4 1]);

% TACs
C = C(:,mask);

%% Define fit options

% fit option 
opt.Decay      = log(2)/109.8;
opt.TimeStep   = 1.0;
opt.PbrParam   = [1 0 0];
opt.ScanTime   = scant;
opt.LowerBound = [0 0  0  0  0 -10];
opt.UpperBound = [1 5  5  5  0.5 50];
opt.MaxIter    = 100;
opt.PrmSens    = [1 1 1 1 0 1];
opt.Initials   = [0.1 0.1 0.1 0.1 0 0]';

% weighting factor
switch Dopt.weightType
    case 0
        wt = ones(size(C));
    case 1
        wt = diag(dt.*exp(-log(2)/109.8*t/60))*(1./max(C,mean(C(:))*1e-3));
end
wt = wt./repmat(sum(wt),[length(dt) 1]);

% fitting with time delay
mm = 1:length(dt);
% estimate delay?
dflag = Dopt.delayFlag;
if dflag
    opt.PrmSens(6)    = 1;
    n_add = 1;
else
    opt.PrmSens(6)    = 0;
    n_add = 0;
end

% Model 0: 0T0P
Opt{1} = opt;
Opt{1}.OptName  = '0T0P';
Opt{1}.PrmSens(1:5)  = [0,0,0,0,0];
Opt{1}.Initials = [1 0.0 0.0 0.0 0.0 0]';
Opt{1}.MaxIter  = 1;

% Model 1: 1T3P
Opt{2} = opt;
Opt{2}.OptName  = '1T3P';
Opt{2}.PrmSens(1:5)  = [1,1,1,0,0];
Opt{2}.Initials = [0.01 0.01 0.01 0.0 0.0 0]';

% Model 2: 2T4P
Opt{3} = opt;
Opt{3}.OptName  = '2T4P';

%% Parametric Imaging

% all voxels
nn = [1,2,3];
for l = 1:length(nn)
    tic;
    n = nn(l);
    Opt_n = Opt{n};
    [Cd,tt] = finesample(scant(mm,:), Ca(mm), Opt_n.TimeStep);
    Cd=Cd(:);
    [kfit, cfit] = kfit_2tcm(C(mm,:), Cd(1:int32(scant(mm(end),2))), scant(mm,:), Opt_n.Initials, Opt_n, wt(mm,:));
    tdel = kfit(6,:);
    AIC = model_select(C,cfit,length(find(Opt_n.PrmSens))+n_add,wt);
    toc;
    % sae data
    save(sprintf('%s/voxfit_%s',foldOut,Opt_n.OptName),'kfit','tdel','Opt_n','mask','AIC');
    
    % model selecion
    if l==1
        AIC_prev = AIC;
        Ks = [kfit(1:5,:); tdel];
        O = zeros(size(tdel));
    else
        idx = AIC<AIC_prev;
        Ks(:,idx) = [kfit(1:5,idx); tdel(idx)];
        O(idx) = n-1;
        AIC_prev(idx) = AIC(idx);
    end

end

% save the final result
[Ks(7,:),Ks(8,:)] = mic2mac(Ks(2:4,:)); % Ki and apparent VD
save(sprintf('%s/voxfit_%s',foldOut,'oms'),'Ks','mask','O');

%% output parametric images
kk    = [1 2 3 4 7 8];
kname = {'vb', 'K1', 'k2', 'k3', 'Ki', 'V0'};
for i = 1:length(kk)
    k = kk(i);
    x = zeros(imgsiz);
    x(mask) = Ks(k,:);
    
    % smoothing
    if Dopt.smoothFlag
        x = K*x(:); 
    end
    x = gaussfilt(x,imgsiz,1.2);
    x = reshape(x,imgsiz);
    
    % save data
    put_data(sprintf('%s/parimg_%s_%s.raw',foldOut,'oms',kname{i}),x);
end

% model order
x = zeros(imgsiz)-1;
x(mask) = O;
put_data(sprintf('%s/parimg_%s_%s.raw',foldOut,'oms','label'),x);


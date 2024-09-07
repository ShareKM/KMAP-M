function Cp = kfit4d(filek, Cp, fitopt, filex, wt, imgsize, frmidx)
%--------------------------------------------------------------------------
%

% total number of voxel
num_vox = prod(imgsize);

% frame index
[num_frm, frmidx] = frame_index(frmidx);

% fitting setting
scant = fitopt.ScanTime;
dt = ( scant(:,2) - scant(:,1) ) / 60; % in minute
num_par = fitopt.NumOfPar;

% target roi
roi_tar = logical(get_data(fitopt.TarRoi, inf, 'int8'));
roi_tar = reshape(roi_tar, [imgsize(1)*imgsize(2), imgsize(3)]);

% get data
xt = get_data(filex); xt(isnan(xt(:))) = 0;
xt = reshape(xt, [num_vox length(xt)/num_vox]);
xt = reshape(xt(:,frmidx), [imgsize(1)*imgsize(2), imgsize(3), num_frm]);
xt_fit = zeros(size(xt));

% input function
if isempty(Cp)
    roi_bld = logical(get_data(fitopt.BldRoi, inf, 'int8'));
    pbr = p2blood_ratio(mean(scant,2), fitopt.PbrPar);
    Cp = mean(xt(roi_bld(:),:),1)'./dt .* pbr;
    Cp = Cp .* exp(opt.Decay.*mean(scant,2)/60);
end
    
% initial kinetic parameters
k = get_data(filek);
if isempty(k)
    kinit = fitopt.KinInit(:)';
    k = zeros(num_vox, num_par);
    k(roi_tar(:),:) = repmat(kinit, [length(find(roi_tar(:))), 1]);
end
k = reshape(k, [imgsize(1)*imgsize(2), imgsize(3), num_par]);

% slice-wise parameter estimation
for z = 1:imgsize(3)
    z
    ct = diag(1./dt) * squeeze(xt(:,z,:))';
    kk = squeeze(k(:,z,:))';
    
    if length(find(roi_tar(:,z)))>0
        cc = zeros(size(ct));
        roi_z = roi_tar(:,z);
        [kk(:,roi_z), cc(:,roi_z)] = fit_fun(kk(:,roi_z), Cp, fitopt, ct(:,roi_z), wt); 
        k(:,z,:) = reshape(kk',[size(kk,2), 1, size(kk,1)]);
    end
    
end

% write into file
put_data(filek, k);

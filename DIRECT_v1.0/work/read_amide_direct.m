function [ROIName, TAC, out] = read_amide_direct(filename, suvtype)

if nargin<2
    suvtype = 'mean';
end

fid = fopen([filename]);
i = 1;
while ~feof(fid)

    % read ROI name
    if i==1
        iskip = 4;    
    else
        iskip = 1;    
    end
    C = textscan(fid, [repmat('%s ',1,5)], 1, 'HeaderLines', iskip);
    temp = C{3};
    ROIName{i} = temp{1};
    
    
    % read mean time activity
    C = textscan(fid, [repmat('%f ',1,14)], 'HeaderLines', 4);
    switch suvtype
        case 'mean'
            TAC(:,i) = C{7};
        case 'max'
            TAC(:,i) = C{11};
    end
    
    % scan time
    out.dt = C{2};
    
    % ROI volume
    out.vol(i) = C{12}(1);
            
    % next ROI
    i = i + 1;
    
end
fclose(fid);

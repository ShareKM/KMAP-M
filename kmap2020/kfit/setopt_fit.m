function opt = setopt_fit(varargin)
%--------------------------------------------------------------------------
% set the option structure for KMAP
%

% check input
if nargin==0 | isempty(varargin{1})
    opt = struct;
else
    opt = varargin{1};
end
if nargin>1 & rem(nargin,2)==0
    error('incorrect inputs')
end

% all items needed for KMAP
opt_items = {'KinType', 'Decay', 'UpperBound', 'LowerBound', 'ParSens', ...
             'MaxIter', 'TimeStep', 'PbrPar', 'ScanTime', 'NumOfPar', ...
             'BldRoi', 'TarRoi'};

% initialization
for i = 1:length(opt_items)
    if ~isfield(opt, opt_items{i})
        opt = setfield(opt, opt_items{i}, []);
    end
end
if nargin<3
    return;
end

% read and assign values for each item
for i = 2:2:nargin
    if ischar(varargin{i})
        if nargin<(i+1)
            error('no input for this item');
        end
        item_i = varargin{i};
        if strcmp(item_i, 'KinType')
            tmp = varargin{i+1};
            if isstr(tmp)
                opt.KinType = tmp;
            else
                error('incorrect input for kinetic model type')
            end
        if strcmp(item_i, 'Decay')
            dk = varargin{i+1};
            if isscalar(dk)
                opt.Decay = dk;
            else
                error('incorrect input for decay constant')
            end
        elseif strcmp(item_i, 'TimeStep')
            td = varargin{i+1};
            if isscalar(td)
                opt.TimeStep = td;
            else
                error('incorrect input for time step')
            end   
        elseif strcmp(item_i, 'UpperBound')
            ub = varargin{i+1};
            if isvector(ub)
                opt.UpperBound = ub;
            else
                error('incorrect input for upper bounds')
            end
        elseif strcmp(item_i, 'LowerBound')
            lb = varargin{i+1};
            if isvector(lb)
                opt.LowerBound = lb;
            else
                error('incorrect input for lower bounds')
            end
        elseif strcmp(item_i, 'ParSens')
            ps = varargin{i+1};
            if isvector(ps)
                opt.ParSens = ps;
            else
                error('incorrect input for active parameters')
            end
        elseif strcmp(item_i, 'Maxiter')
            it = varargin{i+1};
            if isscalar(it)
                opt.MaxIter = it;
            else
                error('incorrect input for maximum iteration number')
            end    
        elseif strcmp(item_i, 'PbrPar')
            pbrp = varargin{i+1};
            if isvector(pbrp)
                opt.PbrPar = pbrp;
            else
                error('incorrect input for whole blood function')
            end
        elseif strcmp(item_i, 'BldRoi')
            broi = varargin{i+1};
            if isvector(broi)
                opt.BldRoi = broi;
            else
                error('incorrect input for blood region')
            end
        elseif strcmp(item_i, 'TarRoi')
            troi = varargin{i+1};
            if isvector(troi)
                opt.TarRoi = troi;
            else
                error('incorrect input for target region')
            end
        else
            error('unknown item');
        end
        
    else
        error('unknown item');        
    end
end
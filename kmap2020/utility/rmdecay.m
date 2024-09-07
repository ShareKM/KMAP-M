function yout = rmdecay(scant,y,radtype)
%FUNCTION redecay(tstart,dt,y,radiotype) is used for removing the
%          effects of radio decay on PET measurements.
%- INPUT -
%  scant    scan time information. if it is a structure, scant can be one
%           of: scant.ts - start scan time;
%               scant.dt - frame durations
%               scant.t  - time at the middle points of frames; t=ts+dt/2;
%           The user can
%           (1) provide both scant.ts and scant.dt;
%           (2) provide scant.dt; It will have a default scant.ts(1) = 0;  
%           (3) provide scant.t;
%           If scant is given as a vector, t = scant.
%
%  y        PET measurements; If t (ts and dt) is a scalar (single frame),
%           y can be a vector or a matrix; If t is a vector, y can be a 
%           multi-dimensional matrix, whose first dimension corresponds to
%           one TAC vector.
%
%  radtype  the radio type for the tracer used.
%
%- OUTPUT -
%       yout        decay-corrected PET measurements
%
%GB@2006.05


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read and check input variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if nargin < 2, error('too less input arguments'), end

%_________________________________________________________________________
% for the argument - scant
%
ts = [];    % start times
dt = [];    % frame durations
if isstruct(scant)
    if isfield(scant,'ts') && isvector(scant.ts)
        ts = scant.ts; ts = ts(:);
    elseif isfield(scant,'dt') && isvector(scant.dt) 
        dt = scant.dt; dt = dt(:); 
        if isempty(ts), ts = cumsum([0; dt(1:end-1)]), end
    elseif isfield(scant,'t') && isvector(scant.t)
        t = scant.t; t = t(:); 
        if isempty(dt), dt = gradient(t);
        elseif isempty(ts), ts = t - dt/2;
        end
    else error('wrong structure-data of scan time.');
    end
elseif isvector(scant)
    dt = scant(:); 
    ts = cumsum([0; dt(1:end-1)]);
else error('wrong structure-data of scan time.');
end

if length(ts)~=length(dt)
    error('incorrect scan time information.');
end

%_________________________________________________________________________
% for the arguments - y
%
if isvector(y)
    y = y(:);
end
if ( length(dt) > 1 ) && ( length(dt)~=size(y,1) )
    error('sizes of input data dismatch');
end

%_________________________________________________________________________
% for the arguments - radtype
%
if nargin < 3
    radtype = 'F-18';
end
switch lower(radtype)
    case {'f-18','f18'}
        halflife = 109.8; % min
    case {'c-11','c11'}
        halflife = 20.4;
    case {'n-13','n13'}
        halflife = 9.98;
    case {'o-15','o15'}
        halflife = 2.05;
    otherwise
        error('no such a radio type')
end


%%%%%%%%%%%%%%%%%%%%%%
%% decay correction
%%%%%%%%%%%%%%%%%%%%%%
%_________________________________________________________________________
% compute the correction coefficients
%
lamda = log(2) / halflife;
c1 = exp( - lamda * ts );
c2 = exp( - lamda * ( ts + dt ) );
cc = lamda * dt ./ ( c1 - c2 );  
%_________________________________________________________________________
% correct
%
if length(dt) == 1
    yout = y .* cc;
else
    sizey = size(y);
    yout = y .* repmat(cc,[1 sizey(2:end)]);
end

return
        
    
    

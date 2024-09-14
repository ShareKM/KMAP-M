% set up matlab path for DIRECT

% look for directory where this file setup.m is installed
if ~(exist('DIRECTdir', 'var'))
	DIRECTdir = which('setup');
	DIRECTdir = fileparts(DIRECTdir);
end

if DIRECTdir(end) ~= filesep % make sure there is a '/' at end of directory
	DIRECTdir = [DIRECTdir filesep];
end

path([DIRECTdir 'demo'], path);       % demo
path([DIRECTdir 'utils'], path);      % utils
path([DIRECTdir 'work'], path);       % work: my own work folder
path([DIRECTdir 'fit'], path);        % fitting

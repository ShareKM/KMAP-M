% set up matlab path for KMAP-M

% look for directory where this file setup.m is installed
if ~(exist('kmapdir', 'var'))
	kmapdir = which('setup');
	kmapdir = fileparts(kmapdir);
end

if kmapdir(end) ~= filesep % make sure there is a '/' at end of directory
	kmapdir = [kmapdir filesep];
end

path([kmapdir 'demo'], path);       % demo
path([kmapdir 'utils'], path);      % utils
path([kmapdir 'fit'], path);        % fitting
path([kmapdir 'tac'], path);        % TAC
path([kmapdir 'data'], path);       % data
path(genpath([kmapdir 'Precompiled_Binaries']), path);  % binary files
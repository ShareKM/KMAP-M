
% look for directory where this file setup.m is installed
if ~exist('curdir', 'var')
	curdir = which('setup');
	curdir = fileparts(curdir);
end

if curdir(end) ~= filesep % make sure there is a '/' at end of directory
	curdir = [curdir filesep];
end

P = path; path(P,[curdir,'blood']);
P = path; path(P,[curdir,'kfit']);
P = path; path(P,[curdir,'utility']);
P = path; path(P,[curdir,'kmap.c']);

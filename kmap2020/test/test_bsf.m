clear; clc;

% scan time setting
dt = [ones(4,1)*20;  ones(4,1)*40;  ones(4,1)*60; ones(4,1)*180; ones(14,1)*300;];
scant = [cumsum([0; dt(1:end-1)]), cumsum(dt)];
t = mean(scant,2);

opt.Decay = log(2)/20.4;
opt.TimeStep = 1.0;

% input function
cp = feng(t);
cp(1) = 0;

% interpolate
res = 1.0;
itype = 'linear';
[blood,tt] = finesample(scant, cp, res, itype);

% interpolant matrix
B = intpsf(scant, res, itype);
u = B * cp;

% check
figure, plot(tt, blood, ':', t, cp, 'o', tt, u)
figure, plot(t,cp,'o', t, inv(B'*B)*B'*u)
clear; clc;
run('/home/gbwang/work/kmap/setup.m');

% scan time setting
dt = [ones(4,1)*20;  ones(4,1)*40;  ones(4,1)*60; ones(4,1)*180; ones(14,1)*300;];
scant = [cumsum([0; dt(1:end-1)]), cumsum(dt)];
t = mean(scant,2);

% optional setting
opt.ScanTime = scant;
opt.Decay = log(2)/109.8;
opt.TimeStep = 1.0;
opt.PbrPar = [1 0 0];
opt.KinType = 'dbif';

% parameters
k0 = [0.05 1 0.5 0.1 0.01 0.75 1.5]'; % = [fv k1 ki k3 k4 fv Kd]

% input function
cp = feng(t);
cwb = cp;
cp = cp.*exp(-opt.Decay*t/60); % for consistance with COMKAT

[c, s] = tac_fun(k0, cp, opt);

figure, plot(t, c)

% test jac computation
if 0
    for i=1:length(k0)
        kk = k0; kk(i) = kk(i) + 1e-5;
        c1 = ktac_2t5pk(kk, cp, scant, opt, cwb);
        figure, plot([(c1-c)/1e-5, s(:,i)])
    end
end


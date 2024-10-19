% This is a demo file to test the Relative Patlak (RP) plot with a real TAC
clear; clc
run('../setup.m');

%% load data
load('/Users/gbwang/Documents/GitHub/DavisKMAP/DavisKMAP-M/demo/demo_realdata/data/demo_2tcm_realdata.mat');

tac = GM;
scant = [cumsum([0; dt(1:end-1)]), cumsum(dt)];
t = mean(scant,2);

% start at t = 18 min;
tx = 20;

%% fit the TAC using the Relative Patlak plot

% used time points
tt = tx:length(t);

% the TACs for RP plot
Cp_RP = Cp(tt);
Ct_RP = tac(tt);
scant_RP = scant(tt,:);

% Relative Patlak plot
[Ki_RP, b_RP, out] = relative_patlak(Cp_RP, Ct_RP, scant_RP);

%% display the fitting result

% the RP plot
figure,
plot(out.xt, out.yt, 'r-.o', out.xt, out.yf, 'b-', 'LineWidth',2);
xlabel('X-axis of the Relative Patlak Plot');
ylabel('Y-axis of the Relative Patlak Plot');
title('Relative Patlak Plot')
set(gca,'FontSize',14);
legend('Noisy','Fitted','Location','Best')
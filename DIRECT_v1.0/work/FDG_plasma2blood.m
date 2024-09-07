function pbr = FDG_plasma2blood(t, Cwb, method);
% a population-based ratio of plasma to whole blood activity for FDG. The
% data is extracted from 
%  https://jnm.snmjournals.org/content/32/7/1433
%
% gbwang@ucdavise.edu, April 25, 2021
%

% time in seconds
if nargin<1 | isempty(t)
    t = 0:3600;
end

% peak time
if nargin<2 | isempty(Cwb)
    t_peak = 10;
    Cwb = zeros(size(t));
else
    [tmp, idx_peak] = max(Cwb);
    t_peak = t(idx_peak);
end

% method
if nargin<3 | isempty(method)
    method = 'Yale';
end

switch method
        
    case 'Todai' % Ohtake et al. JNucIMed 1991;32:1432-143
        % ratio B
        ratio_B = 1.09;

        % ratio A
        x = [0 15 30 45 60 90 120 180 600 3600];
        y = [1.037 1.126 1.071 1.061 1.042 1.040 1.027 1.015 1 0.97];

        for i = 1:length(t)
            ti = t(i);
            if ti<t_peak
                ratio_A(i) = interp1([0 t_peak], [1 y(1)], ti);
            else
                ratio_A(i) = interp1(x,y,ti-t_peak,'linear','extrap');
            end
        end

        % combined ratio
        pbr = ratio_B * ratio_A;
        pbr = reshape(pbr,size(Cwb));
        
    case 'Yale' % Naganawa et al. EJNMMI Phys . 2020 Nov 23;7(1):67
        bpr = 0.97 - 0.06*exp(-0.085*t/60);
        pbr = 1./bpr;

    case 'Scale'
        pbr = 1.09;
        
    otherwise
        pbr = 1;

end

% display
if nargin<1
    figure,plot(t,pbr,'-');
    xlim([0 3600]);
    ylim([0.8 1.5]);
    set(gca,'FontSize',14);
    xlabel('scan time (s)')
    ylabel('Plasma-to-blood ratio')
    
end
    

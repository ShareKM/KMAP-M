function [scant, blood, analytical] = fitanal(k1, k2)

% starting point of each frame
t = [      0
           2
           4
           6
           8
          10
          12
          14
          20
          26
          32
          52
          92
         172
         352
         532
         712
         892
        1072
        1252
        1432
        1612
        1792
        2152
        2512
        2872
        3232
]/60;  % min

tdif =   [2.0000
    2.0000
    2.0000
    2.0000
    2.0000
    2.0000
    2.0000
    6.0000
    6.0000
    6.0000
   20.0000
   40.0000
   80.0000
  180.0000
  180.0000
  180.0000
  180.0000
  180.0000
  180.0000
  180.0000
  180.0000
  180.0000
  360.0000
  360.0000
  360.0000
  360.0000
  360
]/60; % mins
%t=[0:1:2599]'/60; tdif=1*ones(1,length(t))'/60;
blood = (t.*exp(-t)+exp(-t)-(t+tdif).*exp(-t-tdif)-exp(-t-tdif))./tdif;

scant = [t, t+tdif]*60;

% analytical solution for each time frame
expt = exp(-t);
analytical = 1/(k2-1)^2*(t.*expt+2*expt+k2*(-t.*expt-expt)-1/k2*exp(-k2*t));
t = t+tdif;
expt = exp(-t);
analytical = analytical - 1/(k2-1)^2*(t.*expt+2*expt+k2*(-t.*expt-expt)-1/k2*exp(-k2*t));
analytical = -analytical * k1 ./ tdif; 

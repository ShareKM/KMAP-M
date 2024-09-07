function [r,p,mdl,ypred] = corr_plot(x,y, ctype, ptype, ltype)

if nargin<3 | isempty(ctype)
    ctype = 'Pearson'; %'Spearman';
end
if nargin<4
    ptype = 1;
end
if nargin<5
    ltype = 'r-';
end
idx = ~isnan(x(:)) & ~isnan(y(:));


x = x(idx);
y = y(idx);
[r,p] = corr(x(:),y(:),'type',ctype);



ymin = min(y)-std(y);
ymax = max(y)+std(y);
xmin = min(x)-std(x);
xmax = max(x)+std(x);

ylim([ymin ymax]); 
xlim([xmin xmax]);

hold on; box on;
xx = linspace(xmin*1.1,xmax*0.98,100);
mdl = fitlm(x,y);
ypred = predict(mdl,xx(:));

switch ptype
    case 1
        plot(xx,ypred,ltype); 
    otherwise
        plot(xx,xx,ltype)
end

plot(x,y,'o','MarkerSize',8, 'MarkerFaceColor','g'); 

ypred = predict(mdl,x(:));



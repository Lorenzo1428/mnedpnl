clc 
clear
close all
r = 0.2;
xc = 0.7;
yc = xc;
x = xc-r:0.01:r+xc;
yp = sqrt(r^2-(x-xc).^2);
p = plot([0.1,0.2],[0.1,0.2],"o",x,yp+yc,x,-yp+yc)
p(1).XData
xlim([-1,1]);
ylim([-1,1]);
axis square;
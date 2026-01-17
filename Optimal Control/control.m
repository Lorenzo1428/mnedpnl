clc
clear 
close all

Vfunc = load("fh_Vfunc.mat");
V = Vfunc.Vt; 
xn = [0.5,0.5]; %punto iniziale
xc = 0.8; %centro del target
yc = -0.5;
gamma = 0; %peso del controllo

Tf = 1; 
h = 0.05;
dx = 0.05;
epsi = 1e-6;
f = @(x,y,a) [0,0] + [ y/(sqrt(x^2+y^2) + epsi), -x/(sqrt(x^2+y^2) + epsi) ]  + a;
rho = linspace(0,1,4);
rho = rho(2:end);
theta = linspace(0,2*pi,16);
[R,T] = ndgrid(rho,theta);
A = [0, 0 ; R(:).*cos(T(:)) , R(:).*sin(T(:))];
x = -1:dx:1;
N = length(x);
Na = size(A,1);
[X,Y] = ndgrid(x);

interp = @(t,s,Vc) (1-t)*(1-s)*Vc(1) + t*(1-s)*Vc(2) + (1-t)*s*Vc(3) + s*t*Vc(4);


xvel = -0.7;
yvel = 0.3;
radius = 0.1;
l = @(s,x,y,a) (gamma*(a(1)^2 + a(2)^2))*( (x - (xc+xvel*s/Tf))^2 + (y - (yc+yvel*s/Tf))^2 >= radius^2);


theta = linspace(0,2*pi,150);
xp = xc + radius*cos(theta);
yp = yc + radius*sin(theta);
f1 = figure;
Cont = contourf(X,Y,V(:,:,1),20); 
hold on
fi = fill(xp,yp,'r', 'FaceAlpha', 0.4, 'EdgeColor', 'r'); 

p = plot(xn(1),xn(2),"k.-",LineWidth=2.5,MarkerSize=15);
title("Optimal Control");
%legend("V function","Target","Trajectory")
xlim([-1,1]);
ylim([-1,1]);
axis square
axis xy

it1max = 100;
it1 = 0;
arg = 1e3*ones(Na,1);
s = 0;
n = 1;
while (xn(1) - (xc+xvel*s/Tf))^2 + (xn(2) - (yc+yvel*s/Tf))^2 > radius^2 && it1 < it1max
    for m = 1:Na
        xn1 = xn + h*f(xn(1),xn(2),A(m,:));
        if xn1(1) >= 1 || xn1(1) <= -1 || xn1(2) >= 1 || xn1(2) <= -1
            continue;
        else
            [xi,yi] = findcell(xn1,dx);
            Ci = round((xi+1)/dx) +1;
            Cj = round((yi+1)/dx) +1;
            Vc = [V(Ci,Cj,n) V(Ci+1,Cj,n) V(Ci,Cj+1,n) V(Ci+1,Cj+1,n)];
            Vm = interp((xn1(1) - xi)/dx,(xn1(2) - yi)/dx,Vc);
            arg(m) = (Vm + h*l(s,xn(1),xn(2),A(m,:)));
        end   
    end
    s = s+h;
    [~, idx] = min(arg);
    xn = xn + h*f(xn(1), xn(2), A(idx,:));
    arg = 1e3*ones(Na,1);
    xp = (xc + xvel*s/Tf) + radius*cos(theta);
    yp = (yc + yvel*s/Tf) + radius*sin(theta); 
    fi.XData = xp;
    fi.YData = yp;
    p.XData = [p.XData xn(1)];
    p.YData = [p.YData xn(2)];
    %Cont.ZData = V(:,:,n+1);
    drawnow;
    pause(0.05);
    it1 = it1+1;
    n = n+1;
end

function [xi,yi] = findcell(x,dx)
        k = (x(1)+1)/dx;
        l = (x(2)+1)/dx;
        xi = floor(k)*dx -1;
        yi = floor(l)*dx -1;
end
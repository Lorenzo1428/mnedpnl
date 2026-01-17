clc
clear 
close all

xn = [0.3,0.1]; %punto iniziale
xc = -0.8; %centro del target
yc = -0.4;
xc1 = 0.8;
yc1 = 0.8;
lambda = 0; %discount factor
gamma = 1; %peso del controllo

h = 0.05;
dx = 0.05;
epsi = 1e-6;
f = @(x,y,a) [0,0] + 0*[ y/(sqrt(x^2+y^2) + epsi), -x/(sqrt(x^2+y^2) + epsi) ]  + a;
rho = linspace(0,1.1,4);
rho = rho(2:end);
theta = linspace(0,2*pi,16);
[R,T] = ndgrid(rho,theta);
A = [0, 0 ; R(:).*cos(T(:)) , R(:).*sin(T(:))];
x = -1:dx:1;
N = length(x);
Na = size(A,1);
[X,Y] = ndgrid(x);

interp = @(t,s,Vc) (1-t)*(1-s)*Vc(1) + t*(1-s)*Vc(2) + (1-t)*s*Vc(3) + s*t*Vc(4);

tol = 1e-8;
itmax = 500;
radius = 0.1;
l = @(x,y,a) (1 + 0*gamma*(a(1)^2 + a(2)^2))*( (x - xc)^2 + (y - yc)^2 >= radius^2 && (x- xc1)^2 + (y - yc1)^2 >= (radius + 0.1)^2);
V = 100*((X-xc).^2 + (Y-yc).^2 >= radius^2 & (X - xc1).^2 + (Y - yc1).^2 >= (radius + 0.1)^2);
Vm = zeros(N);
Vk = Vm;
err = 1;
it = 0;

while err > tol && it < itmax
    for i = 1:N
        for j = 1:N
            Vm(i,j) = 100;

            if (X(i) - xc)^2 + (X(j) - yc)^2 < radius^2 || (X(i)- xc1)^2 + (X(j) - yc1)^2 < (radius + 0.1)^2
                Vm(i,j) = 0;
                continue; % Salta il calcolo per questo punto
            end

            for m = 1:Na
                xs = [X(i,j), Y(i,j)] + h*f(X(i,j),Y(i,j),A(m,:));
                if xs(1) >= 1 || xs(1) <= -1 || xs(2) >= 1 || xs(2) <= -1
                    continue;
                else
                    [xi,yi] = findcell(xs,dx);
                    Ci = round((xi+1)/dx) +1;
                    Cj = round((yi+1)/dx) +1;
                    Vc = [V(Ci,Cj) V(Ci+1,Cj) V(Ci,Cj+1) V(Ci+1,Cj+1)];
                    Vij = (1/(1-h*lambda))*(interp((xs(1) - xi)/dx,(xs(2) - yi)/dx,Vc)+ h*l(X(i,j),Y(i,j),A(m,:)));
                    Vm(i,j) = min(Vm(i,j),Vij);
                end
            end
            Vk(i,j) = Vm(i,j);
        end
    end
    err = norm(Vk(:)-V(:),inf);
    V = Vk;
    it = it+1;
end
save("ih_Vfunc","V");

theta = linspace(0,2*pi,150);
xp = xc + radius*cos(theta);
yp = yc + radius*sin(theta);
xp1 = xc1 + (radius + 0.1)*cos(theta);
yp1 = yc1 + (radius + 0.1)*sin(theta);
f1 = figure;
contourf(X,Y,V,20); 
hold on;
fill(xp,yp,'r', 'FaceAlpha', 0.4, 'EdgeColor', 'r'); 
legend("target 1")
hold on;
fill(xp1,yp1,'r', 'FaceAlpha', 0.4, 'EdgeColor', 'r'); 
legend("target 2")
hold on;

p = plot(xn(1),xn(2),"k.-",LineWidth=2.5,MarkerSize=15);
title("Optimal Control");
legend("V function","Target","Trajectory")
xlim([-1,1]);
ylim([-1,1]);
axis square
axis xy

it1max = 500;
it1 = 0;
arg = 1e3*ones(Na,1);
while (xn(1) - xc)^2 + (xn(2) - yc)^2 > radius^2 && it1 < it1max
    for m = 1:Na
        xn1 = xn + h*f(xn(1),xn(2),A(m,:));
        if xn1(1) >= 1 || xn1(1) <= -1 || xn1(2) >= 1 || xn1(2) <= -1
            continue;
        else
            [xi,yi] = findcell(xn1,dx);
            Ci = round((xi+1)/dx) +1;
            Cj = round((yi+1)/dx) +1;
            Vc = [V(Ci,Cj) V(Ci+1,Cj) V(Ci,Cj+1) V(Ci+1,Cj+1)];
            Vm = interp((xn1(1) - xi)/dx,(xn1(2) - yi)/dx,Vc);
            arg(m) = (1/(1-h*lambda))*(Vm + h*l(xn(1),xn(2),A(m,:)));
        end   
    end
    [~, idx] = min(arg);
    xn = xn + h*f(xn(1), xn(2), A(idx,:));
    arg = 1e3*ones(Na,1);
    p.XData = [p.XData xn(1)];
    p.YData = [p.YData xn(2)];
    drawnow;
    pause(0.05);
    it1 = it1+1; 
end

function [xi,yi] = findcell(x,dx)
        k = (x(1)+1)/dx;
        l = (x(2)+1)/dx;
        xi = floor(k)*dx -1;
        yi = floor(l)*dx -1;
end
clc
clear 
close all

xn = [-0.5,-0.5]; %punto iniziale
xc = 0.8; %centro del target
yc = -0.5;
gamma = 0.4; %peso del controllo

Tf = 1;
dt = 0.05;
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

s = 0:dt:Tf;
Nt = length(s);
s = s(Nt-1:-1:1);

V = 1*((X-(xc+xvel*Nt/Tf)).^2 + (Y-(yc+yvel*Nt/Tf)).^2 >= radius^2);
Vm = zeros(N);
Vk = Vm;

Vt = zeros(N,N,Nt);
Vt(:,:,Nt) = V;
for n = 1:Nt-1
    for i = 1:N
        for j = 1:N
            Vm(i,j) = 100;

            if (X(i,j) - (xc+xvel*s(n)/Tf))^2 + (Y(i,j) - (yc+yvel*s(n)/Tf))^2 <= radius^2
                Vm(i,j) = 0;
                continue;
            end

            for m = 1:Na
                xs = [X(i,j), Y(i,j)] + dt*f(X(i,j),Y(i,j),A(m,:));
                if xs(1) >= 1 || xs(1) <= -1 || xs(2) >= 1 || xs(2) <= -1
                    continue;
                else
                    [xi,yi] = findcell(xs,dx);
                    Ci = round((xi+1)/dx) +1;
                    Cj = round((yi+1)/dx) +1;
                    Vc = [V(Ci,Cj) V(Ci+1,Cj) V(Ci,Cj+1) V(Ci+1,Cj+1)];
                    Vij = (interp((xs(1) - xi)/dx,(xs(2) - yi)/dx,Vc)+ dt*l(s(n),X(i,j),Y(i,j),A(m,:)));
                    Vm(i,j) = min(Vm(i,j),Vij);
                end
            end
            Vk(i,j) = Vm(i,j);
        end
    end
    V = Vk;
    Vt(:,:,Nt-n) = V;
end
save("fh_Vfunc","Vt");

function [xi,yi] = findcell(x,dx)
        k = (x(1)+1)/dx;
        l = (x(2)+1)/dx;
        xi = floor(k)*dx -1;
        yi = floor(l)*dx -1;
end
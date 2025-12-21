clc
clear 
close all

T = 2;
dx = 0.01;
x = -1:dx:1;
Nx = length(x);
I = 2:Nx-1;
L = I - 1;
R = I + 1;

u0 = @(x) zeros(length(x),1);%1 - x.^2;
dt = dx/2;
Nt = T/dt;

f = @(U) 0.5*U.^2 - 0.5;

%F = @(U,V) f(max(U,0)) + f(min(V,0));
%lax f
F = @(U,V) 0.5*(f(U) + f(V)) + 0.5*(dx/dt)*(U - V);
U = u0(x);

p = plot(x,U);
ylim([-0.5,1.5]);
%waitforbuttonpress;
n = 0;
while n < Nt
    n = n+1;
    dUl = (U(I) - U(L))/dx;
    dUr = (U(R) - U(I))/dx;
    U(I) = U(I) - dt*( F(dUl,dUr) );
    p.YData = U;
    pause(0.01);
    drawnow;
end








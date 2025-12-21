clc
clear 
close all

T = 2;
dx = 0.01;
x = -1:dx:1;
Nx = length(x);

uL = 1;
uR = 0.2;

u0 = @(x) uL*(x < 0) + uR*(x >= 0);
fb = @(u) 0.5*u.*u;
mu = @(u) 0.5*u.*u;
csi = @(u) (1/3)*u.^3;

U = u0(x);
dt = dx/max(abs(uL),abs(uR));
Nt = floor(T/dt);
lambda = dt/dx;

I = 1:Nx;
L = I - 1;
R = I + 1;
L(1) = 1;
R(end) = Nx - 1;

n = 0;

Us = zeros(1,Nx);
M = zeros(1,Nx);

p = plot(x,U,x,M);
xlim([-1,1]);
ylim([-2*dx,2*dx]);
%legend("u","mu")

gudonov
while n < Nt
    n = n+1;
    Us(I) = U(I) - lambda*(F(U(I),U(R),fb) - F(U(L),U(I),fb));
    M(I) = mu(Us(I)) - mu(U(I)) + lambda*(F(U(I),U(R),csi) - F(U(L),U(I),csi));
    U = Us;
    p(1).YData = U;
    p(2).YData = 100*M;
    drawnow;
end

%lax F
% while n < Nt
%     n = n+1;
%     Us(I) = U(I) - lambda*(F1(U(I),U(R),dx,dt,fb) - F1(U(L),U(I),dx,dt,fb));
%     M(I) = mu(Us(I)) - mu(U(I)) + lambda*(F(U(I),U(R),csi) - F(U(L),U(I),csi));
%     U = Us;
%     p(1).YData = U;
%     p(2).YData = M;
%     drawnow;
% end

function flux = F1(U,V,dx,dt,f)
    flux = 0.5*(f(U) + f(V)) + 0.5*(dt/dx)*(U - V);
end

function flux = F(U1,U2,f)
    flux = max(f(U1),f(U2)).*(U1 >= U2) + f(U2).*(U1 <= U2 & U2 < 0) + f(U1).*(U1 < U2 & U1 > 0);  
end
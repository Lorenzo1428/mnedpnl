clc
clear
close all

T = 7;
dx = 0.005;
x = -1:dx:4;
Nx = length(x);

uL = 0;
uR = 0.5;

u0 = @(x) uL*(x < 0) + uR*(x >= 0);
fb = @(u) u.*(1-u);

U = u0(x);
dt = dx;
Nt = floor(T/dt);
lambda = dt/dx;

I = 1:Nx;
L = I - 1;
R = I + 1;
L(1) = 1;
R(end) = Nx - 1;

p = plot(x,U);
xlim([-1,4]);
ylim([-1,2]);

%gudonov
n = 0;
while n < Nt
    if n*dt < 1 || n*dt > 1.5
        n = n+1;
        U(I) = U(I) - lambda*(F(U(I),U(R),fb) - F(U(L),U(I),fb));
        p.YData = U;
        drawnow;
    else 
        M = 400;
        n = n + 1;
        U(M) = 0;
        U(I) = U(I) - lambda*(F(U(I),U(R),fb) - F(U(L),U(I),fb));
        p.YData = U;
        drawnow;
        pause(0.01)
    end
end

function flux = F(U1,U2,f)
    flux = max(f(U1),f(U2)).*(U1 >= U2) + f(U2).*(U1 <= U2 & U2 < 0) + f(U1).*(U1 < U2 & U1 > 0);  
end
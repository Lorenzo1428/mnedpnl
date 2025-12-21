clc
clear 
close all

T = 1;
dx = 0.001;
x = -1:dx:1;
Nx = length(x);

uL = 0.5;
uR = 0.2;

u0 = @(x) uL*(x < 0) + uR*(x >= 0);
fb = @(u) 0.5*u.*u;

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

p = plot(x,U);
xlim([-1,1]);
ylim([-1,2]);

%upwind
% while n < Nt
%     n = n+1;
%     U(I) = U(I) - lambda*( (fb(U(R)) - fb(U(I))).*(U(I) < 0 ) + (fb(U(I)) - fb(U(L))).*( U(I) > 0 ) );
%     p.YData = U;
%     drawnow;
% end

%lax f
% while n < Nt
%     n = n+1;
%     U(I) = 0.5*(U(R) + U(L)) - 0.5*lambda*(fb(U(R)) - fb(U(L)));
%     p.YData = U;
%     drawnow;
% end

%mac cormack
Us = zeros(1,Nx);
while n < Nt
    n = n+1;
    Us(I) = U(I) - lambda*(fb(U(R)) - fb(U(I)));
    U(I) = 0.5*(U(I) + Us(I)) - 0.5*lambda*(fb(Us(I)) - fb(Us(L)));
    p.YData = U;
    drawnow;
end

%lax w
% Up = zeros(1,Nx);
% Um = zeros(1,Nx);
% while n < Nt
%     n = n+1;
%     Up(I) = 0.5*(U(R) + U(I)) - 0.5*lambda*(fb(U(R)) - fb(U(I)));
%     Um(I) = 0.5*(U(L) + U(I)) - 0.5*lambda*(fb(U(I)) - fb(U(L)));
%     U(I) = U(I) - lambda*(fb(Up(I)) - fb(Um(I)));
%     p.YData = U;
%     drawnow;
% end

%gudonov
% while n < Nt
%     n = n+1;
%     U(I) = U(I) - lambda*(F(U(I),U(R),fb) - F(U(L),U(I),fb));
%     p.YData = U;
%     drawnow;
% end
% 
% function flux = F(U1,U2,f)
%     flux = max(f(U1),f(U2)).*(U1 >= U2) + f(U2).*(U1 <= U2 & U2 < 0) + f(U1).*(U1 < U2 & U1 > 0);  
% end


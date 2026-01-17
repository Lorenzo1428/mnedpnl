clc
clear 
close all

T = 1;
dt = 0.01;
dx = 0.01;
Nt = floor(T/dt);
h = 0.01;
x = -1:dx:1;
a = -1:h:1;
L = @(a) 1;%0.5*abs(a)^2;

U0 = x.^2;
U = U0;
U1(:,1) = U0;
for  n = 1:Nt    
    U_new = inf(size(U));
    for i = 1:length(a)
        xs = x - a(i)*dt;
        bound_index = xs >= x(1) & xs <= x(end);
        xs_bound = xs(bound_index);
        pos = (xs_bound-x(1))/dx +1;
        k = floor(pos);
        k = max(1, min(length(x)-1, k));

        w = pos - k;
        Us = (1-w).*U(k) + w.*U(k+1);
        Ui = dt*L(a(i)) + Us;
        U_new(bound_index) = min(U_new(bound_index),Ui);
    end 
    mask_inf = isinf(U_new);
    if any(mask_inf)
        U_new(mask_inf) = U(mask_inf); 
    end

    U_new(1) = 0;
    U_new(end) = 0;
    U = U_new;
    U1(:,n+1) = U;
end

%f = @(x,t) (abs(x) - t*0.5).*(abs(x) > t) + (0.5*abs(x).^2/t).*(abs(x) <= t);
%f = @(x,t) (abs(x) - t).*(abs(x) > t) + (zeros(1,length(x))).*(abs(x) <= t);
f = @(x,t) 1 - abs(x);

err = norm(U1(:,end)' - f(x,1),inf);

p = plot(x,U1(:,1),x,f(x,1),"--");
legend("U approx","U esatta");
title("T = 0");
xlim([x(1) x(end)]);
ylim([0 1]);
for n = 1:Nt
    p(1).YData = U1(:,n+1);
    title("T = " + n*dt);
    drawnow
    pause(0.1)
end
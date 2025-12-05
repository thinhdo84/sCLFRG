function xd=RK4(dynamics, x,u,h,T)
N = round(T/h);
for i =1:N
    k1 = dynamics(x,u);
    k2 = dynamics(x+0.5*h*k1,  u);
    k3 = dynamics(x+0.5*h*k2,  u);
    k4 = dynamics(x+h*k3,  u);
    x = x + (1/6)*h*(k1 + 2*k2+  2*k3 + k4);
end
xd = x;
end
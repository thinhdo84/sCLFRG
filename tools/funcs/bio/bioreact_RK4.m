function xi=bioreact_RK4(xi_p,u,h,Ts, para)
Nstep = round(Ts/h);
xi = xi_p;
for i =1: Nstep
    k1 = bioreact(xi,u, para);
    k2 = bioreact(xi+0.5*h*k1,  u, para);
    k3 = bioreact(xi+0.5*h*k2,  u, para);
    k4 = bioreact(xi+h*k3,  u, para);
    xi = xi   +  (1/6)*h*(k1 + 2*k2+  2*k3 + k4);
end
end
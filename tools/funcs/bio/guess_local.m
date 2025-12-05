function useq = guess_local(K,x,Npred,h,tau, para, umax, umin)
useq = sat(K*x,umax,umin);
nu = size(useq,1);
useq = [useq, zeros(nu,Npred-1)];
for i =1:Npred-1
%     x = bioreact(x,useq(:,i), para)*tau + x;
    x = bioreact_RK4(x,useq(:,i),h,tau, para);
    useq(:,i+1) = sat(K*x,umax,umin);;
end
end

function us = sat(u, up, low)
 us = min(up, max(low, u));
end
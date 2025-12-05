function sys = getCSE(ns,mu,delta,k)
%getCSE return the matrices describing the CSE system
%   See "On complexity reduction in a variable terminal set
%   setpoint-tracking MPC scheme" for more details
%   ns: the number of strings
%   mu, delta, k is the system parameters ...


%%%%% temporary matrices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mtmp = -2*eye(ns); Mtmp(1,1)=1; Mtmp(ns,ns)=1;
Ntmp = -[zeros(ns-1,1)',0;
eye(ns-1),zeros(ns-1,1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kc=1/2*(Mtmp + 2*Ntmp + (Mtmp + 2*Ntmp)')*k;
Dc = zeros(ns,2); Dc(1,1) = 1; Dc(ns,2) = -1;
Mc = mu * eye(ns);
Lc = delta*eye(ns);
A = [zeros(ns,ns), eye(ns,ns);
    -Mc^(-1)*Kc, -Mc^(-1)*Lc];
B = [zeros(ns,2); Mc^(-1)*Dc];
sys.A = A;
sys.B = B;
C = zeros(2, 2*ns);
C(1,1) = 1; C(2,2)=1;
sys.C = C;
end

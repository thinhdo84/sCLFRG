function gamma_list = findLevelSetSCLF(V_para, a, b, c)
% This code returns the "largest" level set:
% {x: V_para(x-c)<= Gamma}
% V_para(x) is a convex parameterized SCLF
N = numel(V_para.V_list);
level_list = zeros(1,N);
gamma_list = zeros(1,N);
for i = 1:N
    V = V_para.V_list{i};
    if isa(V,'Lyapunov_quad')
        gamma_list(i) = (a'*c - b)^2/(a'*V.P^(-1)*a);
    elseif isa(V,'Lyapunov_2p')
         M = size(V.F,1);
         rho = -(a' * c - b)/norm(V.F * ((V.F' * V.F)\a));
%          gamma_list(i) =( rho * M^(1/(2*V.p)-1/2)  )^(2*V.p);
 gamma_list(i) =( rho * M^(1/(2*V.p)-1/2)  )^(2*V.p);
        
    else
        error('undefined class of Lyapunov function')
    end
    
end

end
function Phi = CLFCondition(Vx,sys, beta, kappa,Bm,a)
nx = size(sys.A,1);
X=sdpvar(nx,1);
constr=Vx.eval(X)<=kappa;

%Compute Vx's gradient
N=Vx.N;
cnt =1;
for i = 1:N
    if Vx.alpha_list(i)>1e-12
        Vi = Vx.V_list{i};
        if cnt == 1
            if isa(Vi,'Lyapunov_quad')
                dVx = 2*(Vi.P*X)';
            end
            if isa(Vi,'Lyapunov_2p')
                dVx_tmp = [];
                for k=1:nx
                    dVk=0;
                    for j=1:Vi.nc
                        dVk=dVk + 2*Vi.p*Vi.F(j,k)*(Vi.F(j,:)*X)^(2*Vi.p-1);
                    end
                    dVx_tmp = [dVx_tmp, dVk];
                end
                dVx=dVx_tmp;
            end
            cnt = cnt+1;
        else
            if isa(Vi,'Lyapunov_quad')
                dVx = dVx+ 2*(Vi.P*X)';
            end
            if isa(Vi,'Lyapunov_2p')
                dVx_tmp = [];
                for k=1:nx
                    dVk=0;
                    for j=1:Vi.nc
                        dVk=dVk + 2*Vi.p*Vi.F(j,k)*(Vi.F(j,:)*X)^(2*Vi.p-1);
                    end
                    dVx_tmp = [dVx_tmp, dVk];
                end
                dVx=dVx+dVx_tmp;
                
            end
            cnt = cnt+1;
        end
    end
end
assign(X,zeros(nx,1)+0.0001);
objective= dVx*sys.A*X  -   supportBall(sys.B'*dVx', Bm,a)  + beta * Vx.eval(X);
options=sdpsettings('solver','+ipopt','ipopt.max_iter',6000,'usex0',1);
% options=sdpsettings('usex0',1);
% options=sdpsettings('solver','sdpt3');
% options=sdpsettings('solver','cplex','usex0',1);
res=optimize(constr,-objective,options);
disp(res.info)
% objective
Phi = value(objective);
% disp('--------')
% disp(Vx.eval(value(X)))
% disp('--------***')
% disp(value(X))
% disp('--------***')

end

%% functions
function val = supportBall(x,Bm,a)
val = a * sqrtm(x'*Bm'*Bm*x);

end
function B = find_transform(Q)
[U,S,V] = svd(Q);
B = ((U*sqrt(S)*V')')^-1;
end

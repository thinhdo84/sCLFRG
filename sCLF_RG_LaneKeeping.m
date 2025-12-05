clc
clear 
close all
addpath(genpath('tools'))
addpath(genpath('funcs'))
%% Model
% Lane Keeping Using Lyapunov Function-Based Reference Governor: An Optimization-Free Approach
% Xiao Li et al 2025
l1 = 5; l2=2;
w = 3.5;
v = 10;
delta_max = 30*pi/180;
xie = [0; 0]; ue = 0;
f_field = @(xt)[v*sin(xt(2)); 0];
g_field = @(xt) [0; v/l1];
LaneKeeping = @(xi,u) [v*sin(xi(2)); v*u/l1];
% Input constraint set
umax = tan(delta_max);
umin = -tan(delta_max);
U = Polyhedron('ub',umax,'lb',umin);
% state constraint set S in (5)
bs = 0.5*(w-l2);
AS = [1 0;
     -1 0;
      1 l1;
     -1 -l1];
S = Polyhedron('A',AS,'b',bs*ones(4,1));
%% Linearize the model around xie,ue
A = [0, v*cos(xie(2));
    0,          0];
B = [0 ;v/l1];
if rank(ctrb(A,B))~=2
    error uncontrollable
else
   disp('The system is locally controllable!') 
end
%% Formulating the MPC controller
%%%%%%  Ref: A Quasi-Infinite Horizon Nonlinear Model Predictive Control Scheme with Guaranteed Stability
%%%%%%  H Chen, F Allgöwer 1998
% Find the terminal ingredients
% NMPC design local state feedback tuning
Q = diag([10 30]);
R = 20;
K = -lqr(A,B,Q,R);
AK = A+B*K;
% pause
disp('______________________________')
% Choose kappa to satisfy the condition
kappa = -max(real(eig(AK)))*0.855;
if kappa >= -max(real(eig(AK)))
    error('kappa >= -max(real(eig(AK)))')
end
P=lyap((AK+kappa*eye(2))',Q+K'*R*K-0.001*eye(2));
disp(P)
disp('______________________________')
% Input constraints in the state space 
clc
Su = Polyhedron('A', [K; -K], 'b',[umax; -umin]);
S_full = S.intersect(Su); % intersection of input and state constraint
alpha_1 = Ellipse_inner_fixed_shape(S_full,P);
% Compute the bound by bisection
FcnPhi = @(x,AKt,Kt) LaneKeeping(x,Kt*x) - AKt*x;
%  Set lower and upper bounds for alpha for bisection
alpha = alpha_1;
alpha_ub = alpha_1;
alpha_lb = 0;
max_iter = 200;                      % Max iterations for bisection
opt = optimset('MaxFunEvals',10000,'MaxIter',10000,'Display','off'); % Options for fmincon
% take some random points for intial guess inside the constraint set
x2_init = cprnd(20,S.A,S.b)';
alph_pre = -100;
% start bisection
for i= 1:max_iter
    fval = 1e50;
    % inner for loop to check multiple initializations for fmincon
    for j = 1:size(x2_init,2)
        pts = random_points_on_ellipsoid(P,alpha,[0;0],25);
        ell_approx = Polyhedron('V',pts');
        [x2,val] = fmincon(@(xt)-(xt'*P*FcnPhi(xt,AK,K) - kappa*xt'*P*xt)/(xt'*P*xt),...
            x2_init(:,j),ell_approx.A,ell_approx.b,[],[],[],[],[],opt);
        % Take minimum value
        fval = min(fval,val);
    end
    % Check condition: Is optimal value of (7) nonpositive  and update upper and lower bounds for alpha  
    fval = -fval; % Correct sign due to maximization
    if fval>=0 % condition not satisfied
        alpha_ub = alpha;
    else % condiction satisfied
        alpha_lb = alpha;
        fprintf('i= %d, alpha = %.3f\n', i, alpha)
        if abs(alph_pre-alpha)/alpha<=1e-6
            disp('converge')
            break;
        end
        alph_pre = alpha;
    end
    alpha = (alpha_lb + alpha_ub)/2; % bisection
    
end
alpha = alpha_lb; % Take best known feasible value for alpha
disp(['alpha = ',num2str(alpha)])
%% Create the MPC law as a function
clear opts
[nx, nu] = size(B);
tau = 0.02; % sampling time
Npred = 250;
% create NMPC optimizer
solver=casadi.Opti();
% variables and parameters
xt=solver.variable(nx,Npred+1);
ut=solver.variable(nu,Npred);
xinit=solver.parameter(nx,1);
objective=0;       % initialize the objective
% add constraints to the solver object
solver.subject_to(xt(:,1)==xinit);% the state sequence has to start from the initial state
for k=1:Npred
    % EAS model
    solver.subject_to(...
        xt(:,k+1)== tau*LaneKeeping(xt(:,k),ut(:,k)) + xt(:,k) ... % here for some reason xe is skipped
        );
    % input bounds
    solver.subject_to(umin-ue<=ut(:,k));
    solver.subject_to(ut(:,k)<=umax-ue);
    solver.subject_to(S.A*xt(:,k)<=S.b);
    objective=objective...
        +tau*(xt(:,k))'*Q*(xt(:,k))...
        +tau*(ut(:,k))'*R*(ut(:,k));
end
% add to the objective the terminal cost
objective=objective+(xt(:,Npred+1))'*P*(xt(:,Npred+1));
solver.subject_to(xt(:,Npred+1)'*P*xt(:,Npred+1)<=alpha);
% solver.subject_to(Xf0.A*xt(:,Npred+1)<=Xf0.b); %:) Inner approx. of Xf
% add the objective to the solver object
solver.minimize(objective);
opts.ipopt.print_level=0;
opts.print_time=0;
opts.ipopt.sb='yes';
opts.ipopt.max_iter = 2000;
solver.solver('ipopt',opts);

omega = @(x)solveMPC(x,xinit,ut,solver);
%% Collect points from the grids
tic;
cost_g = @(x,u)(x'*Q*x + u'*R*u); % stage cost

npts = 40; % grid size

x1max = 1;
x2max = 0.4;
x1train =linspace(-x1max,x1max,npts) ;
x2train = linspace(-x2max,x2max,npts);
gX = [];
gU = [];
g_cost= [];
xtrain = [];
cnt = 1;
for x1 = x1train
    for x2 = x2train
        X = [x1;x2];
        if all(S.A*X<=S.b)
            xtrain = [xtrain, X];
        end
    end
end
for i =1:size(xtrain,2)
    X = xtrain(:,i);
    umpc = omega(X);
    if ~isnan(umpc)
        gX = [gX, X];
        gU = [gU, umpc];
        g_cost = [g_cost, cost_g(gX(:,end),gU(:,end))];
        disp(cnt)
        cnt = cnt+1;
    else
        error('infeasible at interior')
    end
end
xtrainBound =  SamplePolytopeBoundary(10,S.A,S.b,'Nsample',500);
for i =1:size(xtrainBound,2)
    X = xtrainBound(:,i);
    umpc = omega(X);
    if ~isnan(umpc)
        gX = [gX, X];
        gU = [gU, umpc];
        g_cost = [g_cost, cost_g(gX(:,end),gU(:,end))];
        disp(cnt)
        cnt = cnt+1;
    else
        error('infeasible at boundary')
    end
end
dataTime = round(toc,5);
% beep
fprintf('Data collected in %f (s)\n',dataTime)
M = size(gX,2);
%% basis function setup
V = {};
V{1}=Lyapunov_quad(P); % quadratic function
W0 =[1 0;
     0 10;
     1 5]*0.5;
eta_list = 1;
ang_list = 0:0.2:pi;
for eta=eta_list
    for ang =ang_list
        Rota = [cos(ang), -sin(ang); 
            sin(ang), cos(ang)];
            V{end+1}=Lyapunov_2p(W0*Rota'*eta,4); % gather some 2p-norm functions
    end
end
N = numel(V);
l_list = ones(1,N);
fprintf('size of the basis function set %d\n',N)
%% LP matrix arrangement
Ain_eq = [];
bin_eq = [];
for i =1:M
    disp(i)
    [atmp, btmp] = getLinearConstraint_NL(V,LaneKeeping,gX(:,i),gU(:,i),g_cost(i));
    Ain_eq =  [Ain_eq; atmp];
    bin_eq = [bin_eq; btmp];
end
Ain_eq = round(Ain_eq,15);
bin_eq = round(bin_eq,15);
%% Solve the linear program
clc
tic
alp2 = SolveCLF_new(Ain_eq,bin_eq);
toc
rela_tol = 0.5e-8/100; % discard the coefficient smaller than this value
tol = sum(alp2)*rela_tol;
disp('________________________________________________')
fprintf('number of non-zero coefficients: %d\n',numel(alp2(alp2>tol)))
tol = sum(alp2)*rela_tol;
idx_clf=find(alp2>tol);
V_full = Lyapunov_para(V(idx_clf),alp2(idx_clf)); % the resulted subCLF
%% sCLF-induced Controller
gain_ctrl=-0.5*R^(-1);
ubar = @(xt) sat_scalar(gain_ctrl*g_field(xt)'*V_full.grad(xt)',umax,-umax);

%% simulation test
tsim = 30;
tau = 0.005;
Nsim = round(tsim/tau);
tt = linspace(0,tsim,Nsim);
RefLift = [1; 0];
lv_set_input = 1.3; 
X0 = [-0.725; -0.0]; % initial condition
clc
ref = zeros(1,Nsim);
reflist=[0, 0.5 ,-2];
switchReftime=floor(linspace(1,Nsim,size(reflist,2)+1));
for k = 1:numel(switchReftime)-1
    for i = switchReftime(k):switchReftime(k+1)
        ref(:,i) = reflist(:,k);
    end
    
end
Xsim= zeros(nx,Nsim+1);
Xsim(:,1)=X0;
usim=zeros(nu,Nsim);

eta_erg =0.1;
kappa = 25.0;
ref_virtual = X0(1);
ref_list = zeros(1,Nsim);
px = zeros(1,Nsim+1);
clc
for i = 1:Nsim
    xcurrent = Xsim(:,i);
    tic
    for j = 1:size(S.A,1)
        % get the safe hyperplanes wrt the state constraints
        a_safe = S.A(j,:); b_safe = S.b(j);
        % find the safest level wrt the state constraints
        % Do not forget to lift the space from ref space to state space
%         lvl_list = findLevelSetSCLF(V_full, a_safe', b_safe, xcurrent);
        lvl_list = findLevelSetSCLF(V_full, a_safe', b_safe,  RefLift*ref_virtual);
        if j == 1
            Gamma_state = lvl_list * V_full.alpha_list';
        else
            Gamma_state = min(lvl_list * V_full.alpha_list',Gamma_state);
        end
    end
    Gamma = min(Gamma_state,lv_set_input);
    Delta =  Gamma - V_full.eval(Xsim(:,i)-RefLift*ref_virtual); % dynamic safety margin
    rho_field = (ref(:,i) - ref_virtual)/max(norm(ref(:,i) - ref_virtual),eta_erg); % attraction field
    ref_virtual = ref_virtual + tau*kappa*Delta*rho_field;
    usim(:,i) =  ubar(Xsim(:,i)-RefLift*ref_virtual);
    if isnan(usim(1,i)) || Delta< 0
        error NAN
    end
    comptime(i)=toc;
    Xsim(:,i+1) = RK4(LaneKeeping,Xsim(:,i),usim(:,i)+ue,tau,tau);
    ref_list(:,i)= ref_virtual;
    px(i+1)=px(i)+v*tau;
    Gamma_list(i) = Gamma;
end


%%
Slvsamp = 1;
idx_lv_samp=floor(linspace(1,Nsim,Slvsamp));
close all
figure
plot(px,Xsim(1,:),'b','linewidth',2)
hold on
% yline(bs,'b--')
% yline(-bs,'b--')
plot(px(1:Nsim),ref_list(1:Nsim),'linestyle','--', 'linewidth',2.5,...
    'color',[0, 161, 48]/255)
plot(px(1:Nsim),ref(1:Nsim),'r-.', 'linewidth',1.25)
yline(w/2,'k-','linewidth',2)
yline(-w/2,'k-','linewidth',2)
% 
for i = floor(linspace(1,Nsim,100))
Draw_rectangle(px(i),Xsim(1,i),Xsim(2,i),l1,l2,'b',...
    'edgealpha',0,'facealpha',0.05)
% axis equal
end

% ylim([-(w/2),w/2])

xlim([min(px),max(px)])
xlabel('x (m)')
ylabel('y (m)')


S_aug_time = Polyhedron('A',[[zeros(size(S.A,1),1),S.A]; [[1;-1], zeros(2,2)]],...
    'b', [S.b;max(tt);0]);
idx_time_samp=floor(linspace(1,Nsim,500));

close all
figure
subplot(3,2,[1 2])
plot(px(idx_time_samp),Xsim(1,idx_time_samp),'b','linewidth',2)
hold on
plot(px(idx_time_samp),ref_list(idx_time_samp),'linestyle','--', 'linewidth',2.5,...
    'color',[0, 161, 48]/255)
plot(px(idx_time_samp),ref(idx_time_samp),'r-.', 'linewidth',1.25)
yline(w/2,'k-','linewidth',2)
yline(-w/2,'k-','linewidth',2)
% 
for i = floor(linspace(1,Nsim,100))
Draw_rectangle(px(i),Xsim(1,i),Xsim(2,i),l1,l2,'b',...
    'edgealpha',0,'facealpha',0.05)
end

% ylim([-(w/2),w/2])

xlim([min(px),max(px)])
xlabel('x (m)')
ylabel('y (m)')

subplot(3,2,4)
plot(tt(idx_time_samp),usim(idx_time_samp))
hold on
grid on
plot([0,max(tt)],umax*[1 1],'linewidth',1,'color','k')
plot([0,max(tt)],-umax*[1 1],'linewidth',1,'color','k')
ylim([-1 1]*1.1*umax)

ax=subplot(3,2,6);
plot(tt(idx_time_samp), Xsim(1,idx_time_samp),'b','linewidth',1.5)
hold on
plot(tt(idx_time_samp), ref_list(1,idx_time_samp),'--','color',[0, 161, 48]/255,...
    'linewidth',1.5)
plot(tt(idx_time_samp),ref(idx_time_samp),'r-.', 'linewidth',1.25)


subplot(3,2,[3 5])
plot3(tt(idx_time_samp), Xsim(1,idx_time_samp),Xsim(2,idx_time_samp),'b','linewidth',2)
hold on
plot3(tt(idx_time_samp), ref_list(1,idx_time_samp),tt(idx_time_samp)*0,'--','color',[0, 161, 48]/255,...
    'linewidth',2)
plot(S_aug_time,'alpha',0.01,'edgealpha',0.3)
grid off


list_snaps=floor(Nsim./[25, 14.5, 2.8, 1.38]);
figure

for j = 1:numel(list_snaps)
% 13,4,1.8,1.3
subplot(2,2,j)
idx = list_snaps(j);
plot(S,'alpha',0.025)
hold on

xcurrent =Xsim(:,idx);
for j = 1:size(S.A,1)
    % get the safe hyperplanes wrt the state constraints
    a_safe = S.A(j,:); b_safe = S.b(j);
    lvl_list = findLevelSetSCLF(V_full, a_safe', b_safe,  RefLift*ref_list(idx));
    if j == 1
        Gamma_state = lvl_list * V_full.alpha_list';
    else
        Gamma_state = min(lvl_list * V_full.alpha_list',Gamma_state);
    end
end
Gamma = min(Gamma_state,lv_set_input);
cttmp = ComputeContour2D(500,V_full,Gamma);
plot(cttmp(1,:)+ref_list(idx),cttmp(2,:), 'color', [0, 161, 48]/255,'linewidth',1)
plot(Xsim(1,1:idx),Xsim(2,1:idx),'b','linewidth',2)
scatter(Xsim(1,idx),Xsim(2,idx),85,'markeredgecolor','b','marker','x','linewidth',1.5)
scatter(ref(1,idx),0,35,'filled','markeredgecolor',[255, 86, 74]/255,'marker','d',...
    'linewidth',0.5,'markerfacecolor',[255, 86, 74]/255)
scatter(ref_list(idx),0,25,'filled','markerfacecolor',[0, 161, 48]/255)


set(gca,'TickLabelInterpreter','latex','fontsize',14)
xlabel('$\xi_1$ (m)', 'interpreter','latex','fontsize',22)
ylabel('$\xi_2$ (rad)', 'interpreter','latex','fontsize',22)
end

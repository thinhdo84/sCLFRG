clc
clear
close all
addpath(genpath('tools'))
addpath(genpath('funcs'))
addpath(genpath('data'))
% to try with
% D:\PhD_data\TP_works\Quanser\Courseware Resources\Fundamental_Labs\Pendulum\State_Space_Modeling\Documentation
%% load data
load('backup_14-May-2024_19h52.mat')
%% Filter out some CLF
idx_filter = find(V_full_iter{end}.alpha_list>0.1e-4);
V_full=Lyapunov_para(V_full_iter{end}.V_list(idx_filter),...
    V_full_iter{end}.alpha_list(idx_filter));
fprintf('the number of basis functions used %d\n',numel(idx_filter));
%% Gather up the basis
order_p = [];
% this is because of the full rowrank condition ...
% Otherwise I would just have to recompute the SCLF again
alp_list = [];
F_list = [];
V_list_new = {};
alpha_new = [];
for i = 1:numel(V_full.V_list)
    if isa(V_full.V_list{i}, 'Lyapunov_2p')
        %     order_p = [order_p, V_full.V_list{i}.p];
        % lucky enough, all p = 2;
        alp_list = [alp_list, 1];
        F_list =[F_list ;
            V_full.V_list{i}.F*V_full.alpha_list(i)^(1/(2*V_full.V_list{i}.p))];
    else
        V_list_new{end+1} =  V_full.V_list{i};
        alpha_new = [alpha_new, V_full.alpha_list(i)];
    end
end
V_list_new{end+1} = Lyapunov_2p(F_list,2);
alpha_new = [alpha_new, 1];
V_full_new = Lyapunov_para(V_list_new,alpha_new);
%%
% Xtest = rand(6,1);
% V_full.eval(Xtest) - V_full_new.eval(Xtest)

V_full = Lyapunov_para({V_list_new{1},V_list_new{end}},alpha_new([1 5]));
%% grid the space
x1_range = linspace(-1,1,10);
x3_range = linspace(-1,1,10);
x2_range = linspace(-0.5,0.5,10);
x4_range = linspace(-0.5,0.5,10);
%% Random verification based on support functional of the input set
% clc
U = Polyhedron('ub', umax-ue, 'lb', umin-ue);
g_field = @(x) [0, 0;
    -sin(x(5)), epsilon*cos(x(5));
    0, 0;
    cos(x(5)), epsilon*sin(x(5));
    0, 0;
    0, 1];
f_field = @(x) [x(2);
    -sin(x(5));
    x(4);
    cos(x(5))-1;
    x(6);
    0];
% NSample = 100000;
%
% % alpha = 0.1;
% lv_alp = 480;
% checkres = checkInvariance(V_full, error_center, lv_alp, f_field, g_field, U, NSample);
% if checkres
%     disp('Verfied invariance')
% else
%     warning('The level set is not positively invariant')
% end

%% Verification based on boundary sampling
% the idea is to sample exterior points, and project them into the boundary

% Formulate the control based on the SCLF
ut_max =umax-ue;
ut_min = umin-ue;
gain_ctrl=-0.5*R^(-1);
ubar = @(xt) sat_scalar(gain_ctrl*g_field(xt)'*V_full.grad(xt)',ut_max,ut_min);

error_center = zeros(6,1);
% Sample some intial conditions
clc
% Options for fmincon solver
opts = optimoptions('fmincon', 'Display', 'none', ...
    'Algorithm', 'interior-point', 'MaxFunctionEvaluations', 1e4);
lv_set_input = 5; %150 was not ok, 100 seems fine, just take 80 for care, 10 is good, 2 is good
[ub, lb] = bounding_box_levelset(V_full, error_center, lv_set_input);
% Expanded bounding box
gam = 1;
SampleBox = Polyhedron('ub', ub+gam, 'lb', lb-gam);
NSample = 400;
Samples = cprnd(NSample,SampleBox.A,SampleBox.b)'; % NSample vector stacked side-by-side
fprintf('%d random points sampled...\n', NSample)
Samples = [error_center,Samples];
Xbnd = [];

checkBnd = false;
falseCnt = 0;
checkInv = true;
for i =1:NSample
    X = Samples(:,i);
    disp(i)
    if V_full.eval(X)>lv_set_input
        f_obj = @(lambda) -lambda;
        constraint = @(v) deal([], V_full.eval(v*X)-lv_set_input);
        [lambda_star, ~] = fmincon(f_obj, 1e-4, [], [], [], [], 0, 1, constraint, opts);
        xb = lambda_star*X;
        Xbnd = [Xbnd,xb];
        checkBnd = V_full.grad(xb)*(f_field(xb)+g_field(xb)*ubar(xb))<=0;
        if ~checkBnd
            falseCnt = falseCnt +1;
            checkInv = false;
            warning('Invariance condition violated')
            break;
        end
    end
end
%% RG with SCLF
clc
close all
% some state constraints
sc =[-0,6;
    6,6];
coeff = [sc(1,:)', [1;1]]^(-1)*sc(2,:)';
a_HP = [-coeff(1), 1]'; b_HP = coeff(2);

V2D_1 = Lyapunov_quad([3 0.5;0.5,1.70]);
V2D_2 = Lyapunov_2p([3 0.5; 1 2; 0 0],2);

xt2s = [2.6; 2.7];
V_full_test = Lyapunov_para({V2D_1,V2D_2},1);
listlv = findLevelSetSCLF(V_full_test, a_HP, b_HP, xt2s);




% Grid for plotting level set
[x1, x2] = meshgrid(linspace(1, 5, 300), ...
    linspace(1, 6, 300));
zz = zeros(size(x1));

for i = 1:numel(x1)
    x = [x1(i); x2(i)];
    zz(i) = V2D_2.eval(x-xt2s);
end








figure
drawHyperPlane2D(a_HP,b_HP,sc(1,:))
hold on
Ellipsoid_plot2D(V2D_1.P,listlv(1),xt2s,'linestyle','--','color','b','linewidth',2)
contour(x1, x2, zz, [1 1]*listlv(2), 'r', 'LineWidth', 2);  % Level set
%% Avoidance scenario
ObsA = [3 3 1 0;
    2 6 6 2]';
ObsA = Polyhedron('V',ObsA);
[PA,cA]=Ellipse_outer(ObsA);

ObsB = [11 10 8 7;
    2 5 5 2]';
ObsB = Polyhedron('V',ObsB);
[PB,cB]=Ellipse_outer(ObsB);

Obs{1}.c = cA;Obs{1}.P = PA;
Obs{2}.c = cB;Obs{2}.P = PB;

%%
Waypoints=[0 4 4 10 12;
    0 2 4 6 3];
xcurr = [6;1];
xcurr = [2.16;1.08];

[a_safe,b_safe]=getSafeHyperPlane(xcurr,PA, cA) ;

clc
figure
axis square
hold on
Ellipsoid_plot2D(PA,1,cA,'linewidth',0.5,'color','k');
Ellipsoid_plot2D(PB,1,cB,'linewidth',0.5,'color','k');
plot(Waypoints(1,:),Waypoints(2,:),'s--','linewidth',1)
title('obstacles and path')

drawHyperPlane2D(a_safe,b_safe,[xcurr(1)-2 xcurr(1)+2]);
%% simulation
tau= 0.02;
h = tau/1;
tsim = 500;
Nsim = round(tsim/tau);
tt = linspace(0,tsim,Nsim);
RefLift = [1 0 0 0 0 0;
    0 0 1 0 0 0]';
X0 = [5;0;-2;0;-pi/12;0];
ref = zeros(2,Nsim);
reflist=[0 3.5 4 10 12;
    0.5 1 4 6  5];
switchReftime=floor(linspace(1,Nsim,size(reflist,2)+1));

for k = 1:numel(switchReftime)-1
    for i = switchReftime(k):switchReftime(k+1)
        ref(:,i) = reflist(:,k);
    end

end
Xsim= zeros(nx,Nsim+1);
Xsim(:,1)=X0;
usim=zeros(nu,Nsim);

% the linear constraint in R^6
a_constr = RefLift * a_HP;
b_constr = b_HP;

eta_erg = 0.2;
kappa = 3;
ref_virtual = X0([1 3]);

ref_list = zeros(2,Nsim);
clc
for i = 1:Nsim
    tic;
    xcurrent = Xsim(:,i);
    for j = 1:numel(Obs)
        % get the safe hyperplanes wrt the obstacles in Obs
        [a_safe, b_safe] = getSafeHyperPlane(xcurrent([1,3],:),Obs{j}.P, Obs{j}.c) ;
        % find the safest level wrt the state constraints
        % Do not forget to lift the space from ref space to state space
        lvl_list = findLevelSetSCLF(V_full, RefLift*a_safe, b_safe, RefLift*ref_virtual);
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
    %

    usim(:,i) =  ubar(Xsim(:,i)-RefLift*ref_virtual);
    if isnan(usim(1,i))
        error NAN
    end
    comptime(i) = toc;
    Xsim(:,i+1) = RK4(dyna,Xsim(:,i),usim(:,i)+ue,h,tau);
    ref_list(:,i)= ref_virtual;
end
%%
figure
hold on
plot(U,'alpha',0.05,'linewidth',0.5)
plot(usim(1,:),usim(2,:),'linewidth',3,'color','b')

figure
subplot(2,1,1)
plot(tt,usim(1,:))

subplot(2,1,2)
plot(tt,usim(2,:))
%%
darkgreen = [6, 140, 42]/255;
close all
idx = floor(linspace(1,Nsim,500));
figure
subplot(2,1,1)
plot(tt(idx),ref_list(1,idx),'--','linewidth',2,'color',darkgreen)
hold on
plot(tt(idx),Xsim(1,idx),'b-','linewidth',1)
plot(tt(idx),ref(1,idx),'r-.')
ylim([0 13])


subplot(2,1,2)
plot(tt(idx),ref(2,idx),'r-.')
hold on
plot(tt(idx),ref_list(2,idx),'--','linewidth',2,'color',darkgreen)
plot(tt(idx),Xsim(3,idx),'b-','linewidth',1)
ylim([-2.2 6.5])
legend('a','b','v')
%%
close all
d= 0.3;
col.patchalpha = 1;
col.patchcolor = [1 1 1]*100/255;
col.patchcolor = [110, 143, 156]/255;
col.linecolor = 'none';
st = floor(Nsim/300);
cnt_plot = 1;

figure
% drawHyperPlane2D(a_HP,b_HP,sc(1,:),'linestyle','-','color','k','linewidth',2)
hold on
while cnt_plot<=Nsim
    plot_vtol(Xsim(1,cnt_plot),Xsim(3,cnt_plot),Xsim(5,cnt_plot),d,col);
    cnt_plot = cnt_plot + st;
    % drawnow
    % pause(0.1)
end
plot(Xsim(1,1:Nsim),Xsim(3,1:Nsim),'linewidth',2)
% grid off
% ylim([1.5 6.5])
hold on
Ellipsoid_plot2D(PA,1,cA,'linewidth',0.5,'color','k');
Ellipsoid_plot2D(PB,1,cB,'linewidth',0.5,'color','k');
plot(reflist(1,:),reflist(2,:),'s--','linewidth',1)






%% snapshots
close all
d= 1.05;
col.patchalpha = 1;
col.patchcolor = [1 1 1]*100/255;
col.patchcolor = [110, 143, 156]/255;
col.linecolor = 'none';
drawBox = Polyhedron('ub', [20 20],'lb', [-5 -5]);
snapsIdx = floor([1 Nsim/3.5, Nsim/2, Nsim/1.35, Nsim]);
figure
for i = 1:numel(snapsIdx)
    idx = snapsIdx(i);
    subplot(1,numel(snapsIdx),i)
    hold on

    plot(reflist(1,:),reflist(2,:),'+--','linewidth',1,'color','r')

    plot(Xsim(1,1:idx),Xsim(3,1:idx),'linewidth',1,'color','b')

    Ellipsoid_patch2D(PA,1,cA,'linewidth',0.5,'edgecolor','k','facecolor','k','facealpha',0.5);
    Ellipsoid_patch2D(PB,1,cB,'linewidth',0.5,'edgecolor','k','facecolor','k','facealpha',0.5);
    xcurrent = Xsim(:,idx);

    for j = 1:numel(Obs)
        % get the safe hyperplanes wrt the obstacles in Obs
        [a_safe, b_safe] = getSafeHyperPlane(xcurrent([1,3],:),Obs{j}.P, Obs{j}.c) ;
        drawHyperPlane2D(a_safe,b_safe,[-5 20],'k-','linewidth',0.2);
        unSafeSet = Polyhedron('A', [drawBox.A;-a_safe'], 'b', [drawBox.b;-b_safe]);
        plot(unSafeSet, 'alpha',0.1,'edgealpha',0)
    end

    hold on
    plot_vtol(Xsim(1,idx),Xsim(3,idx),Xsim(5,idx),d,col);
    hold on
    xlim([-2 15])
    ylim([-3 7])
    axis square
end
%%



cnt = 0;

clear Fvid vid
vid = VideoWriter('Vtol_sclf_obs.mp4', 'MPEG-4');

idx_ani = floor(linspace(1,Nsim,150));
sizeFrame=[1920 1920]/1;
h=figure('units','pixels','Position',[50 100, sizeFrame ]);
for i = idx_ani
    cnt = cnt+1;

    hold on

    plot(reflist(1,:),reflist(2,:),'+--','linewidth',1,'color','r')

    plot(Xsim(1,1:i),Xsim(3,1:i),'linewidth',1,'color','b')

    Ellipsoid_patch2D(PA,1,cA,'linewidth',0.5,'edgecolor','k','facecolor','k','facealpha',0.5);
    Ellipsoid_patch2D(PB,1,cB,'linewidth',0.5,'edgecolor','k','facecolor','k','facealpha',0.5);
    xcurrent = Xsim(:,i);

    for j = 1:numel(Obs)
        % get the safe hyperplanes wrt the obstacles in Obs
        [a_safe, b_safe] = getSafeHyperPlane(xcurrent([1,3],:),Obs{j}.P, Obs{j}.c) ;
        drawHyperPlane2D(a_safe,b_safe,[-5 20],'k-','linewidth',0.2);
        unSafeSet = Polyhedron('A', [drawBox.A;-a_safe'], 'b', [drawBox.b;-b_safe]);
        plot(unSafeSet, 'alpha',0.1,'edgealpha',0)
    end

    hold on
    plot_vtol(Xsim(1,i),Xsim(3,i),Xsim(5,i),d,col);
    hold on
    xlim([-2 15])
    ylim([-3 7])
    axis square
    drawnow
    Fvid(cnt) = getframe(gcf);
    hold off
    clf
end

vid.FrameRate=10;
Ask2proceed('Save video? [Y]')
open(vid)
writeVideo(vid,Fvid);
close(vid)
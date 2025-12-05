
% ===============================================
% MPC on SO(3) using rotation vector formulation
% (CasADi + Opti in MATLAB)
% ===============================================
import casadi.*

% Parameters
J = diag([1.0, 1.2, 1.5]);
dt = 0.05;
N = 20;  % Horizon

% Declare optimization variables
opti = Opti();

X = opti.variable(6, N+1); % [phi; omega] state
U = opti.variable(3, N);   % torque

phi0 = opti.parameter(3,1);     % initial rotation vector error
omega0 = opti.parameter(3,1);   % initial angular velocity

% Initial state
x0 = [phi0; omega0];
opti.subject_to(X(:,1) == x0);

% Dynamics
for k = 1:N
    phi_k = X(1:3,k);
    omega_k = X(4:6,k);
    tau_k = U(:,k);

    % Dynamics: phi_dot = omega, omega_dot = J^{-1}(tau - omega x J omega)
    phi_next = phi_k + dt * omega_k;
    omega_cross = [0, -omega_k(3), omega_k(2);
                   omega_k(3), 0, -omega_k(1);
                  -omega_k(2), omega_k(1), 0];
    omega_dot = J \ (tau_k - omega_cross * J * omega_k);
    omega_next = omega_k + dt * omega_dot;

    opti.subject_to(X(:,k+1) == [phi_next; omega_next]);
end

% Cost function: penalize rotation error and torque
cost = 0;
for k = 1:N
    cost = cost + X(1:3,k)'*X(1:3,k) + 0.1 * U(:,k)'*U(:,k);
end
opti.minimize(cost);

% Input bounds
opti.subject_to(-1 <= U <= 1);
opts.ipopt.print_level=0;
opts.print_time=0;
opts.ipopt.sb='yes';
opts.ipopt.max_iter = 2000;
% Solver
opti.solver('ipopt',opts);

% Simulation
x_hist = zeros(6, 100);
phi_val = [0.4; 0.3; -0.2];
omega_val = [0.5; -0.3; 0.2];
x_hist(:,1) = [phi_val; omega_val];
torque = zeros(3, 100);
for t = 1:99
    opti.set_value(phi0, phi_val);
    opti.set_value(omega0, omega_val);

    try
        sol = opti.solve();
    catch
        warning('MPC failed at step %d', t);
        break;
    end

    % Apply control
    tau = sol.value(U(:,1));

    % Integrate forward (Euler)
    omega_cross = [0, -omega_val(3), omega_val(2);
                   omega_val(3), 0, -omega_val(1);
                  -omega_val(2), omega_val(1), 0];
    omega_dot = J \ (tau - omega_cross * J * omega_val);
    omega_val = omega_val + dt * omega_dot;
    phi_val = phi_val + dt * omega_val;

    x_hist(:,t+1) = [phi_val; omega_val];
    torque(:,t) = tau;
end

%% Plot results
time = (0:99)*dt;
figure;
plot(time, x_hist(1:3,:)');
xlabel('Time [s]');
ylabel('Rotation vector error');
legend('\phi_1','\phi_2','\phi_3');
title('Rotation Error under MPC');

figure;
plot(time, x_hist(4:6,:)');
xlabel('Time [s]');
ylabel('Angular velocity [rad/s]');
legend('\omega_1','\omega_2','\omega_3');
title('Angular Velocity under MPC');

figure
plot(time, torque(1,:))
hold on
plot(time, torque(2,:))
plot(time, torque(3,:))
%% SO(3) representation
R_log = zeros(3, 3, numel(x_hist(1,:)));

for i = 1:numel(x_hist(1,:))
    phi = x_hist(:,i);
    phihat = hat(phi);
    R_log(:,:,i) = expm(hat(phi));
end
%%
close all
figure;
view(3)
axis equal;
grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
xlim([-1.5 1.5]); ylim([-1.5 1.5]); zlim([-1.5 1.5]);
% title('Rotation Frame Evolution (SO(3))');
hold on;

% World axes
quiver3(0,0,0,1,0,0,'r','LineWidth',1);
quiver3(0,0,0,0,1,0,'g','LineWidth',1);
quiver3(0,0,0,0,0,1,'b','LineWidth',1);

for k = 1:5:100
    Rk = R_log(:,:,k);
    origin = [0;0;0];

    % Body axes (in world frame)
    Xb = Rk(:,1);
    Yb = Rk(:,2);
    Zb = Rk(:,3);

    % Clear previous frame
%     cla;
    quiver3(0,0,0,1,0,0,'r','LineWidth',1); % world frame
    quiver3(0,0,0,0,1,0,'g','LineWidth',1);
    quiver3(0,0,0,0,0,1,'b','LineWidth',1);

    % Plot body frame
    quiver3(0,0,0,Xb(1),Xb(2),Xb(3),'r','LineWidth',2);
    quiver3(0,0,0,Yb(1),Yb(2),Yb(3),'g','LineWidth',2);
    quiver3(0,0,0,Zb(1),Zb(2),Zb(3),'b','LineWidth',2);
    
    % Add text
%     text(1.1, 0, 0, 'World X');
%     text(0, 1.1, 0, 'World Y');
%     text(0, 0, 1.1, 'World Z');

%     title(sprintf('Time: %.2f s', time(k)));
%     drawnow;
end
Ellipsoid_plot3D(eye(3),1,[0;0;0],16,'facecolor','none','edgecolor', [1 1 1]*150/255);
%%
% ============================
% Helper Functions
% ============================
function S = hat(w)
    S = [  0   -w(3)  w(2);
          w(3)   0   -w(1);
         -w(2) w(1)    0 ];
end

function v = vee(S)
    if norm(S + S', 'fro') > 1e-8
        error('Input must be skew-symmetric');
    end
    v = [S(3,2); S(1,3); S(2,1)];
end
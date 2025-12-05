% ============================
% Attitude Tracking on SO(3)
% ============================
clear; clc;

% Parameters
J = diag([1.0, 1.2, 1.5]);   % Inertia matrix
dt = 0.01;                   % Time step
T = 5;                       % Total simulation time
N = round(T / dt);          % Number of steps

% Gains
kR = 8;     % Rotation gain
kW = 2.5;   % Angular velocity gain

% Desired constant rotation matrix Rd
Rd = eye(3);

% Initial conditions
R = expm(hat([0.5; 0.5; 0]));  % Initial rotation
omega = [0.5; -0.3; 0.2];      % Initial angular velocity

% Logging
R_log = zeros(3, 3, N);
omega_log = zeros(3, N);
time = zeros(1, N);

for k = 1:N
    t = (k-1)*dt;
    
    % Error rotation matrix
    Re = Rd' * R;
    
    % Log error (axis-angle)
    phi = vee(logm(Re));
    
    % Control torque
    tau = -kR * phi - kW * omega;
    
    % Dynamics
    omega_dot = J \ (tau - cross(omega, J*omega));
    omega = omega + dt * omega_dot;
    R = R * expm(dt * hat(omega));  % SO(3) integration using exponential map
    
    % Log
    R_log(:,:,k) = R;
    omega_log(:,k) = omega;
    time(k) = t;
end

% Plot angular velocity
figure;
plot(time, omega_log);
xlabel('Time [s]');
ylabel('Angular Velocity [rad/s]');
legend('\omega_1','\omega_2','\omega_3');
title('Body Angular Velocity');
%%

% ============================
% Animate Rotation Frame
% ============================
figure;
view(3)
axis equal;
grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
xlim([-1.5 1.5]); ylim([-1.5 1.5]); zlim([-1.5 1.5]);
title('Rotation Frame Evolution (SO(3))');
hold on;

% World axes
quiver3(0,0,0,1,0,0,'r','LineWidth',1);
quiver3(0,0,0,0,1,0,'g','LineWidth',1);
quiver3(0,0,0,0,0,1,'b','LineWidth',1);

for k = 1:10:N
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
    text(1.1, 0, 0, 'World X');
    text(0, 1.1, 0, 'World Y');
    text(0, 0, 1.1, 'World Z');

    title(sprintf('Time: %.2f s', time(k)));
    drawnow;
end

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

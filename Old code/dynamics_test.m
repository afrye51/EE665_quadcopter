% SS derivation from: https://www.kth.se/polopoly_fs/1.588039.1550155544!/Thesis%20KTH%20-%20Francesco%20Sabatino.pdf
% Some constants from: http://lup.lub.lu.se/luur/download?func=downloadFile&recordOId=8847641&fileOId=8859343

% x = [phi; = y_angle
%     theta;= z_angle
%     psi;  = x_angle
%     p;    = phi_dot
%     q;    = theta_dot
%     r;    = psi_dot
%     u;    = x_dot
%     v;    = y_dot
%     w;    = z_dot
%     x;
%     y;
%     z];

% Define constants
g = 9.8;            % m/s^2
m = 0.5;            % kg
Ix = 460 * 10^-6;   % N-m-s^2
Iy = 460 * 10^-6;   % N-m-s^2
Iz = 920 * 10^-6;   % N-m-s^2
motor_radius = 0.1; % m

%%
% Define state matrices
A_unstab = [0 0 0 1 0 0 0 0 0 0 0 0;...
    0 0 0 0 1 0 0 0 0 0 0 0;...
    0 0 0 0 0 1 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0 0 0;...
    0 -g 0 0 0 0 0 0 0 0 0 0;...
    g 0 0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 1 0 0 0 0 0;...
    0 0 0 0 0 0 0 1 0 0 0 0;...
    0 0 0 0 0 0 0 0 1 0 0 0];

B_ctrl = [0 0 0 0;...
    0 0 0 0;...
    0 0 0 0;...
    0 1/Ix 0 0;...
    0 0 1/Iy 0;...
    0 0 0 1/Iz;...
    0 0 0 0;...
    0 0 0 0;...
    1/m 0 0 0;...
    0 0 0 0;...
    0 0 0 0;...
    0 0 0 0];

B_unctrl = [0 0 0 0;...
    0 0 0 0;...
    0 0 0 0;...
    0 1/Ix 0 0;...
    0 0 1/Iy 0;...
    0 0 0 1/Iz;...
    0 0 0 0;...
    0 0 0 0;...
    1/m 0 0 0;...
    0 0 0 0;...
    0 0 0 0;...
    0 0 0 0];

C = eye(size(A));

% check observability and controllability
Obs = obsv(A_unstab, C);
rank(Obs)
Ctr = ctrb(A_unstab, B_ctrl);
rank(Ctr)
Ctr = ctrb(A_unstab, B_unctrl);
rank(Ctr)

% Initialization
dt = 0.1;
k = 1;
time = 10;
steps = round(time / dt);
x = zeros(size(A_unstab, 1), steps);
x_dot = zeros(size(A_unstab, 1), steps);
x(:,k) = [0.1; 0.1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 10];

x_uns = zeros(size(A_unstab, 1), steps);
x_dot_uns = zeros(size(A_unstab, 1), steps);
x_uns(:,k) = [0.1; 0.1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 10];

goal = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 10];

p_r = -2;
p_i1 = 0.5i;
p_i2 = 1i;
p = [p_r, p_r, p_r, p_r, p_r+p_i1, p_r-p_i1, p_r+p_i1, p_r-p_i1,...
    p_r+p_i2, p_r-p_i2, p_r+p_i2, p_r-p_i2];
k_ctrl = place(A_unstab, B_ctrl, p);
eig(A_unstab)
A_stab = (A_unstab - B_ctrl * k_ctrl);
eig(A_stab)


%%
% Thrusts from each motor, [FL, FR, BR, BL] in N
% Assuming motor directions are: [CW, CCW, CW, CCW] (From the top)
y_CW = [1 0 1 0];
y_CCW = [0 1 0 1];
r_right = [1 0 0 1];
r_left = [0 1 1 0];
p_fwd = [0 0 1 1];
p_back = [1 1 0 0];
u_rotors = [1, 1, 1, 1];

% Total motor force
f_m = sum(u_rotors);
% Torque about x (Roll)
t_x = sum(u_rotors.*p_fwd - u_rotors.*p_back);
% Torque about y (Pitch)
t_y = sum(u_rotors.*r_right - u_rotors.*r_left);
% Torque about z (Yaw)
t_z = sum(u_rotors.*y_CW - u_rotors.*y_CCW);

%%
u_val = 0.00001;
% [f_sum, f_y, -f_x, yaw]
u = [0; 0; u_val; 0];

% A is not invertible, so this method fails
%A_dt = expm(A*dt);
%B_dt = A^-1*(A_dt-eye(size(A)))*B;
%C_dt = C;

%%
for k = 1:steps
    
    u(:,k) = -k_ctrl * (x(:,k) - goal);
    x_dot(:,k) = A_unstab*x(:,k) + B_ctrl*u(:,k);
    x(:,k+1) = x_dot(:,k) * dt + x(:,k);
    
    x_dot_uns(:,k) = A_unstab*x_uns(:,k);
    x_uns(:,k+1) = x_dot_uns(:,k) * dt + x_uns(:,k);
end

%% Plotting
plt_x = x(10,:);
plt_y = x(11,:);
plt_z = x(12,:);

plt_x_uns = x_uns(10,:);
plt_y_uns = x_uns(11,:);
plt_z_uns = x_uns(12,:);

figure(1);

% Plot Obsv, Stab system
subplot(221);
plot3(plt_x, plt_y, plt_z);
title('Controllable, Stabilized system');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
grid();

% Plot Obsv, Stab system
subplot(222);
plot3(plt_x, plt_y, plt_z);
title('Controllable, Unstable system');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
grid();

% Plot Obsv, Stab system
subplot(223);
plot3(plt_x, plt_y, plt_z);
title('Uncontrollable, Stabilized system');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
grid();

% Plot Obsv, Stab system
subplot(224);
plot3(plt_x, plt_y, plt_z);
title('Uncontrollable, Unstable system');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
grid();
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
A = [0 0 0 1 0 0 0 0 0 0 0 0;...
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

B = [0 0 0 0;...
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
Obs = obsv(A, C);
rank(Obs)
Ctr = ctrb(A, B);
rank(Ctr)

% Initialization
dt = 0.1;
k = 1;
time = 10;
steps = round(time / dt);
x = zeros(size(A, 1), steps);
x_dot = zeros(size(A, 1), steps);
x(:,k) = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 10];

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
% Torque about x (Pitch)
t_x = sum(u_rotors.*p_fwd - u_rotors.*p_back);
% Torque about y (Roll)
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
    x_dot(:,k) = A*x(:,k);% + B*u;
    x(:,k+1) = x_dot(:,k) * dt + x(:,k);
end

plt_x = x(10,:);
plt_y = x(11,:);
plt_z = x(12,:);

figure(1);
plot3(plt_x, plt_y, plt_z);
title('Quadrotor position in 3-space');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
grid();
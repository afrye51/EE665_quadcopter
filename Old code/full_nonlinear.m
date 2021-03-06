% SS derivation from: https://www.kth.se/polopoly_fs/1.588039.1550155544!/Thesis%20KTH%20-%20Francesco%20Sabatino.pdf
% Some constants from: http://lup.lub.lu.se/luur/download?func=downloadFile&recordOId=8847641&fileOId=8859343
% Wind is not included at this time

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
dt = 0.001;
k = 1;
time = 20;
steps = round(time / dt);
u = zeros(4, steps);
e = zeros(size(A, 1), steps);
x = zeros(size(A, 1), steps);
lin_error = zeros(size(A, 1), steps);
x_dot = zeros(size(A, 1), steps);
x(:,k) = [1; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 10];

x_bar = zeros(size(A, 1), steps);
x_dot_bar = zeros(size(A, 1), steps);
x_bar(:,k) = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 10];

goal = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 10];

p_r = -2;
p_i1 = 0.5;
p_i2 = 1;
p = [p_r, p_r, p_r, p_r, p_r+p_i1, p_r-p_i1, p_r+p_i1, p_r-p_i1,...
    p_r+p_i2, p_r-p_i2, p_r+p_i2, p_r-p_i2];
k_ctrl = place(A, B, p);

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
u(:,k) = [9.9; 0; 0; 0];

% A is not invertible, so this method fails
%A_dt = expm(A*dt);
%B_dt = A^-1*(A_dt-eye(size(A)))*B;
%C_dt = C;

%%
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
for k = 1:steps
    % Nonlinear system state
    phi = x(1,k);
    theta = x(2,k);
    psi = x(3,k);
    p = x(4,k);
    q = x(5,k);
    r = x(6,k);
    u_ = x(7,k);
    v = x(8,k);
    w = x(9,k);
    
    F = u(1,k);
    Tx = u(2,k);
    Ty = u(3,k);
    Tz = u(4,k);
    
    Fwx = 0;
    Fwy = 0;
    Fwz = 0;
    Twx = 0;
    Twy = 0;
    Twz = 0;
    
    Ax_true = [p + r*cos(phi)*tan(theta) + q*sin(phi)*tan(theta);...
        q*cos(phi) - r*sin(phi);...
        r*cos(phi)/cos(theta) + q*sin(phi)/cos(theta);...
        (Iy-Iz)/Ix*r*q + (Tx+Twx)/Ix;...
        (Iz-Ix)/Iy*p*r + (Ty+Twy)/Iy;...
        (Ix-Iy)/Iz*p*q + (Tz+Twz)/Iz;...
        r*v - q*w - g*sin(theta) + Fwx/m;...
        p*w - r*u_ + g*sin(phi)*cos(theta) + Fwy/m;...
        q*u_ - p*v + g*cos(phi)*cos(theta) + (Fwz -F)/m;...
        w*(sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta))...
            - v*(cos(phi)*sin(psi) + cos(psi)*sin(phi)*sin(theta))...
            + u_*cos(psi)*cos(theta);...
        v*(cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))...
            - w*(sin(phi)*cos(psi) + sin(psi)*cos(phi)*sin(theta))...
            + u_*cos(psi)*cos(theta);...
        w*cos(psi)*cos(theta) - u_*sin(theta) + v*cos(theta)*sin(phi)];
    
    % Linear controller estimate
    %u(:,k) = -k_ctrl * (x(:,k) - goal);
    %x_dot_bar(:,k) = A*x(:,k) + B*u(:,k);
    %x_bar(:,k+1) = x_dot_bar(:,k) * dt + x_bar(:,k);
    
    e(:,k) = x(:,k) - goal;
    u(:,k) = -k_ctrl * e(:,k) + [g*m 0 0 0]';
    x_dot(:,k) = Ax_true + B*u(:,k);
    x(:,k+1) = x_dot(:,k) * dt + x(:,k);
    lin_error(:,k) = Ax_true - A*x(:,k);
end

plt_x = x(10,:);
plt_y = x(11,:);
plt_z = x(12,:);

%plt_x_bar = x_bar(10,:);
%plt_y_bar = x_bar(11,:);
%plt_z_bar = x_bar(12,:);

figure(1);
plot3(plt_x, plt_y, plt_z);
hold on
%plot3(plt_x_bar, plt_y_bar, plt_z_bar);
title('Quadrotor position in 3-space');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
%legend('True state', 'Estimated state');
grid();
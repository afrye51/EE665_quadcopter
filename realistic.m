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

C_unobs = [1 0 0 0 0 0 0 0 0 0 0 0;...
     0 1 0 0 0 0 0 0 0 0 0 0;...
     0 0 1 0 0 0 0 0 0 0 0 0;...
     0 0 0 1 0 0 0 0 0 0 0 0;...
     0 0 0 0 1 0 0 0 0 0 0 0;...
     0 0 0 0 0 1 0 0 0 0 0 0;...
     0 0 0 0 0 0 1 0 0 0 0 0;...
     0 0 0 0 0 0 0 1 0 0 0 0;...
     0 0 0 0 0 0 0 0 1 0 0 0];

 
B_true = [0 0 0 0;...
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

C_true = [1 0 0 0 0 0 0 0 0 0 0 0;...
     0 1 0 0 0 0 0 0 0 0 0 0;...
     0 0 1 0 0 0 0 0 0 0 0 0;...
     0 0 0 1 0 0 0 0 0 0 0 0;...
     0 0 0 0 1 0 0 0 0 0 0 0;...
     0 0 0 0 0 1 0 0 0 0 0 0;...
     0 0 0 0 0 0 1 0 0 0 0 0;...
     0 0 0 0 0 0 0 1 0 0 0 0;...
     0 0 0 0 0 0 0 0 1 0 0 0;...
     0 0 0 0 0 0 0 0 0 1 0 0;...
     0 0 0 0 0 0 0 0 0 0 1 0;...
     0 0 0 0 0 0 0 0 0 0 0 1];

C = eye(size(A_unstab));
% check observability and controllability
% Reference goal rank
fprintf('Controllable system\n');
fprintf('Number of states is: %d\n', size(A_unstab, 1));
fprintf('Rank of Controllability Matrix is: %d\n', rank(ctrb(A_unstab, B_ctrl)));
fprintf('Rank of Observability Matrix is: %d\n', rank(obsv(A_unstab, C)));

fprintf('\n\nUnobservable system\n');
fprintf('Number of states is: %d\n', size(A_unstab, 1));
fprintf('Rank of Controllability Matrix is: %d\n', rank(ctrb(A_unstab, B_ctrl)));
fprintf('Rank of Observability Matrix is: %d\n', rank(obsv(A_unstab, C_unobs)));

n_obs = rank(obsv(A_unstab, C_unobs));
Ctr = ctrb(A_unstab, B_ctrl);
Obs = obsv(A_unstab, C_unobs);
extra_rows =   [0 0 0; 0 0 0; 0 0 0;...
                0 0 0; 0 0 0; 0 0 0;...
                0 0 0; 0 0 0; 0 0 0;...
                1 0 0; 0 1 0; 0 0 1];
rank(extra_rows)
Q = horzcat(Obs(1:12,1:9), extra_rows);
rank(Q)
P = Q^-1;
A_reobs = P*A_unstab*(P^-1);
B_reobs = P*B_ctrl;
C_reobs = C_unobs*(P^-1);
A_reobs = A_reobs(1:n_obs,1:n_obs);
B_reobs = B_reobs(1:n_obs,:);
C_reobs = C_reobs(:,1:n_obs);

fprintf('\nKalman-decomposed system\n');
fprintf('Number of states is: %d\n', size(A_reobs, 1));
fprintf('Rank of Controllability Matrix is: %d\n', rank(ctrb(A_reobs, B_reobs)));
fprintf('Rank of Observability Matrix is: %d\n', rank(obsv(A_reobs, C_reobs)));

% Initialization
dt = 0.01;
k = 1;
time = 20;
steps = round(time / dt);

u = zeros(size(B_reobs, 2), steps);
e = zeros(size(A_reobs, 1), steps);

x = zeros(size(C_true, 1), steps);
x_dot = zeros(size(C_true, 1), steps);
y = zeros(size(C_true, 1), steps);
x(:,k) = [1; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1];

x_uns = zeros(size(A_reobs, 1), steps);
x_dot_uns = zeros(size(A_reobs, 1), steps);
y_uns = zeros(size(C_reobs, 1), steps);
x_uns(:,k) = [1; 1; 0; 0; 0; 0; 0; 0; 0];

x_bar = zeros(size(A_reobs, 1), steps);
x_dot_bar = zeros(size(A_reobs, 1), steps);
y_bar = zeros(size(C_reobs, 1), steps);
x_bar(:,k) = [1; 1; 0; 0; 0; 0; 0; 0; 0];

goal = [0; 0; 0; 0; 0; 0; 0; 0; 0];

p = [-2, -3, -4, -5, -6, -7, -8, -9, -10];
k_ctrl = place(A_reobs, B_reobs, p);
A_stab = (A_reobs - B_reobs * k_ctrl);

L_ctrl = place(A_reobs, C_reobs, p);
A_est = (A_stab - L_ctrl * C_reobs);

fprintf('Eigenvalues of systems\n');
fprintf('Unstabilized system Eigenvalues: ');
eig(A_reobs)
fprintf('\nStabilized system Eigenvalues: ');
eig(A_stab)
fprintf('\nEstimated system Eigenvalues: ');
eig(A_est)

% sys_cs = ss(A_stab, B_ctrl, C, 0);
% sys_cu = ss(A_unstab, B_ctrl, C, 0);
% 
% figure();
% step(sys_cs)
% figure();
% step(sys_cu)
% 
% figure();
% pzmap(sys_cs)
% figure();
% pzmap(sys_cu)

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
u(:,k) = [0; 0; 0; 0];

% A is not invertible, so this method fails
%A_dt = expm(A*dt);
%B_dt = A^-1*(A_dt-eye(size(A)))*B;
%C_dt = C;

%%
for k = 1:steps
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
    
    % Create control signal from current state estimate
    u(:,k) = -k_ctrl * (x_bar(:,k) - goal);
    e(:,k) = x_bar(:,k) - goal;
    u(:,k) = -k_ctrl * e(:,k);
    
    % True system
    y(:,k) = C_true * x(:,k);
    x_dot(:,k) = Ax_true + B_true*(u(:,k) - [g*m 0 0 0]');
    x(:,k+1) = x_dot(:,k) * dt + x(:,k);
    
    %Estimate
    u(:,k) = -k_ctrl * e(:,k);
    y_bar(:,k) = C_reobs * x_bar(:,k);
    x_dot_bar(:,k) = A_reobs*x_bar(:,k) + B_reobs*u(:,k) + L_ctrl*(y(1:9,k) - y_bar(:,k));
    x_bar(:,k+1) = x_dot_bar(:,k) * dt + x_bar(:,k);
    
    y_uns(:,k) = C_reobs * x_uns(:,k);
    x_dot_uns(:,k) = A_reobs*x_uns(:,k);
    x_uns(:,k+1) = x_dot_uns(:,k) * dt + x_uns(:,k);
end

%% Plotting
plt_x = y(7,:);
plt_y = y(8,:);
plt_z = y(9,:);

plt_x_bar = y_bar(7,:);
plt_y_bar = y_bar(8,:);
plt_z_bar = y_bar(9,:);

plt_x_uns = y_uns(7,:);
plt_y_uns = y_uns(8,:);
plt_z_uns = y_uns(9,:);

% Plot connection points (A point every n seconds on every plot for
% comparison)
points = [1, 11, 21, 41, 71, 101];


figure(1)
subplot(2, 2, 1);
% Plot Stab system
plot3(plt_x, plt_y, plt_z, 'linewidth', 2, 'DisplayName', 'True State');
title('Linear Velocity');
xlabel('x (m/s)');
ylabel('y (m/s)');
zlabel('z (m/s)');
grid();
hold on

plot3(plt_x_bar, plt_y_bar, plt_z_bar, 'linewidth', 2,...
    'DisplayName', 'Estimated State');
plot_markers(points, 0, plt_x, plt_y, plt_z);
plot_markers(points, 0, plt_x_bar, plt_y_bar, plt_z_bar);
legend();

subplot(2, 2, 2);
% Plot Stab system
plot3(y(10,:), y(11,:), y(12,:), 'linewidth', 2,...
    'DisplayName', 'Cartesian Location');
hold on
plot_markers(points, 0, y(10,:), y(11,:), y(12,:));
title('System position');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
grid();
legend();

subplot(2, 2, [3, 4]);
% Plot Stab system
title('Estimator Error vs. Time');
xlabel('time steps');
ylabel('error');
hold on
t = (0:steps - 1) * dt;
est_error = x(1:9,:) - x_bar;
for i = 1:9
    plot(t(1:101), est_error(i,1:101), 'linewidth', 2);
end
for i = 1:9
    plot_markers(points, 0, t(1:101), est_error(i,1:101));
end
grid();
legend('phi', 'theta', 'psi', 'p', 'q', 'r', 'u', 'v', 'w');

function plot_markers(namedPts, named, x, y, z)
num_named_pts = size(namedPts);
color = ['k', 'b', 'r', 'g', 'c', 'm', 'y'];
symbol = ['o', '+', '*', 'd', 's', 'p'];
if named == 1
    switch nargin
        case 4
            for j = 1:num_named_pts(2)
                plot(x(namedPts(j)), y(namedPts(j)),...
                    strcat(color(1),symbol(j)), 'linewidth', 2,...
                    'DisplayName' ,strcat('Point',  int2str(j)));
            end
        case 5
            for j = 1:num_named_pts(2)
                plot3(x(namedPts(j)), y(namedPts(j)), z(namedPts(j)),...
                    strcat(color(1),symbol(j)), 'linewidth', 2,...
                    'DisplayName', strcat('Point',  int2str(j)));
            end
        otherwise
            return
    end
else
    switch nargin
        case 4
            for j = 1:num_named_pts(2)
                plot(x(namedPts(j)), y(namedPts(j)),...
                    strcat(color(1),symbol(j)), 'linewidth', 2,...
                    'HandleVisibility','off');
            end
        case 5
            for j = 1:num_named_pts(2)
                plot3(x(namedPts(j)), y(namedPts(j)), z(namedPts(j)),...
                    strcat(color(1),symbol(j)), 'linewidth', 2,...
                    'HandleVisibility','off');
            end
        otherwise
            return
    end
end
end
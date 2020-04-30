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
dt = 0.1;
k = 1;
time = 10;
steps = round(time / dt);

x = zeros(size(A_reobs, 1), steps);
x_dot = zeros(size(A_reobs, 1), steps);
y = zeros(size(C_reobs, 1), steps);
x(:,k) = [1; 1; 0; 0; 0; 0; 0; 0; 0];

x_uns = zeros(size(A_reobs, 1), steps);
x_dot_uns = zeros(size(A_reobs, 1), steps);
y_uns = zeros(size(C_reobs, 1), steps);
x_uns(:,k) = [1; 1; 0; 0; 0; 0; 0; 0; 0];

x_bar = zeros(size(A_reobs, 1), steps);
x_dot_bar = zeros(size(A_reobs, 1), steps);
y_bar = zeros(size(C_reobs, 1), steps);
x_bar(:,k) = [1; 1; 0; 0; 0; 0; 10; 0; 0];

goal = [0; 0; 0; 0; 0; 0; 0; 0; 0];

p_r = -2;
p_i1 = 0.5i;
p_i2 = 1i;
p = [p_r, p_r, p_r, p_r+p_i1, p_r-p_i1, p_r+p_i1, p_r-p_i1,...
    p_r+p_i2, p_r-p_i2];
p = [-2, -3, -4, -5, -6, -7, -8, -9, -10];
k_ctrl = place(A_reobs, B_reobs, p);
A_stab = (A_reobs - B_reobs * k_ctrl);

p_r = -2;
p_i1 = 0.5i;
p_i2 = 1i;
p = [p_r, p_r, p_r, p_r+p_i1, p_r-p_i1, p_r+p_i1, p_r-p_i1,...
    p_r+p_i2, p_r-p_i2];
p = [-2, -3, -4, -5, -6, -7, -8, -9, -10];
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
u = [0; 0; u_val; 0];

% A is not invertible, so this method fails
%A_dt = expm(A*dt);
%B_dt = A^-1*(A_dt-eye(size(A)))*B;
%C_dt = C;

%%
for k = 1:steps
    
    % Create control signal from current state estimate
    u(:,k) = -k_ctrl * (x_bar(:,k) - goal);
    
    % True system
    y(:,k) = C_reobs * x(:,k);
    x_dot(:,k) = A_reobs*x(:,k) + B_reobs*u(:,k);
    x(:,k+1) = x_dot(:,k) * dt + x(:,k);
    
    %Estimate
    y_bar(:,k) = C_reobs * x_bar(:,k);
    x_dot_bar(:,k) = A_reobs*x_bar(:,k) + B_reobs*u(:,k) + L_ctrl*(y(:,k) - y_bar(:,k));
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

figure();

% Plot Stab system
plot3(plt_x, plt_y, plt_z);
title('Controllable, Stabilized system');
xlabel('x (m/s)');
ylabel('y (m/s)');
zlabel('z (m/s)');
grid();
hold on

plot3(plt_x_bar, plt_y_bar, plt_z_bar);
legend('True System', 'Estimated System');


close all
clear
clc

rng(0);

%% System Setup
sys_name = "double_integrator_noisy";
sys = get_system(sys_name);

%% Params and Buffers
T_ini = 2;
N = 5;
n = 2; 

Q = eye(sys.ny);
R = eye(sys.nu) * 0.1;

% penalty_type = "homogeneous";
penalty_type = "full";
lambda_g = 10;
lambda_y = 1000;

max_traj_len = 50;
tol = 1e-5;

%% Simulation Using Random Input Sequence (data for the Hankel matrix)
x_data = repmat(sys.xs, 1, T_ini+1);
u_data = repmat(sys.us, 1, T_ini);
y_data = repmat(sys.ys, 1, T_ini);
for i = 1:50
    u = (-1 + rand(sys.nu, 1) * 2) * 0.1;
    u_data = [u_data, u];
    y_data = [y_data, sys.model.h(x_data(:, end), u)];
    x_data = [x_data, sys.model.f(x_data(:, end), u)];
end

%% DeePC Controller Setup
fprintf("- Initializing DeePC Controllers ...\n");
data.u = u_data;
data.y = y_data;
params.Q = Q;
params.R = R;
params.T_ini = T_ini;
params.N = N;
params.n = n;
ldeepc_ctrl = DeePC(sys.constraints, sys.yf, sys.uf, params, data);
fprintf("\t Nominal DeePC initialized!\n");

params.lambda_g = lambda_g;
params.lambda_y = lambda_y;
rldeepc_ctrl = RobustDeePC(sys.constraints, sys.yf, sys.uf, params, data, penalty_type);
fprintf("\t Robust DeePC initialized!\n");

%% DeePC Simulation
fprintf("\n- Simulating DeePC (N = %d) ...\n", N);
x = repmat(sys.xs, 1, T_ini+1);
u = repmat(sys.us, 1, T_ini);
y = repmat(sys.ys, 1, T_ini);
    
for k = 1:max_traj_len
    [u_i, info] = ldeepc_ctrl.solve(u(:, end-T_ini+1:end), y(:, end-T_ini+1:end));
    if info.errorcode ~= 0
        error("Problem occurred: %s\n", yalmiperror(info.errorcode));
    end
    u = [u, u_i];
    y = [y, sys.model.h(x(:, end), u_i)];
    x = [x, sys.model.f(x(:, end), u_i)];

    if norm(y(:, end) - sys.yf, inf) < tol
        break;
    end
end
j = io_traj_cost(y, sys.yf, Q, u, sys.uf, R);
fprintf("\t trjectory cost = %.3f\n", j(T_ini+1));
fprintf("- Done!\n");

%% RDeePC Simulation
fprintf("\n- Simulating Robust DeePC (N = %d) ...\n", N);
x_r = repmat(sys.xs, 1, T_ini+1);
u_r = repmat(sys.us, 1, T_ini);
y_r = repmat(sys.ys, 1, T_ini);
pred_r = {};
    
for k = 1:max_traj_len
    [u_i, info] = rldeepc_ctrl.solve(u_r(:, end-T_ini+1:end), y_r(:, end-T_ini+1:end));
    if info.errorcode ~= 0
        error("Problem occurred: %s\n", yalmiperror(info.errorcode));
    end
    if mod(k, 5) == 1
        pred_r{length(pred_r)+1} = [y_r(:, end), info.Y];
    end
    u_r = [u_r, u_i];
    y_r = [y_r, sys.model.h(x_r(:, end), u_i)];
    x_r = [x_r, sys.model.f(x_r(:, end), u_i)];

    if norm(y_r(:, end) - sys.yf, inf) < tol
        break;
    end
end
j_r = io_traj_cost(y_r, sys.yf, Q, u_r, sys.uf, R);
fprintf("\t trjectory cost = %.3f\n", j_r(T_ini+1));
fprintf("- Done!\n");

%% Plot
cnstr_settings = {
    'LineWidth', 1.5, ...
    'Color', 'black'
};
traj_settings = {
    'LineWidth', 1, ...
    'Marker', '.', ...
    'MarkerSize', 8 ...
};
pred_traj_settings = {
    'LineWidth', 0.8, ...
    'Marker', '.', ...
    'MarkerSize', 6, ...
    'Color', '#999999'
};

figure(1);
% Constraints
y_cnstr = sys.y_max * [1, 1, -1, -1, 1; 1, -1, -1, 1, 1];
plot(y_cnstr(1,:), y_cnstr(2,:), cnstr_settings{:});
hold on;
% Trajectories
plot(y(1, :), y(2,:), traj_settings{:});
plot(y_r(1, :), y_r(2,:), traj_settings{:});
for i = 1:5
    y_pred = pred_r{i};
    plot(y_pred(1, :), y_pred(2,:), pred_traj_settings{:});
end
xlabel("y_1");
ylabel("y_2");
grid on;
legend(["Constraints", "DeePC", "RDeePC", "RDeePC-pred"]);
function sys = double_integrator_noisy()
    sys.name = "double-integrator-noisy";
    %% Noise Statistics
    n_mean = 0;
    n_std = 0.1;
    %% State-space Matrices
    sys.nx = 2;
    sys.nu = 1;
    sys.ny = 2;

    sys.model.A = [1, 1; 0, 1];
    sys.model.B = [0; 1];
    sys.model.C = [1, 0; 0, 1];
    sys.model.D = 0;
    sys.model.f = @(x, u) sys.model.A*x + sys.model.B*u;  % dynamics
    sys.model.h = @(x, u) sys.model.C*x + sys.model.D*u + (n_std * randn(sys.ny, 1) + n_mean);  % measurement
    
    %% Constraints
    sys.u_max = 1;
    sys.y_max = 4;
    % sys.constraints.Hx = [eye(sys.nx); -eye(sys.nx)];
    % sys.constraints.hx = ones(2*sys.nx, 1) * 4;
    sys.constraints.Hu = [eye(sys.nu); -eye(sys.nu)];
    sys.constraints.hu = ones(2*sys.nu, 1) * sys.u_max;
    sys.constraints.Hy = [eye(sys.ny); -eye(sys.ny)];
    sys.constraints.hy = ones(2*sys.ny, 1) * sys.y_max;
    
    %% Initial Condition
    % If the second state is chosen not to be 0, be careful constructing 
    % the initial trajectory
    sys.xs = [-3.50; 0];
    sys.us = 0;
    sys.ys = sys.model.h(sys.xs, sys.us);
    
    %% Equilibrium
    sys.xf = [0; 0];
    sys.uf = 0;
    sys.yf = 0;
end
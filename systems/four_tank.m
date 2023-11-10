function sys = four_tank()
    sys.name = "four-tank";
    %% State-space Matrices
    sys.model.A = [0.921,     0, 0.041,     0;
                       0, 0.918,     0, 0.033;
                       0,     0, 0.924,     0;
                       0,     0,     0, 0.937];
    sys.model.B = [0.017, 0.001;
                   0.001, 0.023;
                       0, 0.061;
                   0.072,     0];
    sys.model.C = [1, 0, 0, 0;
                   0, 1, 0, 0];
    sys.model.D = zeros(2, 2);
    sys.model.f = @(x, u) sys.model.A*x + sys.model.B*u;  % dynamics
    sys.model.h = @(x, u) sys.model.C*x + sys.model.D*u;  % measurement
    
    sys.nx = 4;
    sys.nu = 2;
    sys.ny = 2;
    
    %% Constraints
    sys.constraints.Ax = [eye(sys.nx); -eye(sys.nx)];
    sys.constraints.bx = [1.2; 1.2; 2; 2; 0; 0; 0; 0];
    sys.constraints.Au = [eye(sys.nu); -eye(sys.nu)];
    sys.constraints.bu = [1.2; 2; 0; 0];
    sys.constraints.Ay = [eye(sys.ny); -eye(sys.ny)];
    sys.constraints.by = [1.2; 1.2; 0; 0];
    
    %% Initial Condition
    map = [eye(sys.nx)-sys.model.A, -sys.model.B; sys.model.C, sys.model.D];
    sys.ys = [0.25; 0.25];
    temp = map \ [zeros(sys.nx, 1); sys.ys];
    sys.us = temp(end-sys.nu+1:end);
    sys.xs = temp(1:sys.nx);
    
    %% Equilibrium
    sys.yf = [1; 1];
    temp = map \ [zeros(sys.nx, 1); sys.yf];
    sys.uf = temp(end-sys.nu+1:end);
    sys.xf = temp(1:sys.nx);
end
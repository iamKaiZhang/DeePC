classdef RobustDeePC < handle
    properties
        nu  % input dimension
        ny  % output dimension

        % Equilibrium
        yf
        uf

        % Hankel matrices
        Up
        Yp
        Uf
        Yf

        ioconstraints
        params
        yalmip_optimizer
    end

    methods
        function obj = RobustDeePC(constraints, yf, uf, params, data, penalty)
            obj.yf = yf;
            obj.uf = uf;
            
            Q = params.Q;
            R = params.R;
            T_ini = params.T_ini;
            N = params.N;
            n = params.n;  % the guess of system's order
            obj.nu = size(R, 1);
            obj.ny = size(Q, 1);

            U = make_hankel(data.u, T_ini+N);
            Y = make_hankel(data.y, T_ini+N);
            % Persistently excitation check
            r = rank([U; Y], 1e-6);
            if r < obj.nu * (T_ini+N) + n
                error("\t PE condition is not satisfied!\n");
            end
            fprintf("\t PE condition is satisfied!\n");
            obj.Up = U(1:T_ini*obj.nu, :);
            obj.Uf = U(T_ini*obj.nu+1:end, :);
            obj.Yp = Y(1:T_ini*obj.ny, :);
            obj.Yf = Y(T_ini*obj.ny+1:end, :);

            obj.ioconstraints = constraints;
            obj.params = params;

            obj.init_optimizer(params.lambda_g, params.lambda_y, penalty);
        end

        function init_optimizer(obj, lambda_g, lambda_y, penalty)
            T_ini = obj.params.T_ini;
            N = obj.params.N;
            Q = obj.params.Q;
            R = obj.params.R;
            Hu = obj.ioconstraints.Hu;
            hu = obj.ioconstraints.hu;
            Hy = obj.ioconstraints.Hy;
            hy = obj.ioconstraints.hy;

            ng = size(obj.Up, 2);

            u_var = sdpvar(obj.nu, N, 'full');
            y_var = sdpvar(obj.ny, N, 'full');
            u_ini = sdpvar(obj.nu, T_ini, 'full');
            y_ini = sdpvar(obj.ny, T_ini, 'full');
            sigma = sdpvar(obj.ny * T_ini, 1, 'full');
            g = sdpvar(ng, 1, 'full');

            if penalty == "homogeneous"
                PI = pinv([obj.Up; obj.Yp; obj.Uf]) * [obj.Up; obj.Yp; obj.Uf];
                g_proj = (eye(ng) - PI) * g;
                objective = lambda_g * (g_proj'*g_proj);
            else
                objective = lambda_g * (g'*g);
            end
            objective = objective + lambda_y * (sigma'*sigma);
            constraints = [
                obj.Up * g == reshape(u_ini, [], 1), ...
                obj.Yp * g == reshape(y_ini, [], 1) + sigma, ...
                obj.Uf * g == reshape(u_var, [], 1), ...
                obj.Yf * g == reshape(y_var, [], 1), ...
            ];
            for i = 1:N
                objective = objective ...
                    + (y_var(:, i)-obj.yf)'*Q*(y_var(:, i)-obj.yf) ...
                    + (u_var(:, i)-obj.uf)'*R*(u_var(:, i)-obj.uf);
                constraints = [
                    constraints, ...
                    Hu * u_var(:, i) <= hu, ...
                    Hy * y_var(:, i) <= hy, ...
                ];
            end

            opts = sdpsettings('verbose', 1, 'solver', 'mosek', 'mosek.MSK_DPAR_INTPNT_QO_TOL_REL_GAP', 1e-8);
            obj.yalmip_optimizer = optimizer( ...
                constraints, objective, opts, ...
                {u_ini, y_ini}, ...  % params: initial condition
                {u_var, y_var, objective} ...  % outputs
            );
        end

        function [u, info] = solve(obj, u_ini, y_ini)
            [out, errorcode] = obj.yalmip_optimizer({u_ini, y_ini});
            [U, Y, objective] = out{:};
            u = U(:, 1);
            info = struct('errorcode', errorcode, 'objective', objective, 'U', U, 'Y', Y);
        end
    end

end
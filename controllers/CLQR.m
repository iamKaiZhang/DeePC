classdef CLQR < handle
    properties
        xf  % state equilibrium
        uf  % input equilibrium
        K
        project
    end

    methods
        function obj = CLQR(sys, Q, R)
            obj.xf = sys.xf;
            obj.uf = sys.uf;
            obj.K = -dlqr(sys.model.A, sys.model.B, Q, R);
            Au = sys.constraints.Au;
            bu = sys.constraints.bu;

            % optim-prob for constraints projection
            nu = size(Au, 2);
            u_var = sdpvar(nu, 1, 'full');
            u_raw = sdpvar(nu, 1, 'full');
            objective = (u_var - u_raw)' * (u_var - u_raw);
            constraints = (Au * u_var <= bu);

            opts = sdpsettings('verbose',1,'solver','mosek');
            obj.project = optimizer(constraints, objective, opts, u_raw, u_var);
        end

        function u = solve(obj, x)
            u = obj.uf + obj.project(obj.K * (x - obj.xf));
        end
    end
end

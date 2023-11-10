function j = io_traj_cost(y, yf, Q, u, uf, R)
% -------------------------------------------------------------------------
% This function computes the closed-loop cost associated to an input-state
% trajectory.
% 
% INPUTS
%   - y    [ny,*]   output trajectory
%   - u    [nu,*]   input trajectory
%   - Q    [ny,ny]  cost applied to output
%   - R    [nu,nu]  cost applied to input
%
% OUTPUTS
%   - j    [1,*]    cost associated to each input-output pair
% -------------------------------------------------------------------------
% Get number of input-output pairs in the trajectory
N = size(y, 2);
% Compute stage cost for each state and input
stageCost = zeros(1, N);
for i = 1:N
    stageCost(i) = (y(:,i)-yf)'*Q*(y(:,i)-yf)+(u(:,i)-uf)'*R*(u(:,i)-uf);
end
% Compute cumulated sum to get value function
j = cumsum(stageCost,'reverse');
end
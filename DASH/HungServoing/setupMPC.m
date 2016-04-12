function miqp = setupMPC(mpc)%Set up MPC Matricies

N = 10; % Horizon depth
miqp.N = N;
Q_seed = eye(N);
S_seed = zeros(N); 
T_seed = zeros(N,1);
for i=1:N
    S_seed = S_seed + tril(ones(N),-i+1);
    T_seed(i) = i;
end

Q_bar = cell2mat(arrayfun(@(n) mpc.Q*n, Q_seed, 'UniformOutput',false));
S_bar = cell2mat(arrayfun(@(n) mpc.A^abs(n-1)*mpc.B*(n~=0), S_seed, 'UniformOutput',false));
T_bar = cell2mat(arrayfun(@(n) mpc.A^n, T_seed, 'UniformOutput',false));

H = 2*S_bar' * Q_bar * S_bar;
F = 2*T_bar' * Q_bar * S_bar;

% Format as MILP
miqp.U_init  = zeros(N,1);
% v is the vector to minimize and is [U, z]', where U is [u(0) ... u(N-1)],
% where N is the horizon depth
miqp.U_ind = 1:N;   % Parts of v that are integers, all of U, integer constraints
miqp.f_J = @(x_k) [x_k'*F, 1]; % Modified cost function vector
miqp.A_ineq = @(u_bar_k) [u_bar_k'*H, -1]; % A*v <= b
miqp.b_ineq = @(u_bar_k) .5 * u_bar_k' * H * u_bar_k; % A*v <= b
miqp.lb = [-ones(1,N), 0]';  % Lower bound on [U_bar, z]
miqp.ub = [ones(1,N), Inf]'; % Upper bound on [U_bar, z]
miqp.dt = mpc.dt;

end
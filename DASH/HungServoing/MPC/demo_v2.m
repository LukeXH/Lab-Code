clear
close all

% Dynamic Properties, discrete form
I = 1;
c = .1;
dt= .1;
A = dt*[0, 1;
        0, -c/I] + eye(2);
B = dt*[0;1/I];

% Control Parameters
Q = diag([10,.1]);


%Set up MPC Matricies
N = 10; % Horizon depth
Q_seed = eye(N);
S_seed = zeros(N); 
T_seed = zeros(N,1);
for i=1:N
    S_seed = S_seed + tril(ones(N),-i+1);
    T_seed(i) = i;
end

Q_bar = cell2mat(arrayfun(@(n) Q*n, Q_seed, 'UniformOutput',false));
S_bar = cell2mat(arrayfun(@(n) A^abs(n-1)*B*(n~=0), S_seed, 'UniformOutput',false));
T_bar = cell2mat(arrayfun(@(n) A^n, T_seed, 'UniformOutput',false));

H = 2*S_bar' * Q_bar * S_bar;
F = 2*T_bar' * Q_bar * S_bar;

% Format as MILP
% v is the vector to minimize and is [U, z]', where U is [u(0) ... u(N-1)],
% where N is the horizon depth
U_ind = 1:N;   % Parts of v that are integers, all of U, integer constraints
f_J = @(x_k) [x_k'*F, 1]; % Modified cost function vector
A_ineq = @(u_bar_k) [u_bar_k'*H, -1]; % A*v <= b
b_ineq = @(u_bar_k) .5 * u_bar_k' * H * u_bar_k; % A*v <= b
lb = [-ones(1,N), 0]';  % Lower bound on [U_bar, z]
ub = [ones(1,N), Inf]'; % Upper bound on [U_bar, z]


% Set up MILP
x_init  = [.1, 0]';
U_init  = zeros(N,1);
options = optimoptions(@intlinprog,'Display','off'); % Suppress iterative display
% First run of MILP
A_mat = A_ineq(U_init);
b_mat = b_ineq(U_init);
[vLinVar,fval,exitFlagInt,output] = intlinprog(f_J(x_init),U_ind,...
                                        A_mat,b_mat,...
                                        [],[],...
                                        lb,ub,options);
thediff = 1e-4;
iter = 1; % iteration counter
U_star = vLinVar(U_ind); % the x variables
truequadratic = b_ineq(U_star);
z_slack = vLinVar(U_ind(end)+1:end); % slack variable value
history = [truequadratic,z_slack];

% Loop Through MILP to satisfy MIQP
while abs((z_slack - truequadratic)/truequadratic) > thediff % relative error
    A_mat = [A_mat;A_ineq(U_star)];
    b_mat = [b_mat;truequadratic];
    % Solve the problem with the new constraints
    [vLinVar,fval,exitFlagInt,output] = intlinprog(f_J(x_init),U_ind,...
                                        A_mat,b_mat,...
                                        [],[],...
                                        lb,ub,options);
%     U_star = (U_star+vLinVar(U_ind))/2; % Midway from the previous to the current
    U_star = vLinVar(U_ind); % the x variables
    truequadratic = b_ineq(U_star);%
    z_slack = vLinVar(U_ind(end)+1:end);
    history = [history;truequadratic,z_slack];
    iter = iter + 1;
end

X =S_bar*U_star + T_bar*x_init;
X = [X(1:2:end-1), X(2:2:end)];
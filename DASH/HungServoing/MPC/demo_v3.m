clear
close all

% Dynamic Properties, discrete form
I = 1;
c = .1;
dt= .2;
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
x_init  = [pi, 0]';
U_init  = zeros(N,1);
% v is the vector to minimize and is [U, z]', where U is [u(0) ... u(N-1)],
% where N is the horizon depth
U_ind = 1:N;   % Parts of v that are integers, all of U, integer constraints
f_J = @(x_k) [x_k'*F, 1]; % Modified cost function vector
A_ineq = @(u_bar_k) [u_bar_k'*H, -1]; % A*v <= b
b_ineq = @(u_bar_k) .5 * u_bar_k' * H * u_bar_k; % A*v <= b
lb = [-ones(1,N), 0]';  % Lower bound on [U_bar, z]
ub = [ones(1,N), Inf]'; % Upper bound on [U_bar, z]

options = optimoptions(@intlinprog,'Display','off'); % Suppress iterative display
[vLinVar,~,~,~] = solverMIQP(f_J(x_init),U_ind,[],[],[],[],lb,ub,A_ineq,b_ineq, U_init, options);
U_star = vLinVar(U_ind); % the x variables

X = A*x_init + B*U_star(1);
% X = S_bar*U_star + T_bar*x_init;
% X = [X(1:2:end-1), X(2:2:end)];

%%% Run a simulation
figure(100)
clf
Xs = [x_init'];
for i = 2:30
    [vLinVar,~,~,~] = solverMIQP(f_J(Xs(end,:)'),U_ind,[],[],[],[],lb,ub,A_ineq,b_ineq, U_init, options);
    U_star = vLinVar(U_ind); % the x variables

    X = A*Xs(end,:)' + B*U_star(1);
    Xs = [Xs; X'];
    
    Y = S_bar*U_star + T_bar*Xs(end,:)';
    plot(dt*(i:(i+N-1)),Y(1:2:end-1),'o')
    hold on
end
plot(dt*(1:i),Xs(:,1),'k', 'Linewidth',2)
xlabel('time (sec)', 'Fontsize', 14)
ylabel('x', 'Fontsize', 14)
title('MPC Response To Second Order System', 'Fontsize', 16)
plot(dt*(1:i), 0*(1:i),'r--','Linewidth', 1.5)
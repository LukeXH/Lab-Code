function [K, K_horizon] = MPC_K_gain(horizon, A, B, Q, R )
%MPC_K_GAIN Summary of this function goes here
%   Detailed explanation goes here

N = horizon; % Horizon depth
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

if ~exist('R')
H = 2*S_bar' * Q_bar * S_bar;
else
    R_seed = eye(N);
    R_bar = cell2mat(arrayfun(@(n) R*n, R_seed, 'UniformOutput',false));
    H = 2*(R_bar + S_bar' * Q_bar * S_bar);
end

F = 2*T_bar' * Q_bar * S_bar;

% Feedback Gain
K_horizon = -H\F';
K = K_horizon(1,:);
disp('Feedback Gain:')
disp(K)

end


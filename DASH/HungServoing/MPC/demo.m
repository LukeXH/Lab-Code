clear
close all

% Dynamic Properties
I = diag([.1,1]);
k = 1*[-1, 1;
        1,-1];
c = 1*[0, 0;
      -1, 0];
A = [zeros(2), eye(2);
     I\k, I\c];
B = [0,0,(I*[0;1])']';

% Control Parameters
Q = diag([.1,1,.1,1]);


%Set up MPC
N = 2; % Horizon depth
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

% Feedback Gain
K_horizon = -H\F';
K = K_horizon(1,:);
disp('Feedback Gain:')
disp(K)
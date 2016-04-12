% This simulates the same cube system as runCube_SO3, but without the constraint and with gravity.

% Store our path so we can reset it at the end,
% then add our dependencies to the path
old_path = path;
addpath('SFDyn', 'visualization', 'visualization/groupTheory', 'MPC')

% Begin model setup.
dtool = DynTool;

% Used variables and placeholders
[th_x,w_x] = dtool.addCoord('body_th_x');
[th_y,w_y] = dtool.addCoord('body_th_y');
[th_z,w_z] = dtool.addCoord('body_th_z');
[psi1,dpsi1] = dtool.addCoord('leg1_th_y'); % In neg y-axis
[psi2,dpsi2] = dtool.addCoord('leg2_th_y'); % In pos y-axis


% Parameters of cube
L_a = 1;      % Anchor bar length 
M   = 1;      % body Mass, kg
s   = .1;     % body Cube side length, m
m   = .025;    % Leg mass
L_l = .05;    % Leg length, m
I   = M*s^2/6;% Moment of inertia along cube face axis, symmetric
I   = diag([I,I,I]);
grav= 9.81;
c_z = 0.01;%.001;

% MPC properties
mpc.dt= .1;
mpc.A = mpc.dt*[0, 1;
                0, -c_z/I(1,1)] + eye(2);
mpc.B = mpc.dt*[0;c_z*2.6/I(1,1)];
% Control Parameters
mpc.Q = diag([1,.1]);

% Body position and velocity group
h0    = SE3( [0, 0, -s/2, 0, 0, 0] );   % From anchor on face to CoM
g_b   = SE3( [0, 0, 0, th_x, th_y, th_z] ) * h0;
g_c_a = [0, 0, 0, w_x, w_y, w_z].';    % Angular Velo of anchor pivot
g_c_b = h0.invAdj() * g_c_a;  % Angular Velo of body/CoM

% Right transform from body to leg 1, Assuming leg's pivot is on side of cube at center of mass
h1a     = SE3([0, -s/2 - L_l, 0, pi/4, 0, 0]);
h1b     = SE3([ L_l*sin(psi1), -L_l*cos(psi1), 0, 0, 0, psi1]);
h1      = h1a*h1b;
alpha1  = [0, 0, 0, 0, 0, dpsi1].';
g_L1    = g_b*h1;
g_c_L1  = simplify( h1.invAdj()*(g_c_b + alpha1) );

% Right transform from body to leg 2
h2a     = SE3([0, s/2 + L_l, 0, -pi/4, 0, 0]);
h2b     = SE3([ -L_l/2*sin(psi2), L_l/2*cos(psi2), 0, 0, 0, psi2]);
h2      = h2a*h2b;
alpha2  = [0, 0, 0, 0, 0, dpsi2].';
g_L2    = g_b*h2;
g_c_L2  = simplify( h2.invAdj()*(g_c_b + alpha2) );

% Energy
T_body = .5*M*(g_c_b(1:3).'*g_c_b(1:3)) + (g_c_b(4:6).' * I * g_c_b(4:6));
T_leg1 = .5*m*(g_c_L1(1:3).'*g_c_L1(1:3)) + (L_l^2*m) * (g_c_L1(6).'*g_c_L1(6));
T_leg2 = .5*m*(g_c_L2(1:3).'*g_c_L2(1:3)) + (L_l^2*m) * (g_c_L2(6).'*g_c_L2(6));
P_body = grav * (M*g_b.z + m*(g_L1.z + g_L2.z));
% P_body = sym(0);

% Add Energy terms to object
dtool.addKE(simplify(T_body));
dtool.addKE(simplify(T_leg1));
dtool.addKE(simplify(T_leg2));
dtool.addPE(simplify(P_body));
% accel_fcn = dtool.genAccelFcn;

sfdyn = dtool.genSFDyn();
accel_fcn = sfdyn.gen_accel_fcn();

%%% Set-up MPC
miqp = setupMPC(mpc);

%%% Simulate
[t,q] = simCube(accel_fcn, miqp, c_z);

%%% Display and Animate
animate

% Reset the path before exiting
path(old_path)

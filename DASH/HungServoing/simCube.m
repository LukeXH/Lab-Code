function [t,q] = simCube(accel_fcn, miqp, c)

tspan   = [0,20];
time    = tspan(1):miqp.dt:tspan(2);
q       = zeros(10*(length(time)-1),10);
t       = zeros(10*(length(time)-1),1);
U_init  = miqp.U_init;
q0      = zeros(10,1);
q0(3)   = pi/2;
q0(6:10)= [0;... % w0_x
           0;... % w0_y
           0;... % w0_z
           0;... % dpsi_leg1
           0];   % dpsi_leg2
x0      = [q0(3); q0(8)];

%%%  Set up options
ode_options = odeset('AbsTol',1e-9,'RelTol',1e-9);
mpc_options = optimoptions(@intlinprog,'Display','off'); % Suppress iterative display

%%% Run simulation
disp('Running Simulation...')
for i = 1:length(time)-1
    % Get Control input, trying to return to zero;
    e = x0 - [servotraj(time(i)/10);0];
    [vLinVar,~,~,~] = solverMIQP(miqp.f_J(e),miqp.U_ind,...
                            [],[],[],[],...
                            miqp.lb,miqp.ub,...
                            miqp.A_ineq,miqp.b_ineq,...
                            U_init,...
                            mpc_options);
    U_star = vLinVar(miqp.U_ind); % the x variables
    % Translate to desired velocities
    vel_d = 2*pi * [(1-2*(U_star(1)>0)); (1-2*(U_star(1)>-1))];
%     vel_d = 5*pi * [-U_star(1); -U_star(1)];
    % Put into Dynamics Function
    dyn_fcn = @(t, x) [x(6:10);...
                      accel_fcn(x(1:5), x(6:10), [zeros(3,1); velo_control(x(9:10), vel_d)]) ...
                      - diag([0,0,c,0,0])*x(6:10)];
                  
	% Run Simulation
%     sol = ode45(dyn_fcn, [time(i), time(i+1)]-time(i), q0, ode_options);
    sol = ode45(dyn_fcn, [0, miqp.dt], q0, ode_options);
    t(10*(i-1)+1:10*i)   = linspace(time(i), time(i+1), 10);
%     q(10*(i-1)+1:10*i,:) = deval( sol, linspace(time(i), time(i+1), 10) )';
    q(10*(i-1)+1:10*i,:) = deval( sol, linspace(0, miqp.dt, 10) )';
    
    % Variables set for next run
    U_init = U_star;
    q0   = q(10*i,:)';
    x0   = [q0(3); q0(8)];
    
end
disp('Ending Simulation...');

end
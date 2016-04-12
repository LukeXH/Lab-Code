% Script to animate the cube simulation.
% Must be run after runCube_ryan.

% Plot the z theta coordinate
figure(3)
subplot(2,1,1)
plot(t, q(:,3))
subplot(2,1,2)
plot(t, q(:,8))

figure(4)
clf
plot(t, q(:,3), 'LineWidth',2)
xlabel('time', 'Fontsize', 16)
ylabel('\theta_z', 'Fontsize', 18)
title('Change in Position Over Time', 'Fontsize', 18)

% Store our path so we can reset it at the end,
% then add our dependencies to the path
old_path = path;
addpath('../visualization', '../visualization/groupTheory')

%%% Animate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bodyColors = [[171 175 166]/300;
              [171 175 166]/300;
              [171 175 166]/300;
              [171 175 166]/300;
              [216 90 26]/230;
              [216 90 26]/230;
              [216 90 26]/230;
              [216 90 26]/230];


cube_side = 0.1; % m
leg_len   = 0.05; % m
body = CubeClass(cube_side);
leg1 = CubeClass([leg_len,.1*cube_side]);
leg2 = CubeClass([leg_len,.1*cube_side]);
body.colors = bodyColors;

% Create a figure 
disp('Setting up visualization...')
h.figure = figure(19000);
clf
plot3([0,0],[0,0],[0,1], 'k', 'Linewidth',2)    % Anchor Line
axis([-1, 1, -1, 1, -1, 1]/4)
xlabel('x')
ylabel('y')
zlabel('z')
axis square

% Put the shapes into a plot
body.plot
leg1.plot
leg2.plot

% The range over which
% t = linspace(tspan(1), tspan(2),800);
X = q';%deval(sol, t);

% Video output!
% vidout = VideoWriter('tmp.avi');
% open(vidout);


for i = 1:length(t)
    % Cube position
    body.resetFrame
    g_body = SE3( [zeros(1,3), X(1:3,i)'] ) * SE3( [0, 0, -cube_side/2, 0, 0, 0] );
    xyz(:,i) = g_body.xyz;
    body.globalMove( g_body );
    
    
    % leg1 position
    leg1.resetFrame
    psi1 = X(4,i); % Leg rotation, around z-axis offset from body frame in neg y-axis direction
    h1a = SE3([0, -(cube_side+.1*cube_side)/2 - leg_len, 0, pi/4, 0, 0]);
    h1b = SE3([ leg_len/2*sin(psi1), -leg_len/2*cos(psi1), 0, 0, 0, psi1+pi/2]);
    leg1.globalMove( g_body*h1a*h1b );
    
    % leg2 position
    leg2.resetFrame
    psi2 = X(5,i); % Leg rotation, around z-axis offset from body frame in pos y-axis direction
    h2a = SE3([0, (cube_side+.1*cube_side)/2 + leg_len, 0, -pi/4, 0, 0]);
    h2b = SE3([ -leg_len/2*sin(psi2), leg_len/2*cos(psi2), 0, 0, 0, psi2+pi/2]);
%     h2 = SE3([ leg_len/2*sin(psi2), (cube_side+.1*cube_side)/2 + leg_len - leg_len/2*cos(psi2), 0, 0, 0, psi2+pi/2]);
    leg2.globalMove( g_body*h2a*h2b);
    
    % Update data
    body.updatePlotData
    leg1.updatePlotData
    leg2.updatePlotData

    % Draw figure
    drawnow
    
    % Write this frame to the video file
%     frame = getframe;
%     writeVideo(vidout, frame);
end % End for Loop

% close(vidout);

% Reset the path before exiting
path(old_path)

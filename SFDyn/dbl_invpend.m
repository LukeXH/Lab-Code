% This simulates a double inverted pendulum, using the same coordinate system as the model in SymExtension

clear all
close all
clc

dtool = DynTool;

[phi1,dphi1] = dtool.addCoord('p1');
[phi2,dphi2] = dtool.addCoord('p2');

m1 = 1;
m2 = 1;
l1 = 1;
l2 = 1;
g = 9.81;

I1 = m1 * l1^2/4;
I2 = m2 * l2^2/4;

pos1 = l1/2 * [ sin(phi1), -cos(phi1) ];
pos2 = 2*pos1 + l2/2 * [ sin(phi2), cos(phi2) ];

vel1 = jacobian(pos1, [phi1; phi2]) * [dphi1; dphi2];
vel2 = jacobian(pos2, [phi1; phi2]) * [dphi1; dphi2];

dtool.addKE(simplify(m1/2 * dot(vel1, vel1) + I1/2 * dphi1^2));
dtool.addKE(simplify(m2/2 * dot(vel2, vel2) + I2/2 * dphi2^2));

dtool.addPE(simplify(g * m1 * pos1(2)));
dtool.addPE(simplify(g * m2 * pos2(2)));

sfdyn = dtool.genSFDyn;
accel_fcn = sfdyn.gen_accel_fcn;

ode_fcn = @(t, x) [x(3:4); accel_fcn(x(1:2), x(3:4), [0; 0])];

sol = ode45(ode_fcn, [0, 15], [0; pi/4; 0; 0]);

figure(11001)
clf
plotyy(sol.x,sol.y(1,:),sol.x,sol.y(3,:))
title('Link 1')
figure(11002)
clf
plotyy(sol.x,sol.y(2,:),sol.x,sol.y(4,:))
title('Link 2')

function mainTwoLink()
%THIS IS FOR HW 2 PROBLEMS 4.3 INVOLVING ACKERMAN STEERING CAR
% Written by Lucas Hill

% Set up symbols
syms L
q = sym('q',[3,1]);     % [x, y, theta]
dq= sym('dq',[3,1]);
a = sym('a',[1,1]);
da= sym('da',[1,1]);

g  = SE2( q );
% gh1 = g_1, from joint (center of car) to front link
h1 = SE2( [-L/2*cos(a(1)), L/2*sin(a(1)), -a(1)]);
% gh2 = g_2, from joint (center of car) to rear link
h2 = SE2( [ L/2*cos(a(1)), L/2*sin(a(1)), a(1)]);

%
invAdj = @(x,y) SE2Circ( simplify( inv(x.g)*y.g*x.g ) );

% Drag Matrix
C = diag([L/2, L, 2/3*(L/2)^3]);
% C = diag([2, 1, .1]);

% Force Vectors in g frame
g1_circ = invAdj(h1, SE2Circ(dq + [0;0;-da(1)]) );
g2_circ = invAdj(h2, SE2Circ(dq + [0;0;da(1)]) );

%%% LOW REYNOLDS SWIMMING %%%
% Set Constraints in vector form
F_a = [cos(a), sin(a), 0;...
      -sin(a), cos(a), 0;...
            0,   -L/2, 1] * C*g1_circ.param;
F_b = [cos(a),-sin(a), 0;...
       sin(a), cos(a), 0;...
            0,    L/2, 1] * C*g2_circ.param;
% raw_cons = simplify(C*g1_circ.param + C*g2_circ.param) == zeros(3,1);
raw_cons = simplify(F_a + F_b) == zeros(3,1);

% Pfaffian Constraints
P = equationsToMatrix(raw_cons, [dq.', da.']);
A_L = simplify(inv(P(1:3,1:3)))*P(1:3,end);
disp('Low Reynolds Swimming')
pretty(-A_L)


%%% HIGH REYNOLDS SWIMMING %%%
% Added Mass Matrix
syms w m rho% Width and mass of elipse joint liquid density

% Link Inertia Tensor
I = diag([m, m, m/5*(L^2 + w^2)/4]);
% Added Mass
M_added = diag(rho*[pi*w^2/4, pi*L^2/4, (L^2/4 - w^2/4)^2]);

% Kinetic Energy
KE = .5*g1_circ.param.'*(I+M_added)*g1_circ.param + .5*g2_circ.param.'*(I+M_added)*g2_circ.param;

F = jacobian( jacobian(KE, [dq;da]), [dq;da]);

A_H = simplify( F(1:3,1:3) \ F(1:3,4) );

disp('High Reynolds Swimming')
pretty(-A_H)


end % End Function
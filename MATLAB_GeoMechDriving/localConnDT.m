function local_connections = localConnDT()
%LOCALCONNDT Summary of this function goes here
%   Detailed explanation goes here

local_connections.Kiwi          = kiwiDrive();
local_connections.Differential   = diffDrive();
local_connections.AckermanRear  = ackermanDriveRear();
local_connections.AckermanFront = ackermanDriveFront();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LC = kiwiDrive()
%THIS IS AN EXTENSION:  FINDING THE LOCAL CONNECTION FOR A KIWIDRIVE
%(3 omni-wheels at 120deg spacing)
%     ===
%      x
%    y_|
% 
% \\       //
        

% Set up symbols
L = sym('L');
R = sym('R');
q = sym('q',[3,1]);     % [x, y, theta]
dq= sym('dq',[3,1]);
a = sym('a',[3,1]);
da= sym('da',[3,1]);

g  = SE2( q );
% Extrinsic transforms, parameterizations
% gh1 = g_1, from joint (center of car) to front wheel
h1 = [L, 0, 0];
% gh2 = g_2, from joint (center of car) to left back wheel
h2 = [-L/2, sqrt(3)*L/2, 2/3*pi ];
% gh3 = g_3, from joint (center of car) to right back wheel
h3 = [-L/2, -sqrt(3)*L/2, -2/3*pi ];

% Extrinsic inverse adjoint
invAd = @(h)[ cos(h(3)), sin(h(3)), h(1)*sin(h(3)) - h(2)*cos(h(3));...
             -sin(h(3)), cos(h(3)), h(1)*cos(h(3)) + h(2)*sin(h(3));...
                      0,         0,                              1];

% Right body velocity of g_1, g_2, g_3
g_1_r = invAd(h1)*dq;
g_2_r = invAd(h2)*dq;
g_3_r = invAd(h3)*dq;

% Set Constraints in vector form
raw_cons = [g_1_r(2) - R*da(1) == 0;...
            g_2_r(2) - R*da(2) == 0;...
            g_3_r(2) - R*da(3) == 0];

% Pfaffian Constraints
P = equationsToMatrix(raw_cons, [dq.', da.']);
A = simplify(inv(P(1:3,1:3)))*P(1:3,4:6);
disp('Kiwi Drive Connection')
pretty(-A)
LC = -A;

end  % End Function kiwiDrive

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LC = diffDrive()
%Returns the local connections for a 2-wheeled differential drive
% ||   x   ||
% || y_|   ||
% ||       ||

% Set up symbols
W = sym('W');
R = sym('R');
q = sym('q',[3,1]);     % [x, y, theta]
dq= sym('dq',[3,1]);
a = sym('a',[3,1]);
da= sym('da',[3,1]);

g  = SE2( q );
% Extrinsic transforms, parameterizations
% gh1 = g_1, from joint (center of car) to right wheel
h1 = [0, -W, 0];
% gh2 = g_2, from joint (center of car) to left wheel
h2 = [0, W, 0];

% Extrinsic inverse adjoint
invAd = @(h)[ cos(h(3)), sin(h(3)), h(1)*sin(h(3)) - h(2)*cos(h(3));...
             -sin(h(3)), cos(h(3)), h(1)*cos(h(3)) + h(2)*sin(h(3));...
                      0,         0,                              1];

% Right body velocity of g_1, g_2
g_1_r = invAd(h1)*dq;
g_2_r = invAd(h2)*dq;

% Set Constraints in vector form
raw_cons = [g_1_r(1) - R*da(1) == 0;...
            g_2_r(1) - R*da(2) == 0;...
            dq(2) == 0];

% Pfaffian Constraints
P = equationsToMatrix(raw_cons, [dq.', da.']);
A = simplify(inv(P(1:3,1:3)))*P(1:3,4:5);
disp('Differential Drive Connection')
pretty(-A)
LC = -A;

end  % End Function diffDrive

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LC = ackermanDriveRear()

% Set up symbols
L = sym('L');
R = sym('R');
q = sym('q',[3,1]);     % [x, y, theta]
dq= sym('dq',[3,1]);
a = sym('a',[2,1]);
da= sym('da',[2,1]);    % Rear Wheel speed, Front wheel turning speed

% gh = g_1, from rear wheel (center of car) to front wheel
h = SE2([L, 0, a(2)]);
invAdj_h = @(g)(inv(h.g)*g.g*h.g);

% Right body velocity of g
g_r = SE2Circ(dq.');
% Right body velocity of g_1, the front wheel
g_1_r = SE2Circ( simplify(invAdj_h(g_r)) );

%%% REAR WHEEL DRIVE %%%
% Set Constraints in vector form
raw_cons = [g_r.x - R*da(1) == 0;...
            g_r.y           == 0;...
            g_1_r.y         == 0];

% Pfaffian Constraints
P = equationsToMatrix(raw_cons, [dq.', da.']);
A = simplify(inv(P(1:3,1:3)))*P(1:3,4:5);
disp('Rear-Wheel Drive Ackerman Connection')
pretty(-A)

LC = -A;

end %End Ackerman Rear Driving

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LC = ackermanDriveFront()

% Set up symbols
L = sym('L');
R = sym('R');
q = sym('q',[3,1]);     % [x, y, theta]
dq= sym('dq',[3,1]);
a = sym('a',[2,1]);
da= sym('da',[2,1]);    % Front wheel speed and turning rate

% gh = g_1, from rear wheel (center of car) to front wheel
h = SE2([L, 0, a(2)]);
invAdj_h = @(g)(inv(h.g)*g.g*h.g);

% Right body velocity of g
g_r = SE2Circ(dq.');
% Right body velocity of g_1, the front wheel
g_1_r = SE2Circ( simplify(invAdj_h(g_r)) );        

%%% FRONT WHEEL DRIVE  %%%
% Set Constraints in vector form
raw_cons = [g_1_r.x - R*da(1) == 0;...
            g_r.y             == 0;...
            g_1_r.y           == 0];

% Pfaffian Constraints
P = equationsToMatrix(raw_cons, [dq.', da.']);
A = simplify(inv(P(1:3,1:3)))*P(1:3,4:5);
disp('Front-Wheel Drive Ackerman Connection')
pretty(-A)

LC = -A;

end %End Ackerman Front Driving

end


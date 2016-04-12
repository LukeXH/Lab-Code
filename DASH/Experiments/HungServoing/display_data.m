clear
load('Trial12016-01-27-06.26.52.mat')
quat = [RigidBodybugRotationW, RigidBodybugRotationX, RigidBodybugRotationY, RigidBodybugRotationZ];
xyz  = [RigidBodybugPositionX, RigidBodybugPositionY, RigidBodybugPositionZ];

%%% ALL DATA IS IN "WINDOW" MODE WHERE X IS HEIGHT, Y IS "RIGHT_LEFT" AND Z
%%% IS DEPTH.  THUS MUST ROTATE -90deg AROUND X AXIS.
% correct_trans = [0, 0,-1;
%                  0, 1, 0;
%                  1, 0, 0];
correct_trans = [1, 0, 0;
                 0, 0,-1;
                 0, 1, 0];

cxyz = (correct_trans*xyz')';

r_h = [1,0,0];  % Start Pointing in X Direction
R_h = repmat(r_h, size(quat,1), 1);

head = zeros(size(xyz));
angleX = zeros(size(xyz,1),1);
angleY = zeros(size(xyz,1),1);
angleZ = zeros(size(xyz,1),1);
gamma = zeros(size(xyz,1),1);   % x axis rotation
beta = zeros(size(xyz,1),1);    % y axis rotation
alpha = zeros(size(xyz,1),1);   % z axis rotation

for i = 1:size(quat,1)
    head(i,:) = quatrotate(quat(i,:),r_h);
    q = quaternion( quat(i,2), quat(i,3), quat(i,4), quat(i,1) );
    qRot = quaternion( 0, 0, 0, 1);     % rotate pitch 180 to avoid 180/-180 flip for nicer graphing
    q = mtimes(q, qRot);
    angles = EulerAngles(q,'xyz');
    angleX(i) = -angles(1) * 180.0 / pi;   % must invert due to 180 flip above
    angleY(i) = angles(2) * 180.0 / pi;
    angleZ(i) = -angles(3) * 180.0 / pi;   % must invert due to 180 flip above
    
    % Rotation Vector
    rot_mat = RotationMatrix(q);
    beta(i)    = asin(-rot_mat(3,1));
    alpha(i)   = acos(rot_mat(1,1)/cos(beta(i)));
    gamma(i)   = acos(rot_mat(3,3)/cos(beta(i)));
end
chead = (correct_trans*head')';

figure(1000)
clf
% http://www.mathworks.com/help/matlab/ref/quiverseries-properties.html
quiver3(cxyz(:,1),cxyz(:,2),cxyz(:,3), chead(:,1),chead(:,2),chead(:,3),...
        'MaxHeadSize',.005, 'LineWidth', 1)
hold on
plot3(cxyz(:,1),cxyz(:,2),cxyz(:,3))
axis square

figure(1001)
clf
plot(Time, unwrap( atan2(head(:,2),head(:,1)) ) )

figure(1002)
clf
plot(Time, head(:,3))

figure(1003)
clf
plot3(cxyz(:,1),cxyz(:,2),Time)

figure(1004)
clf
subplot(2,1,1)
nrm = sqrt( chead(:,1).^2 + chead(:,2).^2);
n = 1:10:size(Time,1);
quiver(Time(n), zeros(length(n),1), chead(n,1)./nrm(n),chead(n,2)./nrm(n),...
        'MaxHeadSize',.005, 'LineWidth', 1)
subplot(2,1,2)
plot(Time(n), atan2(chead(n,2)./nrm(n), chead(n,1)./nrm(n)))
    
figure(1005)
clf
subplot(4,1,1)
plot(Time, head(:,1))
hold on
plot(Time, head(:,2),'r')
plot(Time, head(:,3),'g')
subplot(4,1,2)
plot(Time, angleX)
hold on
plot(Time, angleY,'r')
plot(Time, angleZ,'g')
subplot(4,1,3)
plot(Time, gamma)
hold on
plot(Time, beta,'r')
plot(Time, alpha,'g')
subplot(4,1,4)
plot(Time, -atan2(chead(:,2), chead(:,1)))
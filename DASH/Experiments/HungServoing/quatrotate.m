function n = quatrotate( q, r )
%QUATROTATE Rotate a 3D vector r by quaternion q
%   Rotate a 1x3 vector, r, of the form [x,y,z] by the quarternion, q,
%   which is of the form [w,x,y,z].

n = [ 1-2*q(3)^2 - 2*q(4)^2,     2*(q(2)*q(3) + q(1)*q(4)), 2*(q(2)*q(4) - q(1)*q(3));
      2*(q(2)*q(3) - q(1)*q(4)), 1-2*q(2)^2 - 2*q(4)^2,     2*(q(3)*q(4) + q(1)*q(2));
      2*(q(2)*q(4) + q(1)*q(3)), 2*(q(3)*q(4) - q(1)*q(2)), 1-2*q(2)^2 - 2*q(3)^2;
    ]*r.';

end


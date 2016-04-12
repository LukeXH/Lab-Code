% This is a class for the Special Euclidean Group in three dimensions (SE3Circ).
% SE3Circ is a mathematical construct used for keeping track of position
% and orientation in 3D, and performing operations (translations,
% rotations, etc.) in this space.
%
% An SE3Circ object can be created by passing the constructor the following:
%   SE3Circ: Defaults to all zeros
%   SE3Circ(An SE3Circ object)
%   SE3Circ(A 4x4 matrix): Assumed to be SE3Circ
%   SE3Circ([x y z]): Velocities along x, y, and z axes with respect to the origin
%   SE3Circ([x y z gamma beta alpha]): The frames angular velocities about
%       x, y and z (dgamma, dbeta, dalpha) in that order w.r.t. the origin,
%       and then, velocitiesalong x, y, and z, w.r.t the origin.
%
% The SE3Circ object can be manipulated in the following ways:
%   Left multiplication will rotate by Rx, Ry, Rz, and then
%     translate by x, y, z with respect to the origin.
%   Right multiplication will translate by x, y, z, and then
%     rotate by Rz, Ry, Rx with respect to the body.
%   obj.x: Returns dx
%   obj.y: Returns dy
%   obj.z: Returns dz
%   obj.xyz: Returns [dx dy dz]
%   obj.r: Returns the angular velocities as a vector
%   obj.R: Returns the rotation as a matrix
%   obj.g: Returns the 4x4 SE3Circ matrix
%
% Author: Luke Hill

classdef SE3Circ < handle
    properties
        g   % SE3Circ element
    end

    methods
        % Constructor
        function obj = SE3Circ(argument)
            % Defaults if no argument
            if nargin == 0
                obj.g = zeros(4);
                return
            end

            if isa(argument,'SE3Circ') % Is of type SE3Circ
                obj.g = argument.g;
            elseif all([4 4] == size(argument)) % If a 4x4 matrix, assumed to be SE3Circ
                obj.g = argument;
            elseif length(argument) == 3
                % Vector of:
                % dx, dy, dz: velocities along x, y, and z

                % Unpack
                dx = argument(1);
                dy = argument(2);
                dz = argument(3);
                dth1 = 0;
                dth2 = 0;
                dth3 = 0;

                obj.g = [     0, -dth3,  dth2, dx;
                           dth3,     0, -dth1, dy;
                          -dth2,  dth1,     0, dz;
                              0,     0,     0,  0];
                          
            elseif length(argument) == 6
                % Vector of:
                % dx, dy, dz: velocities along x, y, and z
                % dth1, dth2, dth3: angular velocities about x, y, and z

                % Unpack
                dx = argument(1);
                dy = argument(2);
                dz = argument(3);
                dth1 = argument(4);
                dth2 = argument(5);
                dth3 = argument(6);

                obj.g = [     0, -dth3,  dth2, dx;
                           dth3,     0, -dth1, dy;
                          -dth2,  dth1,     0, dz;
                              0,     0,     0,  0];
            else
                class(argument)
                size(argument)
                error('Unknown argument to SE3Circ constructor of type and size listed above')
            end
        end % SE3Circ

        function product = mtimes(a,b) % Matrix multiplication
            % Left multiplication will rotate by Rx, Ry, Rz, and then
            %   translate by x, y, z with respect to the origin.
            % Right multiplication will translate by x, y, z, and then
            %   rotate by Rz, Ry, Rx with respect to the body.

            % If both are SE3Circ classes
            if (isa(a,'SE3Circ') && isa(b,'SE3Circ'))
                product = SE3Circ(a.g*b.g);
            elseif (isa(a,'SE3') && isa(b,'SE3Circ'))
                product = SE3Circ(a.g*b.g);
            elseif (isa(a,'SE3Circ') && isa(b,'SE3'))
                product = SE3Circ(a.g*b.g);
            % Else if a is SE3Circ and the other is a matrix
            elseif (isa(a,'SE3Circ') && isa(b,'double'))
                product = a.g*b;
            % Else if b is SE3Circ and the other is a matrix
            elseif (isa(a,'double') && isa(b,'SE3Circ'))
                product = a*b.g;
            else
                class(a)
                class(b)
                error('Undefined operation * for SE3Circ with types listed above')
            end
        end % mtimes

        function resultant = mrdivide(a,b) % Right matrix division
            % If both are SE3Circ classes
            if (isa(a,'SE3Circ') && isa(b,'SE3Circ'))
                resultant = SE3Circ(a.g/b.g);
            % Else if a is SE3Circ and the other is a matrix
            elseif (isa(a,'SE3Circ') && isa(b,'double'))
                resultant = a.g/b;
            % Else if b is SE3Circ and the other is a matrix
            elseif (isa(a,'double') && isa(b,'SE3Circ'))
                resultant = a/b.g;
            else
                class(a)
                class(b)
                error('Undefined operation / for SE3Circ with types listed above')
            end
        end % mrdivide

        function resultant = mldivide(a,b) % Left matrix division
            % If both are SE3Circ classes
            if (isa(a,'SE3Circ') && isa(b,'SE3Circ'))
                resultant = SE3Circ(a.g\b.g);
            % Else if a is SE3Circ and the other is a matrix
            elseif (isa(a,'SE3Circ') && isa(b,'double'))
                resultant = a.g\b;
            % Else if b is SE3Circ and the other is a matrix
            elseif (isa(a,'double') && isa(b,'SE3Circ'))
                resultant = a\b.g;
            else
                class(a)
                class(b)
                error('Undefined operation \ for SE3Circ with types listed above')
            end
        end % mldividie

        function xOut = x(obj)  % Extract x velocity
            xOut = obj.g(1,4);
        end

        function yOut = y(obj)  % Extract y velocity
            yOut = obj.g(2,4);
        end

        function zOut = z(obj)  % Extract z velocity
            zOut = obj.g(3,4);
        end

        function xyzOut = xyz(obj)  % Extract velocity Vector
            xyzOut = obj.g(1:3,4);
        end

        function rOut = r(obj) % Extract Rotation Vector
            % dth1,2,3: angular velocities rotations about x, y, and z
            % R(3,1) = -sin(beta)
            dth1 = obj.g(3,2);
            dth2 = obj.g(1,3);
            dth3 = obj.g(2,1);
            rOut = [dth1; dth2; dth3];
        end
        
        function rOut = param(obj) % Extract all velocities
           rOut = [obj.xyz; obj.r];
        end

        function ROut = R(obj) % Extract Rotation Matrix
            ROut = obj.g(1:3,1:3);
        end

    end % methods

end % classdef

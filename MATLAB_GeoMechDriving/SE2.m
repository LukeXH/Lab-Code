% This is a class for the Special Euclidean Group in two dimensions (SE2).
% SE2 is a mathematical construct used for keeping track of position
% and orientation in 2D, and performing operations (translations,
% rotations, etc.) in this space.
%
% An SE2 object can be created by passing the constructor the following:
%   SE2: Defaults to all zeros
%   SE2(An SE2 object)
%   SE2(A 3x3 matrix): Assumed to be SE2
%   SE2(A 2x2 matrix): Assumed to be SO2 (Special Orthogonal Group in 3D)
%   SE2([x y]): Translations along x and y axes with respect to the origin
%   SE2([x y theta]): The frame is rotated about z axis (theta)
%       beta, alpha) in that order with respect to the origin, and then
%       translated along x and y with respect to the origin.
%
% The SE2 object can be manipulated in the following ways:
%   Left multiplication will rotate by Rx, Ry, Rz, and then
%     translate by x, y, z with respect to the origin.
%   Right multiplication will translate by x, y, z, and then
%     rotate by Rz, Ry, Rx with respect to the body.
%   obj.x: Returns x
%   obj.y: Returns y
%   obj.xy: Returns [x y]
%   obj.r: Returns the rotation as theta
%   obj.R: Returns the rotation as a matrix
%   obj.g: Returns the 3x3 SE2 matrix
%
% Author: Luke Hill

classdef SE2 < handle
    properties
        g   % SE2 element
    end

    methods
        % Constructor
        function obj = SE2(argument)
            % Defaults if no argument
            if nargin == 0
                x = 0;
                y = 0;
                theta = 0;
                % Rotation matrix (SO3)
                R = [cos(theta), -sin(theta);...
                     sin(theta),  cos(theta)];

                obj.g = [[R; 0 0] [x; y; 1]];
                return
            end

            if isa(argument,'SE2') % Is of type SE2
                obj.g = argument.g;
            elseif all([3 3] == size(argument)) % If a 3x3 matrix, assumed to be SE2
                obj.g = argument;
            elseif all([2 2] == size(argument)) % If a 2x2 matrix, assumed to be SO2
                obj.g = [argument zeros(2,1);
                         0 0 1];
            elseif length(argument) == 2
                % Vector of:
                % x, y: translations along x, and y with respect to the origin

                % Unpack
                x = argument(1);
                y = argument(2);
                theta = 0;
                
                % Rotation matrix (SO3)
                R = [cos(theta), -sin(theta);...
                     sin(theta),  cos(theta)];

                obj.g = [[R; 0 0] [x; y; 1]];
            elseif length(argument) == 3
                % Vector of:
                % x, y: translations along x and ywith respect to the origin
                % theta: rotations about z with respect to the origin

                % Unpack
                x = argument(1);
                y = argument(2);
                theta = argument(3);
                
                % Rotation matrix (SO3)
                R = [cos(theta), -sin(theta);...
                     sin(theta),  cos(theta)];

                obj.g = [[R; 0 0] [x; y; 1]];
            else
                class(argument)
                size(argument)
                error('Unknown argument to SE2 constructor of type and size listed above')
            end
        end % SE2

        function product = mtimes(a,b) % Matrix multiplication
            % Left multiplication will rotate by Rx, Ry, Rz, and then
            %   translate by x, y, z with respect to the origin.
            % Right multiplication will translate by x, y, z, and then
            %   rotate by Rz, Ry, Rx with respect to the body.

            % If both are SE2 classes
            if (isa(a,'SE2') && isa(b,'SE2'))
                product = SE2(a.g*b.g);
            % Else if a is SE2 and the other is a matrix
            elseif (isa(a,'SE2') && isa(b,'double'))
                product = a.g*b;
            % Else if b is SE2 and the other is a matrix
            elseif (isa(a,'double') && isa(b,'SE2'))
                product = a*b.g;
            else
                class(a)
                class(b)
                error('Undefined operation * for SE2 with types listed above')
            end
        end % mtimes

        function resultant = mrdivide(a,b) % Right matrix division
            % If both are SE2 classes
            if (isa(a,'SE2') && isa(b,'SE2'))
                resultant = SE2(a.g/b.g);
            % Else if a is SE2 and the other is a matrix
            elseif (isa(a,'SE2') && isa(b,'double'))
                resultant = a.g/b;
            % Else if b is SE2 and the other is a matrix
            elseif (isa(a,'double') && isa(b,'SE2'))
                resultant = a/b.g;
            else
                class(a)
                class(b)
                error('Undefined operation / for SE2 with types listed above')
            end
        end % mrdivide

        function resultant = mldivide(a,b) % Left matrix division
            % If both are SE2 classes
            if (isa(a,'SE2') && isa(b,'SE2'))
                resultant = SE2(a.g\b.g);
            % Else if a is SE2 and the other is a matrix
            elseif (isa(a,'SE2') && isa(b,'double'))
                resultant = a.g\b;
            % Else if b is SE2 and the other is a matrix
            elseif (isa(a,'double') && isa(b,'SE2'))
                resultant = a\b.g;
            else
                class(a)
                class(b)
                error('Undefined operation \ for SE2 with types listed above')
            end
        end % mldividie

        function xOut = x(obj)  % Extract x Position
            xOut = obj.g(1,3);
        end

        function yOut = y(obj)  % Extract y Position
            yOut = obj.g(2,3);
        end
        
        function xyOut = xy(obj)  % Extract Position Vector
            xyOut = obj.g(1:2,3);
        end

        function rOut = r(obj) % Extract Rotation Vector
            % gamma, beta, alpha: rotations about x, y, and z with respect to the origin
            % R(3,1) = -sin(beta)
            rOut = -asin(obj.g(1,2));
        end

        function ROut = R(obj) % Extract Rotation Matrix
            ROut = obj.g(1:2,1:2);
        end

        function pOut = param(obj) % Extract Parameters
            pOut = [obj.xy; obj.r];
        end
        
    end % methods

end % classdef

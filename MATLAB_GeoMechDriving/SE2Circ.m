% This is a class for the Special Euclidean Group in two dimensions (SE2Circ).
% SE2Circ is a mathematical construct used for keeping track of position
% and orientation in 2D, and performing operations (translations,
% rotations, etc.) in this space.
%
% An SE2Circ object can be created by passing the constructor the following:
%   SE2Circ: Defaults to all zeros
%   SE2Circ(An SE2Circ object)
%   SE2Circ(A 3x3 matrix): Assumed to be SE2Circ
%   SE2Circ(A 2x2 matrix): Assumed to be SO2 (Special Orthogonal Group in 3D)
%   SE2Circ([x y]): Translations along x and y axes with respect to the origin
%   SE2Circ([x y theta]): The frame is rotated about z axis (theta)
%       beta, alpha) in that order with respect to the origin, and then
%       translated along x and y with respect to the origin.
%
% The SE2Circ object can be manipulated in the following ways:
%   Left multiplication will rotate by Rx, Ry, Rz, and then
%     translate by x, y, z with respect to the origin.
%   Right multiplication will translate by x, y, z, and then
%     rotate by Rz, Ry, Rx with respect to the body.
%   obj.x: Returns x
%   obj.y: Returns y
%   obj.xy: Returns [x y]
%   obj.r: Returns the rotation as theta
%   obj.R: Returns the rotation as a matrix
%   obj.g: Returns the 3x3 SE2Circ matrix
%
% Author: Luke Hill

classdef SE2Circ < handle
    properties
        g   % SE2CircCirc element
    end

    methods
        % Constructor
        function obj = SE2Circ(argument)
            % Defaults if no argument
            if nargin == 0
                dx = 0;
                dy = 0;
                dtheta = 0;
                obj.g = [0, -dtheta, dx;...
                         dtheta, 0 , dy;...
                         zeros(1,3)];
                return
            end

            if isa(argument,'SE2Circ') % Is of type SE2Circ
                obj.g = argument.g;
            elseif all([3 3] == size(argument)) % If a 3x3 matrix, assumed to be SE2Circ
                obj.g = argument;
            elseif all([2 2] == size(argument)) % If a 2x2 matrix, assumed to be SO2
                obj.g = [argument zeros(2,1);
                         0 0 0];
            elseif length(argument) == 2
                % Vector of:
                % x, y: translations along x, and y with respect to the origin

                % Unpack
                dx = argument(1);
                dy = argument(2);
                dtheta = 0;
                obj.g = [0, -dtheta, dx;...
                         dtheta, 0 , dy;...
                         zeros(1,3)];
            elseif length(argument) == 3
                % Vector of:
                % x, y: translations along x and ywith respect to the origin
                % theta: rotations about z with respect to the origin

                % Unpack
                dx = argument(1);
                dy = argument(2);
                dtheta = argument(3);
                
                obj.g = [0, -dtheta, dx;...
                         dtheta, 0 , dy;...
                         zeros(1,3)];
            else
                class(argument)
                size(argument)
                error('Unknown argument to SE2Circ constructor of type and size listed above')
            end
        end % SE2Circ

        function product = mtimes(a,b) % Matrix multiplication
            % Left multiplication will rotate by Rx, Ry, Rz, and then
            %   translate by x, y, z with respect to the origin.
            % Right multiplication will translate by x, y, z, and then
            %   rotate by Rz, Ry, Rx with respect to the body.

            % If both are SE2Circ classes
            if (isa(a,'SE2Circ') && isa(b,'SE2Circ'))
                product = SE2Circ(a.g*b.g);
            % Else if a is SE2Circ and the other is a matrix
            elseif (isa(a,'SE2Circ') && isa(b,'double'))
                product = a.g*b;
            % Else if b is SE2Circ and the other is a matrix
            elseif (isa(a,'double') && isa(b,'SE2Circ'))
                product = a*b.g;
            else
                class(a)
                class(b)
                error('Undefined operation * for SE2Circ with types listed above')
            end
        end % mtimes

        function resultant = mrdivide(a,b) % Right matrix division
            % If both are SE2Circ classes
            if (isa(a,'SE2Circ') && isa(b,'SE2Circ'))
                resultant = SE2Circ(a.g/b.g);
            % Else if a is SE2Circ and the other is a matrix
            elseif (isa(a,'SE2Circ') && isa(b,'double'))
                resultant = a.g/b;
            % Else if b is SE2Circ and the other is a matrix
            elseif (isa(a,'double') && isa(b,'SE2Circ'))
                resultant = a/b.g;
            else
                class(a)
                class(b)
                error('Undefined operation / for SE2Circ with types listed above')
            end
        end % mrdivide

        function resultant = mldivide(a,b) % Left matrix division
            % If both are SE2Circ classes
            if (isa(a,'SE2Circ') && isa(b,'SE2Circ'))
                resultant = SE2Circ(a.g\b.g);
            % Else if a is SE2Circ and the other is a matrix
            elseif (isa(a,'SE2Circ') && isa(b,'double'))
                resultant = a.g\b;
            % Else if b is SE2Circ and the other is a matrix
            elseif (isa(a,'double') && isa(b,'SE2Circ'))
                resultant = a\b.g;
            else
                class(a)
                class(b)
                error('Undefined operation \ for SE2Circ with types listed above')
            end
        end % mldividie
        
        function resultant = plus(a,b) % Add matrix
            % If both are SE2Circ classes
            if (isa(a,'SE2Circ') && isa(b,'SE2Circ'))
                resultant = SE2Circ(a.g+b.g);
            % Else if a is SE2Circ and the other is a matrix
            elseif (isa(a,'SE2Circ') && isa(b,'double'))
                resultant = a.g+b;
            % Else if b is SE2Circ and the other is a matrix
            elseif (isa(a,'double') && isa(b,'SE2Circ'))
                resultant = a+b.g;
            else
                class(a)
                class(b)
                error('Undefined operation + for SE2Circ with types listed above')
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
            % R(2,1) = sin(beta)
            rOut = obj.g(2,1);
        end

        function ROut = R(obj) % Extract Rotation Matrix
            ROut = obj.g(1:2,1:2);
        end
        
        function pOut = param(obj) % Extract Parameters
            pOut = [obj.xy; obj.r];
        end

    end % methods

end % classdef

% Convenience class for assembling SFDyn objects.
% Allows for easy coordinate creation, energy specification,
% and creates the SFDyn object for the user.

classdef DynTool < handle
	methods (Access = public)
		% The constructor for this class; does some basic initialization
		%
		% Parameters: none
		%
		% Returns:
		%     this A handle to this object
		function this = DynTool
			this.lagrangian = sym(0);
            this.energy = sym(0);
		end % LagrDyn

		% This function adds a new coordinate to this dynamics system.
		%
		% Parameters:
		%     name The symbolic variable name for this parameter
		%
		% Returns:
		%     expr  The symbolic expression for this coordinate
		%     dexpr The symbolic expression for this coordinate's derivative
		function [expr,dexpr] = addCoord(this, name)
			% Error checking
			if nargin < 2
				error('addCoord needs a name parameter.')
			end
			if ~ischar(name)
				error('Input "name" must be a string')
			end

			% Diagnostic message for the user
			disp(['Adding coordinate ' name])

			% Generate symbolic expressions
			expr   = sym(name,        'real');
			dexpr  = sym(['d'  name], 'real');

			% Append this coordinate to the coordinates list
			this.coords(end+1,1)   = expr;
			this.dcoords(end+1,1)  = dexpr;
		end % addCoord

		% Add a new kinetic energy term
		%
		% Parameters:
		%     expr Kinetic energy expression
		%
		% Returns:
		%     expr The kinetic energy expression
		function expr = addKE(this, expr)
			% Diagnostic message for the user
			disp(['Adding kinetic energy expression ' char(expr)])

			% Directly add the expression to the Lagrangian
			this.lagrangian = this.lagrangian + expr;
            
            % Directly add the expression to the Energy Equation too
            this.energy = this.energy + expr;
		end % addKE

		% Add a new potential energy term
		%
		% Parameters:
		%     expr Potential energy expression
		%
		% Returns:
		%     expr The potential energy expression
		function expr = addPE(this, expr)
			% Diagnostic message for the user
			disp(['Adding potential energy expression ' char(expr)])

			% Directly add (well, subtract) the expression from the Lagrangian
			this.lagrangian = this.lagrangian - expr;
            
            % Directly add the expression to the Energy Equation too
            this.energy = this.energy + expr;
		end % addPE
        
        % Add a velocity constraint, A(q)dq = 0
		%
		% Parameters:
		%     expr Velocity constraint term
		%
		% Returns:
		%     expr The velocity constraint 
        function expr = addVeloCon(this, expr)
			% Diagnostic message for the user
			disp(['Adding velocity constraint expression ' char(expr)])

			% Directly add (well, subtract) the expression from the Lagrangian
			this.constraints = expr;
		end % addPE

        %%% Function Generation Section
        % Generates Energy Function
        function nrgfcn = getEnergy(self)
            % Pass out the Lagrangian as a function handle
            nrgfcn = matlabFunction(self.energy, 'vars', {self.coords, self.dcoords});
        end
        
        % Generates Constraint Function
        function confcn = getConstraints(self)
            % Pass out the Lagrangian as a function handle
            confcn = matlabFunction(self.constraints, 'vars', {self.coords});
        end
        
		% Generates and returns the SFDyn object
		function sfdyn = genSFDyn(self)
			% Pass our stored values to SFDyn's constructor, which does all the work
			sfdyn = SFDyn(self.lagrangian, self.coords, self.dcoords);
        end
        
        % Generates and returns the SFDyn object
		function sfdyn = genSFDynCons(self)
			% Pass our stored values to SFDyn's constructor, which does all the work
			sfdyn = SFDyn(self.lagrangian, self.coords, self.dcoords, self.constraints);
		end
	end % methods

	properties (GetAccess = public, SetAccess = private)
		% Coordinates and their derivatives
		coords@sym
		dcoords@sym

		% Lagrangian expression
		lagrangian@sym
        
        % Constraint expression
        constraints@sym
        
        % Energy expression
        energy@sym
	end % properties
end % classdef

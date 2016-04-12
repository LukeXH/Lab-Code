function [v,fval,exitFlag,msgs] = solverMIQP(f,intcon,A,b,Aeq,beq,lb,ub,A_slack,b_slack, v_init, options)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% First run of MILP
Ai = A_slack(v_init);
bi = b_slack(v_init);
dim1_A = size(Ai,1);
dim2_A = size(Ai,2);
dim2_b = size(bi,2);
A_mat = zeros( 99*dim1_A, dim2_A);
b_mat = zeros( 99*dim1_A, dim2_b);
A_mat = [A; Ai; A_mat];
b_mat = [b; bi; b_mat];
offset = size(A,1);
[v,fval,exitFlag,msgs] = intlinprog(f,intcon,...
                                    A_mat(1:offset+dim1_A,:),...
                                    b_mat(1:offset+dim1_A,:),...
                                    Aeq,beq,lb,ub,options);

thediff = 1e-4;
v_star = v(intcon); % the integer variables
truequadratic = b_slack(v_star);
z_slack = v(intcon(end)+1:end); % slack variable value
% history = [truequadratic,z_slack];

% Loop Through MILP to satisfy MIQP
for i = 2:100
    if abs((z_slack - truequadratic)/truequadratic) < thediff % relative error
       break 
    end
    A_mat(offset+dim1_A*i ,:) = A_slack(v_star);
    b_mat(offset+dim1_A*i ,:) = truequadratic;
    % Solve the problem with the new constraints
    [v,fval,exitFlag,msgs] = intlinprog(f,intcon,...
                                        A_mat(1:offset+dim1_A*i,:),...
                                        b_mat(1:offset+dim1_A*i,:),...
                                        Aeq,beq,lb,ub,options);
%     U_star = (U_star+vLinVar(intcon))/2; % Midway from the previous to the current
    v_star = v(intcon); % the x variables
    truequadratic = b_slack(v_star);%
    z_slack = v(intcon(end)+1:end);
%     history = [history;truequadratic,z_slack];
end

end


function dX = dynFun(X, u)
% Ta-daa
dX = zeros(4,1);

% Dynamic Properties
I = diag(.1,1);
k = 1*[ 1, -1;
       -1,  1];
c = 1*[0, 0;
       1, 0];

A = [zeros(2), eye(2);
     I\k, I\c];
 
B = [0,0,0,1]';

dX = A*X + B*u;

end
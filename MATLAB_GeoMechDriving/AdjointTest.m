clear
clc

syms x y th dx dy w

h     = SE2([x,y,th]);
g_c_1 = SE2Circ([dx, dy, w]);

g_c_2a = SE2Circ( simplify(inv(h.g)*g_c_1.g*h.g) );
disp('Long way')
pretty(g_c_2a.param)

g_c_2b = simplify( h.g*g_c_1.param );
disp('Short way? (Wrong Way?)')
pretty(g_c_2b)

Jac = jacobian(g_c_2a.param, g_c_1.param);
disp('Adjoint Matrix, for h-inverse')
pretty(Jac)
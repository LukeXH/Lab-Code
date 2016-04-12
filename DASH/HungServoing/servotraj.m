function theta_z = servotraj( t )

f   = 2*t;
x = mod(f,2);

if x <= 1
   theta_z = 0;%*x;
elseif x <= 2
   theta_z = pi/2;%*(x-1);
end
   
end


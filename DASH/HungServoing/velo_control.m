function u = velo_control(vel, vel_d)
%     t = t/4;
% 	tgt_vel = 3 * (servotraj(t) - pos(4:5));
	u = .05 * (vel_d - vel);
end

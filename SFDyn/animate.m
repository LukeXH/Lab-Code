figure(3)
close 3
figure(3)

T = sol.x(1):.01:sol.x(end);

handle = plot([0, 1, 2], [0, 0, 0]);
axis equal

for cur_t = T
	cur_state = deval(sol, cur_t);
	th1 = cur_state(1);
	th2 = cur_state(2);
	set(handle, 'XData', [0, sin(th1), sin(th1)+sin(th2)])
	set(handle, 'YData', [0, -cos(th1), -cos(th1)+cos(th2)])
	axis([-2 2 -2 2])
	drawnow
end

clear
close all

% addr = 'fe80::9f:b9dd:6fa4:1e7d';
addr = '::1';
t = tcpip(addr, 27015, 'NetworkRole', 'client');
t.ReadAsyncMode = 'continuous';
fopen(t);
tic;
disp('Sending cmds')
while( toc < 100)
if ( 2.5 < mod(toc,4))
    fwrite(t,uint8([0, 255]));  % [right_fwd, left_fwd]
else
    fwrite(t,uint8([255, 0]));  % [right_fwd, left_fwd]
end

pause(0.05);
end
disp('Done sending cmds');
fclose(t);
disp('Exited')
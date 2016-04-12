clear
close all

% addr = 'fe80::9f:b9dd:6fa4:1e7d';
addr = '::1';
t = tcpip(addr, 27015, 'NetworkRole', 'client');
t.ReadAsyncMode = 'continuous';
fopen(t);
fwrite(t,'Hello');
tic;
while(t.BytesAvailable == 0 && toc < 20)
   pause(0.001); 
end
if (toc < 20)
    data = fread(t,t.BytesAvailable);
    disp(data);
else
    disp('No response.')
end
    
fclose(t);
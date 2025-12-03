server = tcpserver("163.1.249.17",30000)

data = read(server,server.NumBytesAvailable,"double");
plot(data);

data = double([1:64]);
write(server, data)
write(server,'hello world','char')
write(server,"hello world",'string')

string = 'hello world';
string = [string 0 0 0 0];
string = string(1:end-mod(length(string),4));
write(server,string,'string');

read(server,server.NumBytesAvailable, "double");

read(server,server.NumBytesAvailable, "string")

configureTerminator(server,"LF")

%% UDP check

u = udpport("IPV4");

write(u,'/i0001',"uint8","163.1.249.37",30000);

string = 'hello world';
string = [string 0 0 0 0];
string = string(1:end-mod(length(string),4));

write(u,string,"string","163.1.249.37",30000);

string = ',';
string = [string 0 0 0 0];
string = string(1:end-mod(length(string),4));


write(u,'/start  ,f  B4  ',"163.1.249.37",30000);




bonsai_oscsend(u,'/start',"163.1.249.37",30000,'i',45);

data2 = ['/start' 0 0 ',i' 0 0 0 0 0 '-'];
write(u,data2,"163.1.249.37",30000);



%% check if you can do two servers in one matlab
server2 = tcpserver("0.0.0.0",20000)
write(server2,"hello world","string")

write(server2, data, "double")

read(server2, 64, "double")

read(server2,12,"double")

writeline(server2, 'hello')
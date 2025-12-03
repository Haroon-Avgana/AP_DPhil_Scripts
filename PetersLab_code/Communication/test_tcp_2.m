client = tcpclient("localhost",30000, 'Timeout', 20)
data = sin(1:64);
plot(data);

write(client,data,"double")

read(client, client.NumBytesAvailable, "double")

read(client, client.NumBytesAvailable,"string")

% configureTerminator(client,"CR",10)
write(client, 'hello')

% test date
date_today = string(datetime('today', 'Format', 'dd_MM_yyyy'));
write(client, date_today)

configureCallback(client,"terminator",@callbackFcn)


%% for check multiple connections in one matlab
client2 = tcpclient("localhost",20000)
read(client2, client2.NumBytesAvailable,"string")
read(client2, 12, "double")
write(client2,data,"double")
read(client2, 64, "double")

readline(client2)
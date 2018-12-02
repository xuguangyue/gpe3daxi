clear
clc
t=0.1:0.1:10;
zt=load('zt.dat');
zt=zt(1:100);
wx=load('wxt.dat');
wx=wx(1:100)/(2*pi);
wz=load('wzt.dat');
wz=wz(1:100)/(2*pi);

figure
plot(t,zt)
figure
plot(t,wx,t,wz)

clear all
close all
xy=load('xymast.dat');
[m,n]=size(xy);
x=xy(1:m/2,1:n);
y=xy((m/2+1):m,1:n);
h=load('depth.dat');

mesh(x,y,-h)
figure(2)
subplot(211)
plot(x(1,:),-h(1,:))
grid
subplot(212)
plot(x(1,:),-h(1,:))
axis([-3000 3000 -15 5])
grid

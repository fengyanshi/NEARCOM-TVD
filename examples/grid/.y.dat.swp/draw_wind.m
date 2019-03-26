clear all
close all
xy=load('xymast.dat');
[m,n]=size(xy);
m1=m/2;
n1=n;
x=xy(1:m/2,1:n);
y=xy((m/2+1):m,1:n);
wxy=load('../Wind/wind01.dat');
wx=wxy(1:m/2,1:n);
wy=wxy(m/2+1:m,1:n);

sk=5;
quiver(x(1:sk:m1,1:sk:n1),y(1:sk:m1,1:sk:n1),wx(1:sk:m1,1:sk:n1),wy(1:sk:m1,1:sk:n1))

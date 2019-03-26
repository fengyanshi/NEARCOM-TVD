
clear all

fdir='/Users/fengyanshi15/tmp5/';
fdir_input='../input/';
X=load([fdir_input 'x_str.txt']);
Y=load([fdir_input 'y_str.txt']);
[dx,tmp]=gradient(X,1,1);
[tmp,dy]=gradient(Y,1,1);
Jaco=dx.*dy;
dep=load([fdir_input 'dep_str.txt']);

fac=150;
xyaxis=[-6000 6000 -5000 3000];
skx=2;
sky=2;

nstart=input('nstart=');
nend=input('nend=');

set(gcf,'units','inches','paperunits','inches','papersize', [7 8],'position',[0 0 7 8],'paperposition',[0 0 7 8]);

icount=0;
for num=nstart:1:nend
icount=icount+1;

fname=sprintf('%.4d',num);
ele=load([fdir 'eta_' fname]);
mask=load([fdir 'mask_' fname]);
u=load([fdir 'u_' fname]);
v=load([fdir 'v_' fname]);
hs=load([fdir 'hs_' fname]);

uu=sqrt(u.^2+v.^2);
[ux,uy]=gradient(u,1,1);
[vx,vy]=gradient(v,1,1);
ww=(uy-vx);

[n,m]=size(ele);

cc=ele;
cc(mask<1)=NaN;
B=-dep;
B(mask==1)=hs(mask==1);
B(mask<1)=NaN;
u(mask<1)=NaN;
v(mask<1)=NaN;
ww(mask<1)=NaN;

clf
subplot(211)
pcolor(X,Y,B),shading interp
%caxis([0 2.0])
colorbar
%demcmap(B);
hold on
%pcolor(X,Y,cc),shading interp
%colormap(jet)
quiver(X(1:sky:end,1:skx:end),Y(1:sky:end,1:skx:end),u(1:sky:end,1:skx:end)*fac,v(1:sky:end,1:skx:end)*fac,0,'Color','y');

xlabel('East (meters) ','fontsize',12,'fontweight','bold');
ylabel('North (meters) ','fontsize',12,'fontweight','bold');
axis(xyaxis);
title(['color: wave height' '   (t=' num2str((num-1)*5) ' min)'])

subplot(212)
pcolor(X,Y,cc),shading interp
%caxis([0 2.0])
colorbar
hold on
quiver(X(1:sky:end,1:skx:end),Y(1:sky:end,1:skx:end),u(1:sky:end,1:skx:end)*fac,v(1:sky:end,1:skx:end)*fac,0,'Color','y');

xlabel('East (meters) ','fontsize',12,'fontweight','bold');
ylabel('North (meters) ','fontsize',12,'fontweight','bold');
axis(xyaxis);
title(['color: surface elevation' '  (t=' num2str((num-1)*5) ' min)'])

%title(['Time = '  num2str((num(1))*5) ' min'])
%M(:,icount)=getframe(gcf);
pause(0.1)
end


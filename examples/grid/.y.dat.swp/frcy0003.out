      parameter(m=251,n=401,m1=126,n2=201)
      real W(m,n),Pa(m,n),x(m,n),y(m,n),Wx(m,n),Wy(m,n)


       open(1,file='xymast.dat')
       do j=1,n
        read(1,*)(x(i,j),i=1,m)        
       enddo
       do j=1,n
        read(1,*)(y(i,j),i=1,m)
       enddo
       close(1)

c ---- parameters
       pi=3.1415926
       theta=20.0*pi/180.

       Rm=50.*1000.

c -- moving speed  and angle
       U0=10.0
       angle=120.0*pi/180.
c ---  maximum Wind       
       Wm=50.0
       
       Po=0.8
       Pm=1.0

c --- landing site
      xland=3.*1000.
      yland=0.*1000.

      dt=3600.
      ktotal_hour=24.
      x_start=-ktotal_hour*U0*cos(angle)
      y_start=-ktotal_hour*U0*sin(angle)

      do k=1,ktotal_hour
       xc=x_start+U0*cos(angle)*(k-1.)*dt
       yc=y_start+U0*sin(angle)*(k-1.)*dt
      
       do j=1,n
       do i=1,m
        r=sqrt((x(i,j)-xc)**2+(y(i,j)-yc)**2)
        A=-((x(i,j)-xc)*sin(theta)+(y(i,j)-yc)*cos(theta))
        B=((x(i,j)-xc)*cos(theta)-(y(i,j)-yc)*sin(theta))
        if (r.le.Rm)then
         Wx(i,j)=r/(Rm+r)*U0*cos(angle)+Wm*sqrt(r)/(Rm)**(1.5)*A
         Wy(i,j)=r/(Rm+r)*U0*sin(angle)+Wm*sqrt(r)/(Rm)**(1.5)*B
         Pa(i,j)=P0+1./4.*(Pm-P0)*(r/Rm)**(3.)
        else
         Wx(i,j)=r/(Rm+r)*U0*cos(angle)+Wm*sqrt(Rm)/(r)**(1.5)*A
         Wy(i,j)=r/(Rm+r)*U0*sin(angle)+Wm*sqrt(Rm)/(r)**(1.5)*B
         Pa(i,j)=Pm-3./4.*(Pm-P0)*(Rm/r)
        endif
       enddo
       enddo
     
      open(2,file='../Wind/wind01.dat')
       do j=1,n
        write(2,100)(Wx(i,j),i=1,m)
       enddo
       do j=1,n
        write(2,100)(Wy(i,j),i=1,m)
       enddo
       close(2)

       enddo
      
100    format(800f12.2)     

       end

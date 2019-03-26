      subroutine prtplt (n)
c
c ======================================================================
c
c   Purpose -
c     plot and provide formatted writes to paper and film,
c     where    = 1: initial mesh and problem data
c            n = 2: time step and cycle info
c              = 3: field variables
c
c   PRTPLT is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          RIPPLE    NEWCYC
c
c
c   PRTPLT calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c          GLOBAL    ripple
c
c ======================================================================
c
c##############################################################
       implicit real*8 (a-h,o-z)
c       include "32bit.h"
c##############################################################
c
c############
       include "comdk1.h" 
c############
	 character*3 name(1)
	real dummy1(imax,jmax),dummy2(imax,jmax),dummy3(imax,jmax)
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c     ------------
      call global
c     ------------
c
	 open(1,file='data_waveu.dat')
	 open(2,file='data_wavev.dat')
	 open(3,file='data_wavep.dat')
	 open(4,file='data_wavee.dat')

      go to (20,120,140), n
c
c.... prtplt (1) write out mesh data
c
   20 continue
 
c     ----------
c
c.... write out mesh information
c
	open(20,file='data_status.dat') 
       open(21,file='data_xi.dat')
       open(22,file='data_yj.dat')
	open(23,file='data_ar.dat')

	write(20,*) imax,jmax
       write (21,99) (x(i),i=1,imax) 
       write (22,99) (y(j),j=1,jmax) 
 
	do  j=1,jmax
	write(23,99) (ar(ij),ij=(j-1)*imax+1,(j-1)*imax+imax)	
	end do	

	write(4,21) t,(x(i),i=1,imax) 

c	 output in   y level

c	write(1,21) t,(x(i),i=1,imax) 
c	write(2,21) t,(x(i),i=1,imax) 
c	write(3,21) t,(x(i),i=1,imax) 

c	 output in   x level

	write(1,21) t,(yj(j),j=1,jmax) 
	write(2,21) t,( y(j),j=1,jmax)   
	write(3,21) t,(yj(j),j=1,jmax)  


 99	format(3x,800(f10.4)) 

	close(20)
	close(21)
	close(22)
	close(23)

      go to 9999
c
c.... prtplt (2)  write time step, cycle information
c
  120 continue
 
 	 open(99,file='cycle.txt')

	 write(*,210)iter,t,delt,ncyc
	 if(kl.eq.8) write(*,211) t,t-delt,deta

	 write(99,210)iter,t,delt,ncyc
	 if(kl.eq.8) write(99,211) t,t-delt,deta

	write(4,21) t,(etah(i),i=1,imax) 
 
c	 output in   y level
 
c	write(1,21) t,(u(ij),ij=imax+1, 2*imax) 
c	write(2,21) t,(v(ij),ij=imax+1, 2*imax)  
c	write(3,21) t,(p(ij),ij=imax+1, 2*imax)  



c	 output in   x level
 
	write(1,21) t,(u(ij),ij=100, (jmax-1)*imax+100, imax) 
	write(2,21) t,(v(ij),ij=100, (jmax-1)*imax+100, imax)  
	write(3,21) t,(p(ij),ij=100, (jmax-1)*imax+100, imax)  




21	 format(3x,800f10.2)

 210  format(3x,5hiter=,i6,3x,5htime=,e15.6,3x,5hdelt=,e12.6,3x
     1	,6hcycle=,i9)
 211   format(17x,5htime=,e15.6,3x,6htimen=,e12.6,3x,5hdeta=,f10.4)

      go to 9999
c
c.... prtplt (3)  write field variables to paper
c
  140 continue
c	write (13,300) prbname
c      write (13,330) iter,t,idt,itc,jtc,delt,ncyc,vchgt
c      write (13,350)
c
c      write (13,260) fvol,vvol,xmv,ymv,tke
c      write (13,310)
c      do 170 i=1,imax
c        do 170 j=1,jmax
c          ij=(j-1)*imax+i
c          imj=ij-1
c          ijm=ij-imax
c          dij=0.
c          if (j.eq.1.or.i.eq.1.or.j.eq.jmax.or.i.eq.imax) go to 160
c          dij=rdx(i)*(ar(ij)*r(i)*u(ij)-ar(imj)*r(i-1)*u(imj))+
c     1         rdy(j)*(at(ij)*ri(i)*v(ij)-at(ijm)*ri(i)*v(ijm))
c          if (ac(ij).gt.0.0) dij=dij/(ac(ij)*ri(i))
c  160     continue
c          write (13,320) i,j,u(ij),v(ij),p(ij),dij,f(ij),
c     1                   nf(ij)
c  170 continue
c
c      write (13,270) ncyc,t,delt,(i,i=2,imax,2)
c      do 180 jjj=1,jmax
c        j=jmax+1-jjj
c        ij=(j-1)*imax+1
c        write (13,280) j,(nf(i),i=ij,ij+imax)
c  180 continue
c
c      write (13,275) ncyc,t,delt,(i,i=2,imax,2)
c


c	output F function and velocity vectors  Qun


c	interval=80
c	kounter=ncyc/interval
c	if(mod(ncyc,interval).eq.0) then

c------------------------------------------------
c	 ASCII	  output	tested O.K. Qun

c	open(25,file='data_f.dat')
c	open(26,file='data_u.dat')
c	open(27,file='data_v.dat')

	  
c	do j=1,jmax
c	write(25,99) (f(ij),ij=(j-1)*imax+1,(j-1)*imax+imax) 	
c	write(26,99) (u(ij),ij=(j-1)*imax+1,(j-1)*imax+imax) 	
c	write(27,99) (v(ij),ij=(j-1)*imax+1,(j-1)*imax+imax) 	
c	end do


c------------------------------------------------
c	 BINARY	  output	 tested O.K. Qun

		
c       write(unit=name(1),fmt=244)100+kounter

c       open(unit=25,file='data_f.'//name(1),access='direct',
c	&  recl=4*(imax))
c       open(unit=26,file='data_u.'//name(1),access='direct',
c	&  recl=4*(imax))
c       open(unit=27,file='data_v.'//name(1),access='direct',
c	&  recl=4*(imax))


c	do  j=1,jmax
c	do  i=1,imax
c	ij=(j-1)*imax+i
c	dummy1(i,j)=f(ij)
c	dummy2(i,j)=u(ij)
c	dummy3(i,j)=v(ij)
c	end do
c	write(25,rec=j) (dummy1(i,j),i=1,imax)
c	write(26,rec=j) (dummy2(i,j),i=1,imax)
c	write(27,rec=j) (dummy3(i,j),i=1,imax)
c	end do
c	close(25) 
c	close(26) 
c	close(27) 



c	end if


244	format(I3) 

  190 continue
c
 9999 return
c
  260 format (10x," fluid volume ",1pe14.6," void volume ",1pe14.6,/,
     &        10x," x momentum ",1pe14.6," y momentum ",1pe14.6,/,
     &        10x," total ke ",1pe14.6)
  270 format (1h1,1x,43hnf field (incl. fictitious cells) for cycle,i6,
     1         2x,2ht=,1pe14.7,2x,5hdelt=,e12.5//5x,32i4)
  275 format (1h1,1x,43hnw field (incl. fictitious cells) for cycle,i6,
     1         2x,2ht=,1pe14.7,2x,5hdelt=,e12.5//5x,32i4)
  280 format (1x,i3,1x,63i2/(5x,63i2))
  290 format (1h1)
  300 format (a80)
  310 format (4x,1hi,5x,1hj,9x,1hu,14x,1hv,15x,1hp,15x,1hd,12x,
     1         1hf,11x,2hnf)
  320 format (2x,i3,3x,i3,6(3x,1pe12.5),3x,i3)
  330 format (1x," iter = ",i5," time = ",1pe12.5," dt",a2,"(",
     &        i2,",",i2,") = ",1pe12.5," cycle = ",i6," vchgt = ",
     &        1pe12.5)
  350 format (1h0)
  360 format (1h ,18x,a80,1x,a10,2(1x,a8))
  370 format (2x,5hnkx= ,i4)
  380 format (2x,8hmesh-x= ,i4,3x,4hxl= ,1pe12.5,3x,4hxc= ,e12.5,3x,
     1       4hxr= ,e12.5,3x,5hnxl= ,i4,3x,5hnxr= ,i4,3x,6hdxmn= ,e12.5)
  390 format (2x,8hmesh-y= ,i4,3x,4hyl= ,1pe12.5,3x,4hyc= ,e12.5,3x,
     1       4hyr= ,e12.5,3x,5hnyl= ,i4,3x,5hnyr= ,i4,3x,6hdymn= ,e12.5)
  400 format (2x,5hnky= ,i4)
  430 format (13x,i3,2x,1pe12.5,3x,e12.5)
      end

      subroutine prtplt2 (n)

c

c ======================================================================

c

c   Purpose -

c     plot and provide formatted writes to paper and film,

c     where    = 1: initial mesh and problem data

c            n = 2: time step and cycle info

c              = 3: field variables

c

c   PRTPLT is called by -

c

c          name      name      name      name      name      name

c        --------  --------  --------  --------  --------  --------

c          RIPPLE    NEWCYC

c

c

c   PRTPLT calls the following subroutines and functions -

c

c          name      source    name      source    name      source

c        --------  --------  --------  --------  --------  --------

c          GLOBAL    ripple

c

c ======================================================================

c

c##############################################################

       implicit real*8 (a-h,o-z)

c       include "32bit.h"

c##############################################################

c

c############

       include "comdk1.h" 
	character*4 fname
	character*10 fdir
	character*1 ntype
	data kounter/0/

	fdir='../result/' 

      call global

      go to (20,120,140), n


   20 continue


c.... write out mesh information

c

	open(2,file=fdir//'data_status.dat') 
        open(3,file=fdir//'data_xi.dat')
        open(4,file=fdir//'data_yj.dat')
	open(23,file=fdir//'data_ar.dat')
	open(24,file=fdir//'data_at.dat')


	write(2,*) imax,jmax,prtdt,nportype
      write (3,99) (x(i),i=1,imax) 
      write (4,99) (y(j),j=1,jmax) 

	do  j=1,jmax
	write(23,99) (ar(ij),ij=(j-1)*imax+1,(j-1)*imax+imax)	
	end do	

	do  j=1,jmax
	write(24,99) (at(ij),ij=(j-1)*imax+1,(j-1)*imax+imax)	
	end do	

	do 901 ii=1,nportype

	write(ntype(1:1),'(I1)')ii

	open(19,file=fdir//'data_pr'//ntype//'.dat')
	do  j=1,jmax
	write(19,99) (pr(ii,ij),ij=(j-1)*imax+1,(j-1)*imax+imax)	
	end do	
	close(19)
	open(19,file=fdir//'data_pbdy'//ntype//'.dat')
	do  j=1,jmax
	write(19,99) (pbdy(ii,ij),ij=(j-1)*imax+1,(j-1)*imax+imax)	
	end do	
	close(19)
901	continue

 99	format(3x,500(f16.4)) 
	close(2)
	close(3)
	close(4)
	close(23)
	close(24)

      go to 9999

c.... prtplt (2)  write time step, cycle information

c

  120 continue

c	write (13,330) iter,t,idt,itc,jtc,delt,ncyc,vchgt
c     write (13,260) fvol,vvol,xmv,ymv,tke

	write(*,331)iter,t,delt,ncyc

 331  format(3x,5hiter=,i6,3x,5htime=,e15.6,3x,5hdelt=,e12.6,3x

     1	,6hcycle=,i9)

      go to 9999

c
c.... prtplt (3)  write field variables to paper

  140 continue

c	output F function and velocity vectors  Qun

	kounter=kounter+1

	n_first=mod(kounter/1000,10)
	n_second=mod(kounter/100,10)
	n_third=mod(kounter/10,10)
	n_fourth=mod(kounter,10)


	write(fname(1:1),'(I1)')n_first
	write(fname(2:2),'(I1)')n_second
	write(fname(3:3),'(I1)')n_third
	write(fname(4:4),'(I1)')n_fourth

	open(25,file=fdir//'data_f.'//fname)
	open(26,file=fdir//'data_u.'//fname)
	open(27,file=fdir//'data_v.'//fname)
	open(28,file=fdir//'data_p.'//fname)
	open(29,file=fdir//'data_h.'//fname)
	open(30,file=fdir//'data_n.'//fname)
	open(31,file=fdir//'data_k.'//fname)

	print*,'printing...',fname

	do j=1,jmax
	write(25,99) (f(ij),ij=(j-1)*imax+1,(j-1)*imax+imax) 	
	write(26,99) (u(ij),ij=(j-1)*imax+1,(j-1)*imax+imax) 	
	write(27,99) (v(ij),ij=(j-1)*imax+1,(j-1)*imax+imax) 
	write(28,99) (p(ij),ij=(j-1)*imax+1,(j-1)*imax+imax) 	
	write(30,99) (xnutt(1,ij),ij=(j-1)*imax+1,(j-1)*imax+imax)
	write(31,99) (ken(1,ij),ij=(j-1)*imax+1,(j-1)*imax+imax) 	
	end do

	write(29,99) (etah(i),i=1,imax) 

	close(25)
	close(26)
	close(27)
	close(28)
	close(29)
	close(30)
	close(31)

c	open(19,file='tmp.dat')
c	do j=1,jmax
c	write(19,99) (rkc(1,ij),ij=(j-1)*imax+1,(j-1)*imax+imax) 	
c	end do
c	close(19)

244	format(I3) 

  190 continue

 9999 return

  260 format (10x," fluid volume ",1pe14.6," void volume ",1pe14.6,/,

     &        10x," x momentum ",1pe14.6," y momentum ",1pe14.6,/,

     &        10x," total ke ",1pe14.6)

  270 format (1h1,1x,43hnf field (incl. fictitious cells) for cycle,i6,

     1         2x,2ht=,1pe14.7,2x,5hdelt=,e12.5//5x,32i4)

  275 format (1h1,1x,43hnw field (incl. fictitious cells) for cycle,i6,

     1         2x,2ht=,1pe14.7,2x,5hdelt=,e12.5//5x,32i4)

  280 format (1x,i3,1x,63i2/(5x,63i2))

  290 format (1h1)

  300 format (a80)

  310 format (4x,1hi,5x,1hj,9x,1hu,14x,1hv,15x,1hp,15x,1hd,12x,

     1         1hf,11x,2hnf)

  320 format (2x,i3,3x,i3,6(3x,1pe12.5),3x,i3)

  330 format (1x," iter = ",i5," time = ",1pe12.5," dt",a2,"(",

     &        i2,",",i2,") = ",1pe12.5," cycle = ",i6," vchgt = ",

     &        1pe12.5)

  350 format (1h0)

  360 format (1h ,18x,a80,1x,a10,2(1x,a8))

  370 format (2x,5hnkx= ,i4)

  380 format (2x,8hmesh-x= ,i4,3x,4hxl= ,1pe12.5,3x,4hxc= ,e12.5,3x,

     1       4hxr= ,e12.5,3x,5hnxl= ,i4,3x,5hnxr= ,i4,3x,6hdxmn= ,e12.5)

  390 format (2x,8hmesh-y= ,i4,3x,4hyl= ,1pe12.5,3x,4hyc= ,e12.5,3x,

     1       4hyr= ,e12.5,3x,5hnyl= ,i4,3x,5hnyr= ,i4,3x,6hdymn= ,e12.5)

  400 format (2x,5hnky= ,i4)

  430 format (13x,i3,2x,1pe12.5,3x,e12.5)

      end

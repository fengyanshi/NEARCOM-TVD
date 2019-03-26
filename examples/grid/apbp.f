	subroutine apbp
c ======================================================================
c
c   Purpose - calculate frictional force induced by viscous effect ap
c             and turbulence effect bp (Liu et al, 1999), fyshi 10/25/02
c             	   
c ======================================================================
       implicit real*8 (a-h,o-z)
c       include "32bit.h"
c##############################################################
c
c############
       include "comdk1.h" 
c############

        do 410 j=2,jm1
        do 410 i=2,im1
        ij=(j-1)*imax+i	
	do 333 ii=1,nportype
	if(pr(ii,ij).gt.0.0.or.pt(ii,ij).gt.0.0)then	  
	  ap(ii,ij)=alphap(ii)*(1.-porosity(ii))**2
     &         /(porosity(ii)*porosity(ii)*porosity(ii))
     &          *xnumo/grav/Dp50(ii)/Dp50(ii)
	  ucc=sqrt(un(ij)*un(ij)+vn(ij)*vn(ij))
	  if(ucc.lt.1.)ucc=1.0
	  rkc(ii,ij)=ucc*Ttypical/(porosity(ii)*Dp50(ii))
          bp(ii,ij)=betap(ii)*(1.+7.5/rkc(ii,ij))*(1.-porosity(ii))
     &      /(porosity(ii)*porosity(ii)*porosity(ii))/grav/Dp50(ii)
	endif
333	continue	
410	continue

	return
	end
      parameter (ibar2=1000,jbar2=1000,meshx=8,meshy=8,
     &           nobd=20,nfrsrfd=20,nprd=200,ntyp=10,nshp=10)
      parameter (mshx=meshx+1, mshy=meshy+1, nxy=ibar2*jbar2)
      parameter (nequib=5*nxy,itrcmax=10000)
      parameter (mn=30,max=3,c_d=0.2,c_s=0.2,sigmak=0.2)
      character*80 prbname
      character*8 dat,jnm,ochn,tim
      character*2 idt
      logical plots,dump,upright,gfnctn,sym,
     &        smooth,conserve,emin,cray,unfrmmsh
C....added for PLIC scheme
C     &	      csf_normal,csf_tangential,slic, plic, unsplit
C....added for PLIC scheme
C      integer vof_particles

      real*8 nxo,mxo,nyo,myo,nxf,mxf,nyf,myf,kappa,jcb,mod1,
     &		nxp,mxp,nyp,myp
C	real*8 nxx,nyy
c
      common /begini/ lcbegin(1),kk,aa,h0,h0r,xcenter,alpha1,c1,
     1		xstart,tstart,tfinish,ybc(2000,40),ubc(2000,40),
     1          tstart_a,tfinish_a,prtdt_a,prtdt_t,twprt_a,
     2          vbc(2000,40),ddux,ddvy,uxmx,vymn,xxxf,
     3          xxxu,xxxv,xxxk,xxxt,xxxp,ttend,xxl,xxt,xxk,
     3		xxlr,xxkr,c1r,crest,
     3		segma,fxstart,t_pad(2000),u_pad(2000),h_pad(2000),
     4		x_pad(2000),ulinear,usecond,areturn

	common /addition2/ det_k(100),ytrough,mod1
	common /addition3/ j_det,istrflag,
     1	nrs,nloc,nanimation,nout(100)
        common /addition4/ xout(100),fux_t(100,250),fuy_t(100,250),
     1	visxx_t(100,250),visxy_t(100,250),fvx_t(100,250),
     2	fvy_t(100,250),visyx_t(100,250),visyy_t(100,250)
        common /addition14/ xpol(10),ypol(10),conc(10)
	common /addition0/ ttmax,yjet,ujet,vjet,wjet,ylow,yup,xxc
     4		xdis,eta0old,eta1old,eta2old,eta2aold,eta0x,eta2x,
     5		etaimax,c01,c02,cim1,crbc,
     5		eta0,eta1,eta1a,eta2,uleft1,vleft1,xmass0,xmass_d,
     5		xsaka,xmass_a,xmass_a1,xmass_a2,cnf,petit,
     5		umean,ratio,uturb,
     5		ssource,tsource,aawave(100),twave(100),cwave(100),
     5		swave(100),nwave,isources,isourcee,jsources,jsourcee,
     5		nsource,islip,ncenter,ibg,ieg,jbg,jeg,
     6		it,nnn,istable,interx,intery,nexit,
     7		istablemx,kemodel,npollutant,nweakref,nopen
c
      common /fldvarr/ un(nxy),vn(nxy),pn(nxy),fn(nxy),tempv(nxy),
     1       um(nxy),vm(nxy),u1u1(nxy),u2u2(nxy),u3u3(nxy),u1u2(nxy),
     1		       u1u1n(nxy),u2u2n(nxy),u3u3n(nxy),u1u2n(nxy),
     2                 tempu(nxy),fsv(nxy),div(nxy),ar(nxy),at(nxy),
     3	rhsx(nxy),rhsy(nxy),rhsxold(nxy),rhsyold(nxy),ftmp(nxy)

	common /addition1/ ac(nxy),curlu(nxy),kappa(nxy),
     4                 streamf(nxy),cvol(nxy),vpot(nxy),tss(nxy),
     5                 tauxx(nxy),tauyy(nxy),tauxy(nxy),ttots(nxy),
     6                 ftilde(nxy),u(nxy),v(nxy),p(nxy),f(nxy),xpn(nxy),
     7                 xnut(nxy),xnutn(nxy),xk(nxy),xkn(nxy),xp(nxy),
     8		       xnuty(nxy),xnutyn(nxy),xep(nxy),xepn(nxy),
     2	xmint,amp1,amp2,amp3,amp4,amp5,realize,eddycoef,
     9	c1e,c2e,c_mu,sege,segk,xmaxxnut,ymaxxnut,c1ea,c1eb,ticf,y_coef,
     8  shtime
c
        common /addition5/ gradv_ct(nxy),gradu_ct(nxy),
     1	gradv_rlu(nxy),gradu_udr(nxy),Utmp(nxy),Vtmp(nxy),
     2	tauxyo(nxy),tauxxo(nxy),tauyyo(nxy),product(nxy),
     3  blevel(nxy),blevelo(nxy),shield(nxy),shieldo(nxy)

      common /porousmk/ pa2(ntyp,nprd),pa1(ntyp,nprd),pb2(ntyp,nprd),
     1     pb1(ntyp,nprd),
     1     pc2(ntyp,nprd),pc1(ntyp,nprd),xpor(ntyp,nxy),ypor(ntyp,nxy),
     2     pd1(ntyp,nprd),nxp(ntyp,nprd),pd2(ntyp,nprd),ipr(ntyp,nprd),
     3     mxp(ntyp,nprd),pe1(ntyp,nprd),nyp(ntyp,nprd),pe2(ntyp,nprd),
     4	   myp(ntyp,nprd),
     8     nportype,pc(ntyp,nxy),pt(ntyp,nxy),pr(ntyp,nxy),
     1     iph(ntyp,nprd),ijpor(ntyp,nxy),
     2     pbdy(ntyp,nxy),
     3     ap(ntyp,nxy),bp(ntyp,nxy),rkc(ntyp,nxy),
     5     xperm(100),xporosity(100),d50(100),
     6     xalpha(100),xbeta(100),xgamma(100),
     7	   xa(100),xxb(100),gc(100),npc(nxy)


      common /fldvari/ nf(nxy),nfold(nxy),nfs(nxy)
c
      common /geomr/ jcb(nxy),alp(nxy),gam(nxy)
c
      common /meshr/x(ibar2),xi(ibar2),delx(ibar2),etah(ibar2),
     1             etatmp(ibar2),rdx(ibar2),rx(ibar2),y(jbar2),
     2		   yj(jbar2),dely(jbar2),rdy(jbar2),xl(mshx),
     3             xc(meshx),dxmn(meshx),yl(mshy),yc(meshy),
     4             dymn(meshy),r(ibar2),ri(ibar2),
     5             xb(5),yb(5),etahn(ibar2)
c
      common /meshi/ nxl(meshx),nxr(meshx),nyl(meshy),nyr(meshy)
c
      common /voidr/ fvol,vvol
c
      common /intvarr/ delt,t,autot,prtdt,twprt,pltdt,twplt,
     1             twfin,xnu,dtend,dmpdt, rhof, vchgt, ttl,
     2             sigma, cyl, gx, gy, ui, vi, utop,
     3             alpha, beta, flgc, xmin, xmax, ymin,
     4             dangle, bond, omcyl, ymax,tmv,tke,
     5             dtvis, dudr,dudl,dudt,dudb,dvdr,dvdl,dvdt,dvdb,
     6             tquit,tbeg,dtsft,dxmin,dymin,psat,tseti,tsetf,
     7             dtmax,uinf(4),vinf(4),con,fcvlim,
     8             erriccg,xmu,cangler,xmv,ymv,opalp,omalp,
     9             canglet,cangleb,canglel,frctn,twdmp,vanleer
C....added for PLIC scheme
C     1             dxmax,dymax
c
      common /intvari/ ibar,jbar,imax,jmax,im1,jm1,nkx,nky,ncyc,
     1             isurf10,icyl,kl,kr,kb,kt,iter,
     2             nocon,nflgc,
     3             i,j,liter,npack,idiv,ibcextrap,
     4             ibcflg,ibcfinal,itc,jtc,conserve,cray,iotty,
     5             ivis,jvis,isft,jsft,itmxiccg,ndump,
     6             nrestart,sym
C....added for PLIC scheme
C     7		   ifin,ired
c
      common /constr/ emf,emf1,em6,em10,ep10,pi,tpi,rpd,
     1                em6p1,em61,fact
c
      common /obstr/ oa2(nshp,nobd),oa1(nshp,nobd),
     1       ob2(nshp,nobd),ob1(nshp,nobd),
     1       oc2(nshp,nobd),oc1(nshp,nobd),xobs(nshp,nxy),
     2       yobs(nshp,nxy),od1(nshp,nobd),nxo(nshp,nobd),
     2       od2(nshp,nobd),mxo(nshp,nobd),oe1(nshp,nobd),
     3       nyo(nshp,nobd),oe2(nshp,nobd),myo(nshp,nobd)
C....added for PLIC scheme
C     4      xanf(nobd), xend(nobd), yanf(nobd), yend(nobd)
c

      common /obstrs/ ar0(nxy),at0(nxy),ac0(nxy)

      common /obsti/ nmean,istart,iend,iinterval,
     1		nobs(nshp),ioh(nshp,nobd),ijobs(nshp,nxy),nobshp
c
      common /mean/ tmstart,tmend,tinterval,
     1	eflux(ibar2),efluxt(ibar2),disp(ibar2),ep(ibar2),ep0(ibar2),
     1	ek(ibar2),et(ibar2),etaa(ibar2)
c
      common /frsrfr/ fa2(nfrsrfd),fa1(nfrsrfd),
     1                fb2(nfrsrfd),fb1(nfrsrfd),fc2(nfrsrfd),
     2                fc1(nfrsrfd),ifh(nfrsrfd),fd1(nfrsrfd),
     3                nxf(nfrsrfd),fd2(nfrsrfd),mxf(nfrsrfd),
     4                fe1(nfrsrfd),nyf(nfrsrfd),fe2(nfrsrfd),
     5                myf(nfrsrfd),flht,frsurf,deltemin,
     6                yfrsrf(nxy),xfrsrf(nxy)
c
      common /frsrfi/ nfrsrf,nsmooth,smooth,emin,
     1                iequib,itmxemin,itskpemn,
     2                upright,ijfr(nxy)
C    3		      nxx(nxy),nyy(nxy)
C....added for PLIC scheme
C     3                csf_normal, csf_tangential,
C     4                slic, plic, unsplit
c
      common /tensr/ tensx(nxy),tensy(nxy),gradrox(nxy),
     1               gradroy(nxy),gradnwx(nxy),gradnwy(nxy)
c
      common /tensi/  gfnctn
c
      common /ghostr/ pbc(4),pbcl,pbcr,pbct,pbcb
c
      common /labeli/ prbname,dat,jnm,ochn,tim,idt
c
      common /graphicr/ vmxfrctn,scale,xf(5),yf(5)
c
      common /graphici/ plots,dump,iout(40),
     1                  ixsymplt,iysymplt,nfrsrfpt,nobspt,
     2                  nvar,ixskip,iyskip,icolor,unfrmmsh
c
      common /tracer/ ntrc_xy(2,itrcmax),trc_xy(2,itrcmax)
      common /tracer1/ trc_vel(2,itrcmax),ntracer,nfree,nextra,npor

      dimension zplt(1),zmat(1,1),pltbgn(1)
 
c
      equivalence (zplt(1),pn(1)),(zmat(1,1),pn(1)),
     1            (pltbgn(1),un(1))
c

c
c	wave drivers
c
	   real*8 xkn1(mn)

       real*8 waveh,waveht,wavet,wavek,wavec,wavel, 
     1	wavew,wavewt,eta10,rc,xle,umk,deta 


       common /wavecase/ waveh,waveht,wavet,wavek,wavec,wavel,
     1	 xkn1,wavew,wavewt,rc,xle,umk,deta,iwave 

        real*8  kappa1,K1,E2,lmr,mu 

	  common/wavecn/wahd,ursell,kappa1,K1,E2,lmr,mu,a0(0:3),
     1      b0(0:3,0:2)

       common /porous/porosity(ntyp),gammaa(ntyp),ca(ntyp),
     1         grav,alphap(ntyp),betap(ntyp),Ttypical,
     1         Dp50(ntyp),nporous(ntyp),xnumo  

       common /smag/xmut(nxy),xmusmag(nxy),turl,
     1	      turr,turcoef,turmax
	   
	   	common/turb/dkl,dkr,dkt,dkb,
     1        sij(nxy),visx_vt(nxy),visy_vt(nxy),
     2 	      prod(max,nxy),ken(max,nxy),ke(max,nxy),
     3         eps(max,nxy),xnutt(max,nxy),length(max+1,nxy),
     4        itur 

	    real*8 ken,ke,length	
       
	          





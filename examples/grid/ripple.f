 
      program ripple
c
c ======================================================================
c
c   Purpose -
c     RIPPLE: - An Eulerian hydrocode in 2 dimensions for the
c               solution of incompressible flows with free
c               surfaces using the volume of fluid (VOF) method
c               * * * Fluid Dynamics Group T-3, LANL * * *
c
c     Copyright, 1990, the regents of the University of California.
c     This software was produced under a U.S. Government contract
c     (W-7405-ENG-36) by the Los Alamos National Laboratory, which
c     is operated by the University of California for the U.S.
c     Department of Energy.  The U.S. Government is licensed to use,
c     reproduce, and to distribute this software.  Permission is
c     granted to the public to copy and use this software without
c     charge, provided that this notice and any statement of
c     authorship are reproduced on all copies.  Neither the government
c     nor the university makes any warranty, express or implied,
c     or assumes any liability or responsibility for the use of this
c     software.
c
c   RIPPLE is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c            none
c
c
c   RIPPLE calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c          SECOND    system     SETUP    ripple        BC    ripple
c         SRFFRCE    ripple     BEGIN    ripple    NEWCYC    ripple
c         IMPLCTP    ripple   TENSION    ripple    VOFADV    ripple
c          PRTPLT    ripple    VTILDE    ripple
c
c ======================================================================
c
c * * * * * * * * * * * * * * * * * * PROBLEM SETUP * * * * * * * * * * * * * * * * * *
c       
c
c      * * * * * * * * MEMORY PARAMETERS: PARAMETER statement variables * * * * * * * *
c       
c      IBAR2 is the upper dimension bound for element ``i'' in common block arrays of the 
c        form ARRAY(i,j) or X(i), where ARRAY is some field variable array and X is some 
c        x-coordinate mesh array.  IBAR2 must be greater than or equal to the total number of 
c        cells in the x (r) direction.
c      JBAR2 is the upper dimension bound for element ``j'' in common block arrays of the 
c        form ARRAY(i,j) or Y(j), where ARRAY is some field variable array and Y is some
c        y-coordinate mesh array.  JBAR2 must be greater than or equal to the total 
c        number cells in the y (z) direction.
c      MESHX is upper dimension bound for common block arrays used in generating
c        mesh quantities in the x (r) direction. 
c      MESHY is upper dimension bound for common block arrays used in generating
c        mesh quantities in the y (z) direction.
c      NOBD is the upper dimension bound for the obstacle conic function coefficient arrays.
c      NFRSRFD is the upper dimension bound for the free surface conic function coefficient 
c        arrays.
c       
c      * * * * * * * * PROBLEM NAME: First line of INPUT file * * * * * * * *
c       
c      PRBNAME is a character string describing the problem to be computed.  It can have a
c        maximum of 80 characters and must occupy the first line of the INPUT file, which
c        consists of PRBNAME and the following six namelists.
c       
c      * * * * * * * * NUMERICAL PARAMETERS: Namelist /NUMPARAM/ * * * * * * * *
c       
c      ALPHA controls the finite difference representation of the advection terms in the
c        momentum equation. 
c      AUTOT is a flag which set to zero preserves the initial time step and set nonzero 
c        causes an automatic time step adjustment during the course of a calculation.  The
c        time step for the next cycle is taken as the minimum of all the time step constraints in
c        subroutine DELTADJ.
c      CON (>0.0, >0.5) is the maximum allowable Courant number (cell width fraction 
c        traversed by the fluid in a time step) in the calculation.
c      CONSERVE is a logical flag for the momentum advection terms.  If CONSERVE=false, 
c        nonconservative finite difference expressions are used everywhere.
c      CRAY is a logical flag set to true for formatted writes to the standard output
c        with a LUN (logical unit number) of 59.  CRAY=false sets the standard output LUN
c        to 6.  A standard output LUN of 59 is typically used for computations on a CRAY
c        supercomputer. The value of CRAY determines the numerical value of IOTTY.
c      DELT (>0.0) is the initial time step of calculation.  After the calculation begins, 
c        DELT can change from its initial value if the automatic time step adjustment flag,
c        AUTOT, is nonzero.
c      DMPDT (>0.0) is the time increment between calls to subroutine TAPOUT, which 
c        writes restart dumps to TAPE7 by performing a binary write of all common blocks.
c      DTMAX (>0.0) is the maximum allowed time step.
c      ERRICCG (> machine epsilon) is the convergence criterion for the ICCG solution of 
c        the PPE.
c      FCVLIM (>0.0, < 0.5) is the Courant number above which a recomputation of the current
c        cycle is performed with a time step reduced by 20%.  Diagnostics informing the
c        user of such an action are written to the standard ouput and the ERRORS file.
c        Assumptions implicit the VOF advection algorithm break down for Courant numbers
c        greater than about 1/2.
c      FRCTN (>0.0, <1.0) is the fraction of the fluid density, RHOF, below which solutions 
c        to the Navier Stokes equations are not obtained.  Pressures and velocities in those
c        cells satisfying rho < FRCTN*RHOF are set to PSAT and zero, respectively.
c      GFNCTN is a logical flag controlling the form of the g(x) in the expression for the
c        volume force due to surface tension.  The function g(i,j) is set to 2*F(i,j) or 1
c        for GFNCTN equal to true or false, respectively.
c      IDIV is a flag set nonzero to enable a divergence correction to the VOF function
c        after advection.  The correction term is added in those cells with a nonzero 
c        divergence to counteract truncation and rounding errors.
c      ITMXICCG (>0) is the maximum number of iterations allowed in the ICCG solution of 
c        the PPE.
c      KB,KT,KL,KR are boundary condition flags for the bottom, top, left, and right
c        boundaries, respectively.  Six different boundary conditions are possible,
c        and each are associated with one of the 6 possible boundary flag values:
c
c                         1.........Rigid Free-Slip
c                         2.........Rigid No-Slip
c                         3.........Continuative Outflow
c                         4.........Periodic
c                         5.........Applied Pressure
c                         6.........Specified Inflow or Outflow
c
c      NPACK is a flag set nonzero to activate a biased packing (vertical, top-to-bottom) 
c        of the VOF function after advection.  This flag should be used with caution since
c        the packing alters the free surface without regard to the current fluid flow field.
c        The fluid is forcibly packed down until all surface cells (F(i,j)<1) lie above
c        all fluid cells.
c      NRESTART is the restart dump number.  NRESTART< 0 results in the start of a 
c        new calculation from t=0.  NRESTART>0 causes RIPPLE to restart from binary
c        dump #NRESTART on TAPE10.
c      PLTDT (>0.0) is the time increment between calls to subroutine PLOTOUT, which 
c        generates graphics if input variable PLOTS is set to true.
c      PRTDT (>0.0) is the time increment between calls to subroutine PRTPLT, which
c        writes various arrays (i.e., velocity, pressure, etc.) to the EDIT file.
c      SYM is a logical flag informing the ICCG solver of the PPE if the operator matrix 
c        to be inverted is symmetric (SYM=true) or nonsymmetric (SYM=false).  SYM needs only
c        to be false when mixed Dirichlet and Neumann boundary conditions are being enforced
c        for the pressure.
c      TWFIN (>0.0) is the problem finish time.
c       
c      * * * * * * * * FLUID PARAMETERS: Namelist /FLDPARAM/ * * * * * * * *
c       
c      CANGLEB (>0.0, <180.0) is the equilibrium contact angle (in degrees) used on 
c        the bottom boundary for wall adhesion effects.
c      CANGLET (>0.0, <180.0) is the equilibrium contact angle (in degrees) used on 
c        the top boundary for wall adhesion effects.
c      CANGLEL (>0.0, <180.0) is the equilibrium contact angle (in degrees) used on 
c        the left boundary for wall adhesion effects.
c      CANGLER (>0.0, <180.0) is the equilibrium contact angle (in degrees) used on 
c        the right boundary for wall adhesion effects.
c      GX is the body acceleration in the positive x (r) direction.
c      GY is the body acceleration in the positive y (z) direction.
c      ICYL is the problem geometry indicator, set to 0 for Cartesian coordinates
c        and 1 for cylindrical coordinates.
c      ISURF10 is the surface tension indicator, set nonzero to impose surface tension
c        on all free surfaces. 
c      PBC(n) is the pressure to be applied at a constant pressure boundary
c        (KB, KT, KL, or KR equal to 5).   The value of index n in array PBC
c        corresponds to a specific boundary: 1=bottom; 2=top; 3=left; 4=right.  
c      PSAT is the constant void pressure. Fluid pressures are scaled relative
c        to this value.
c      RHOF (>0.0) is the fluid density in the F=1.0 regions.
c      SIGMA (>0.0) is the surface tension coefficient.
c      UI is the initial fluid velocity in the positive x (r) direction.
c      UINF(n) is the flow velocity in the positive x (r) direction at a specified 
c        inflow or outflow boundary (KB, KT, KL, or KR equal to 6).   The value of 
c        index n in array UINF corresponds to a specific boundary: 1=bottom; 2=top;
c        3=left; 4=right.  
c      VI is the initial fluid velocity in the positive y (z) direction.
c      VINF(n) is the flow velocity in the positive y (z) direction at a specified 
c        inflow or outflow boundary (KB, KT, KL, or KR equal to 6).   The value of 
c        index n in array VINF corresponds to a specific boundary: 1=bottom; 2=top;
c        3=left; 4=right.  
c      XNU (>0.0) is the fluid coefficient of kinematic viscosity.
c       
c      * * * * * * * * MESH GENERATION: Namelist /MESH/ * * * * * * * *
c       
c      DXMN(m),DYMN(m) (>0.0) are the minimum space increments in the x-direction and 
c        y-direction, respectively, for submesh m.
c      NKX,NKY (>0) are the number of submesh regions in the x- and y-directions,
c        respectively. 
c      NXL(m) (>0) is the number of cells between locations XL(m) and XC(m) in
c        submesh m.
c      NYL(m) (>0) is the number of cells between locations YL(m) and YC(m) in
c        submesh m.
c      NXR(m) (>0) is the number of cells between locations XC(m) and XL(m+1)
c        in submesh m.
c      NYR(m) (>0) is the number of cells between locations YC(m) and YL(m+1)
c        in submesh m.
c      XC(m),YC(m) (>0.0) are the x- and y-coordinates of the convergence point 
c        in submesh m.
c      XL(m) (>0.0) is the location of the left edge of submesh m.  NKX+1 values 
c        of XL(m) are necessary because the right edge of submesh m is equivalent to 
c        the left edge of submesh m+1, XL(m+1).  The values should be given as an
c        increasing sequence of numbers.
c      YL(m) (>0.0) is the location of the bottom edge of submesh m.  NKY+1 values 
c        of YL(m) are necessary because the top edge of submesh m to equivalent to
c        the bottom edge of submesh m+1, YL(m+1).  The values should be given as an
c        increasing sequence of numbers.
c       
c      * * * * * * * * INTERIOR OBSTACLES: Namelist /OBSTCL/ * * * * * * * *
c
c           f(x,y)=oa1*x + oa2*x**2 + ob1*y + ob2*y**2 + oc1 + oc2*x*y
c                  + od1*cos(nxo*x) + od2*sin(mxo*x) + oe1*cos(nyo*y)
c                  + oe2*sin(myo*y)
c       
c      NOBS (>0) is the number of conic functions defining any flow obstacles interior 
c        to the computational mesh.
c      IOH(n) is the interior obstacle indicator flag on cells satisfying f(x,y)<0 
c        (for any (x,y) within the cell) for obstacle function n.  A nonzero value adds
c        obstacles in those cells; a value of 0 subtracts obstacles.
c      MXO(n) is the real number coefficient of the mesh wavenumber in the x-direction 
c        (k_x*m) for the x-direction sine term in obstacle function n.  
c      MYO(n) is the real number coefficient of the mesh wavenumber in the y-direction 
c        (k_y*m) for the y-direction sine term in obstacle function n.
c      NXO(n) is the real number coefficient of the mesh wavenumber in the x-direction 
c        (k_x*m) for the x-direction cosine term in obstacle function n.  
c      NYO(n) is the real number coefficient of the mesh wavenumber in the y-direction 
c        (k_y*m) for the y-direction cosine term in obstacle function n.  
c      OA1(n) is the coefficient of the x term in obstacle function n.
c      OA2(n) is the coefficient of the x**2 term in obstacle function n.
c      OB1(n) is the coefficient of the y term in obstacle function n.
c      OB2(n) is the coefficient of the y**2 term in obstacle function n.
c      OC1(n) is the constant term in obstacle function n.
c      OC2(n) is the coefficient of the x*y term in obstacle function n.
c      OD1(n) is the amplitude of the x-direction cosine term in obstacle function n.  
c      OD2(n) is the amplitude of the x-direction sine term in obstacle function n.  
c      OE1(n) is the amplitude of the y-direction cosine term in obstacle function n.  
c      OE2(n) is the amplitude of the y-direction sine term in obstacle function n.  
c       
c      * * * * * * * * FREE SURFACES: Namelist /FREESURF/ * * * * * * * *
c
c         f(x,y)=fa1*x + fa2*x**2 + fb1*y + fb2*y**2 + fc1 + fc2*x*y
c                + fd1*cos(nxf*x) + fd2*sin(mxf*x) + fe1*cos(nyf*y)
c                + fe2*sin(myf*y)
c
c       
c      NFRSRF (>0) is the number of conic functions used in initializing free surfaces.
c        Fluid initially occupies the entire domain if conic functions are used to specify
c        the free surface.
c      IFH(n) is the fluid indicator flag on cells satisfying f(x,y)<0 (for any (x,y)
c        within that cell) for free surface function n.  A value of 0 adds fluid in those
c        cells; a nonzero value subtracts fluid.
c      MXF(n) is the real number coefficient of the mesh wavenumber in the x-direction 
c        (k_x*m) for the x-direction sine term in free surface function n.  
c      MYF(n) is the real number coefficient of the mesh wavenumber in the y-direction 
c        (k_y*m) for the y-direction sine term in free surface function n.  
c      NXF(n) is the real number coefficient of the mesh wavenumber in the x-direction 
c        (k_x*m) for the x-direction cosine term in free surface function n.  
c      NYF(n) is the real number coefficient of the mesh wavenumber in the y-direction 
c        (k_y*m) for the y-direction cosine term in free surface function n.  
c      FA1(n) is the coefficient of the x term in free surface function n.
c      FA2(n) is the coefficient of the x**2 term in free surface function n.
c      FB1(n) is the coefficient of the y term in free surface function n.
c      FB2(n) is the coefficient of the y**2 term in free surface function n.
c      FC1(n) is the constant term in free surface function n.
c      FC2(n) is the coefficient of the x*y term in free surface function n.
c      FD1(n) is the amplitude of the x-direction cosine term in free surface function n.  
c      FD2(n) is the amplitude of the x-direction sine term in free surface function n.  
c      FE1(n) is the amplitude of the y-direction cosine term in free surface function n.  
c      FE2(n) is the amplitude of the y-direction sine term in free surface function n.  
c      DELTEMIN (>0.0) is the initial iteration time step used in the surface energy 
c        minimization routine SRFEMIN (see input variable EMIN).
c      EMIN is a logical flag set to true for allowing the user-initialized free
c        surface to change shape and evolve toward a minimal surface energy configuration 
c        it before the computation begins.  Surface energy minimization routine SRFEMIN 
c        computes the free surface evolution by iterating ITMXEMIN times with a typical 
c        computational cycle, resetting fluid velocities to zero every ITSKPEMNth iteration.
c        This option is a more general substitute for IEQUIB.
c      FLHT (>0.0, <y_max) is the vertical fluid height in the positive y (z) direction relative
c        to the y=y_min position.  This input variable applies only when NFRSRF=0.
c      IEQUIB is a flag set nonzero for initializing the fluid free surface in an
c        equilibrium meniscus configuration, subject to the contact angle CANGLER and a
c        nonzero fluid height FLHT.  This initialization requires that NFRFSRF=0, and is
c        valid only in the absence of interior obstacles along the left and right portion of
c        the problem boundaries that contain the free surface.
c      ITMXEMIN (>0) is the total number of iterations used in the surface energy 
c        minimization routine SRFEMIN (see input variable EMIN).
c      ITSKPEMN (>0) is the number of iterations to skip in the surface energy minimization 
c        routine SRFEMIN before resetting fluid velocities to zero (see input variable EMIN).
c      NSMOOTH (>0) is the number of times per computational cycle the VOF function F is 
c        smoothed when input variable SMOOTH=true.
c      SMOOTH is a logical variable set to true to use the smoothed VOF function for the
c        calculation of free surface curvature.
c      UPRIGHT is a logical variable set to false for redefining input variable
c      FLHT as the vertical fluid height in the negative y (z) direction relative 
c        to the y=y_max position.  As a result, the fluid is initialized in an 
c        inverted position. This input variable applies only when NFRSRF=0.
c       
c      * * * * * * * * GRAPHICS: Namelist /GRAPHICS/ * * * * * * * *
c      
c      DUMP is a logical flag set to true for a binary dump to TAPE8 of all the plot
c        variables.
c      ICOLOR is an integer color parameter for all graphics:
c
c                   0  -  no color
c                   1  -  color, with blue shaded where f(ij)=1.0
c                         (fluid), white where f(ij)=0.0 (voids),
c                         and combinations of blue and white where
c                         0.0 < f(ij) < 1.0 (surface cells)
c                   2  -  color, with blue shaded everywhere that
c                         f(ij) > 0.0 and white where f(ij)=0.0
c                         (voids).  The boundary between blue
c                         shade and white is the reconstructed
c                         interface (either horizontal or vertical)
c                         based on the surface cell flags in the
c                         nf(ij) array
c                   3  -  color, with blue shaded everywhere that
c                         f(ij) > 0.5 and white where f(ij) < 0.5.
c                         The boundary between blue and white
c                         shade is the f(ij)=0.5 contour.
c               All obstacles are shaded brick red.
c
c      IOUT(i) is the executive plotting array for all available plot variables.  Each plot
c        variable is associated with a specific element i of IOUT:
c                      RIPPLE Plotting Specifications
c
c               iout(1) - mesh                 iout(11) - ac
c               iout(2) - vlcty vctrs          iout(12) - vorticity
c               iout(3) - u                    iout(13) - kappa
c               iout(4) - v                    iout(14) - stream f
c               iout(5) - p                    iout(15) - volume
c               iout(6) - vof                  iout(16) - potential
c               iout(7) - Fsv*Fsv              iout(17) - tauxx
c               iout(8) - vlcty divergence     iout(18) - tauyy
c               iout(9) - ar                   iout(19) - tauxy
c               iout(10) - at                  iout(20) - vof-tilde
c
c               iout(21) -                     iout(31) - grad rho
c               iout(22) -                     iout(32) - unit wall normal
c               iout(23) -                     iout(33) - f=0.5 contour
c               iout(24) -                     iout(34) -
c               iout(25) -                     iout(35) -
c               iout(26) -                     iout(36) -
c               iout(27) -                     iout(37) -
c               iout(28) -                     iout(38) -
c               iout(29) -                     iout(39) -
c               iout(30) - srfc frc vctrs      iout(40) -
c
c        Array IOUT(i) has integer values in the range 0-7, yielding, for plot variable
c        i, the following plot combinations:
c
c                         Plotting Conventions
c                           (Values for IOUT)
c
c                       Plot Type        IOUT(n)
c                       ---------        -------
c                        no plot            0
c                        contour            1
c                          2d               2
c                          3d               3
c                       3d,contour          4
c                       2d,contour          5
c                         2d,3d             6
c                      2d,3d,contour        7
c
c      IXSKIP (>1) is the integer number of cells to skip in the x (r) direction before
c        the next cell data is plotted.  A value of IXSKIP greater than 1 displays data
c        for a given plot variable on a coarser grid that is a subset of the actual
c        computational grid.
c      IXSYMPLT is an integer flag set to 1 for plots symmetric about the x-axis.
c      IYSKIP (>1) is the integer number of cells to skip in the y (z) direction before
c        the next cell data is plotted.  A value of IYSKIP greater than 1 displays data
c        for a given plot variable on a coarser grid that is a subset of the actual
c        computational grid.
c      IYSYMPLT is an integer flag set to 1 for plots symmetric about the y-axis.
c      PLOTS is a logical flag that, if true, permits the generation of graphics during the 
c        calculation.
c      SCALE (>0.0) is a dimensionless multiplier on the length of all vectors in the 
c        vector plots.  The length of any plotted vector is a product of SCALE, the actual
c        vector magnitude, and the minimum time scale in the system (minimum mesh spacing/
c        maximum vector magnitude) at that time.
c      UNFRMMSH is a logical flag enabling a more accurate and inexpensive generation of 
c        surface plots and contours by the DISSPLA plotting routines.  It can be set to
c        true only in those calculations with a uniform mesh (delta_x=delta_y).
c      VMXFRCTN (>0.0) is the plotting threshold for all vector plots.  Vectors are not 
c        plotted  if their magnitude is below the product of VMXFRCTN and the maximum
c        vector magnitude  at that time. 
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c       
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c##############################################################
      implicit real*8 (a-h,o-z)
c      include "32bit.h"
c##############################################################
c
c############
      include "comdk1.h" 
      include "iccgdk.h" 
c############
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... establish diagnostics
c

c     ------------------
       call begin
c     ------------------
c
c.... read input and do the problem setup
c
c     ------------------
      call setup
c     ------------------
c
      if (nrestart.gt.0) go to 10
c
c.... print initial input data
c
c     ------------------
      call prtplt2 (1)
c     ------------------
c
c.... set initial boundary conditions
c
c     ------------------
      call bc
c     ------------------
c
c.... end of problem generation,
c     write out misc info
c
c     --------------------
c      call second (tbeg)
c      call second (tsetf)
c     --------------------
c      tset=tsetf-tseti
c      write (13,100) tset,t,ncyc
c      write (iotty,100) tset,t,ncyc

c.... define length scale and initiation
c
      if (itur.eq.1) then
	call length_scale
	endif
c     -------------------

c
c.... Begin cycle
c
c     ------------
   10 call newcyc
c     ------------
c
c ... calculate surface friction and drag in porous media - fyshi 10/25/02

c     -------------
        call apbp
c     -------------

c.... Compute the tilde velocities
c
c     -------------
      call vtilde
c     -------------
c
c.... Compute the volume forces due to
c     surface tension if we need them;
c     then accelerate interfacial fluid
c     elements with these forces
c
      if (isurf10.ne.0) then
c     -----------------
        call tension
        call srffrce
c     -----------------
      endif
c
c.... Directly solve the pressure Poisson equation
c
c     -------------
      call implctp
c     -------------

c      
c .... Apply sponge layer for open boundary kr=7
c      
c     -------------------
	call sponge_layer
c     -------------------

c
c.... Advect the VOF function
c
c     -------------
      call vofadv
c     -------------
c
c...  Compute surface
c
c     -------------
	call surface
c     -------------


c...  call turbulence model
c
c     -------------
	if (itur.eq.1)then

	call shear	 
	do m=1,max
	call input_klmodel(m)
	call klmodel_m(m)
	end do  
c     -------------

c
c...  update turbulence model
c
	do i=1,imax
	do j=1,jmax
	ij=(j-1)*imax+i
	do m=1,max
	ken(m,ij)=ke(m,ij)
	ke(m,ij)=0.
	end do
	end do
	end do

	endif	

      if(t.le.twfin)go to 10
     
c
  100 format (" End of problem generation at ",1pe12.5," s ",/,
     &        " Entering main loop at time ",1pe12.5," , ",
     &        " cycle ",i5)
      end

 
      subroutine begin
c
c ======================================================================
c
c   Purpose -
c     get job identification, run-time limit, date, time of day,
c     machine, and initialize graphics
c     tape5:  Input file (INPUT)
c     tape7:  Binary restart dump file
c     tape8:  Binary graphics dump file
c     tape9:  Error message output file
c     tape10:  Binary restart read file
c     tape13:  General edit output file
c    6 or 59:  Standard output
c
c
c   BEGIN is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          RIPPLE
c
c
c   BEGIN calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c            MACH    system    SECOND    system    GETJTL    system
c         SETFLSH    system     GPLOT    system   LIBDISP    system
c
c ======================================================================
c
c##############################################################
      implicit real*8 (a-h,o-z)
c      include "32bit.h"
c##############################################################
c
c############
      include "comdk1.h" 
c############
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... open the input and output files
c
c     ------------------------------------------------
      open (unit=1,file="input",form="formatted")
      open (unit=7,file="tape7",form="unformatted")
      open (unit=8,file="tape8",form="unformatted")
      open (unit=9,file="errors",form="formatted")
      open (unit=10,file="tape10",form="unformatted")
      open (unit=13,file="edit",form="formatted")
c     ------------------------------------------------
c
c.... time the problem setup
c     -------------------
c      call second(tseti)
c     -------------------
c
c.... get job time limit, date, time of day, and machine
c     --------------------
c      call getjtl(ttl)
c      call dateh(dat)
c      call timeh(tim)
c      call mach(ochn)
c     --------------------
c
c.... metafile initialization
c     ---------------------------
c      call gplot (1hu,prbname,80)
c     ---------------------------
c
c.... automatic call to gdone upon job termination
c     -----------------------
c      call setflsh
c     -----------------------
c
c.... establish disspla graphics environment
c     -----------------------
c      call libdisp
c     -----------------------
c
      return
      end
 
  
      subroutine global
c
c ======================================================================
c
c   Purpose -
c     Tally global fluid quantities such as fluid volume,
c     momentum, and kinetic energy
c
c   GLOBAL is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          PRTPLT
c
c
c   GLOBAL calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c            none
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
       include "iccgdk.h" 
c############
c
      dimension umom(nxy),vmom(nxy)
      equivalence (umom(1),sc1(1)),(vmom(1),sc2(1))
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      fvol=0.0
      vvol=0.0
      xmv=0.0
      ymv=0.0
      tke=0.0
      do 200 j=1,jm1
        do 200 i=1,im1
          ip=i+1
          jp=j+1
          ij=(j-1)*imax+i
          ipj=ij+1
          ijp=ij+imax
          rhoij=f(ij)*rhof
          rhoipj=f(ipj)*rhof
          rhoijp=f(ijp)*rhof
          rhorf=(delx(ip)*rhoij+delx(i)*rhoipj)/(delx(i)+delx(ip))
          rhotf=(dely(jp)*rhoij+dely(j)*rhoijp)/(dely(j)+dely(jp))
          umom(ij)=rhorf*u(ij)
          xmv=xmv+umom(ij)
          vmom(ij)=rhotf*v(ij)
          ymv=ymv+vmom(ij)
          tke=tke+0.5*(rhorf*u(ij)**2+rhotf*v(ij)**2)
          if (i.eq.1.or.j.eq.1) go to 200
          fvol=fvol+f(ij)*ac(ij)*(x(i)-x(i-1))*
     &          (y(j)-y(j-1))*(1.+cyl*(ri(i)*tpi-1.))
          vvol=vvol+(1.-f(ij))*ac(ij)*(x(i)-x(i-1))*
     &          (y(j)-y(j-1))*(1.+cyl*(ri(i)*tpi-1.))
  200 continue
c
      return
      end


      subroutine newcyc
c
c ======================================================================
c
c   Purpose -
c     begin cycle:  provide monitor print, optional graphics
c                   and/or field variable print, and test for run
c                   termination.  if continuing, increment time
c                   and cycle, move advance time arrays into
c                   time-n arrays, and adjust time step.
c
c   NEWCYC is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          RIPPLE
c
c
c   NEWCYC calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c          PRTPLT    ripple    SECOND    system      KILL    ripple
c          TAPOUT    ripple   DELTADJ    ripple     EXITA    system
c          TIMING    system   PLOTOUT    ripple
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
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
c.... if vof function convection limit exceeded,
c     decrease the timestep and try again
c
      if (flgc.gt.0.5) go to 100
c
c.... print time and cycle data on paper
c
c     ----------------
      call prtplt2 (2)
c     ----------------
c
      if (t+em6.lt.twplt.and.t.lt.twfin) go to 60
      twplt=twplt+pltdt
c
c.... generate graphics
c
c     ----------------
c      call plotout 		!qun
c     ----------------
c
   60 if (ncyc.le.0) go to 70
      if (t+em6.lt.twprt) go to 75
      twprt=twprt+prtdt
c
c.... print field variable data on paper
c
c     ----------------
   70 call prtplt2 (3)
c     ----------------
c
   75 if (t+em6.lt.twdmp) go to 80
      twdmp=twdmp+dmpdt
c
c.... back-up binary tape dump
c
c     ----------------
c      call tapout (1)
c     ----------------
c
c.... check to see if problem finish time surpassed
c
   80 if (t.gt.twfin) go to 130
c
c.... set the advance time arrays into the time-n arrays
c
      do 90 i=1,imax
        do 90 j=1,jmax
          ij=(j-1)*imax+i
          un(ij)=u(ij)
          vn(ij)=v(ij)
          u(ij)=0.0
          v(ij)=0.0
          pn(ij)=p(ij)
          fn(ij)=f(ij)
   90 continue
c
c.... adjust the time step (delt)
c
c     ----------------
  100 call deltadj
c     ----------------
c
c.... advance time
c
      t=t+delt
c
c.... check for vof function convection limit,
c     terminate if it has been exceeded
c
      if (nflgc.lt.100) go to 110
c
      write (9,170) ncyc,t
      write (iotty,170) ncyc,t
c
c.... exit from newcyc and terminate
c
c     --------------------------
c      call kill (iotty,ncyc,"NEWCYC")    ! qun
	stop
c     --------------------------
c
c.... check for pressure convergence failure,
c     terminate if number of failed cycles
c     (noncon) is nonzero
c
  110 if (nocon.lt.5) go to 120
c
      write (9,160) ncyc,t
      write (iotty,160) ncyc,t
c
c.... exit from newcyc and terminate
c
c     --------------------------------
c      call kill (iotty,ncyc,"NEWCYC")
	stop				! qun
c     --------------------------------
c
c.... everything ok, advance cycle
c
  120 ncyc=ncyc+1
c
c.... check cpu time used thus far   ! qun
c
c     -------------------
c      call second (time)
c     -------------------
c
c      tleft=ttl-time-tquit
c
c.... write a restart dump if out of time
c
c     ------------------------------------
c      if (tleft.lt.0.) call tapout (0)
c     ------------------------------------
c
c.... check and monitor cpu (grind) time
c     used, printing out every 25 cycles,
c     before starting a new cycle
c
c      if (mod(ncyc,25).ne.0) go to 9999
c
c     -------------------
c      call second (tend)
c      call second (tendd)
c     -------------------
c
c      grind=1000.*(tend+tend-tendd-tbeg)/float(ibar*jbar)
c      write (13,140) t,ncyc,grind,iter
c      write (iotty,140) t,ncyc,grind,iter
c      go to 9999
c
c.... normal termination (t > twfin)
c
c     --------------------------
c  130 call timing(cpu,sys,xiom)
c     --------------------------
c
c      write (iotty,180)
c      write (13,180)
c      write (iotty,190) cpu,sys,xiom
c      write (13,190) cpu,sys,xiom
c
c.... exit and destroy the dropfile
c
c     --------------
c      call exita(0)
c     --------------
c
c.... check timestep: if too small (delt < dtend),
c     terminate run on next cycle
c
c 9999 if(delt.gt.dtend) go to 20
c      write (iotty,150) ncyc,t
c      write (9,150) ncyc,t
c      twfin=t*0.999
c
c   20 iter=0
c     ------------------
c      call second (tbeg)
c     ------------------
c
 130  continue    ! qun
      return
c                    * * * * error section * * * *
  140 format (2x,"t = ",1pe12.5,"  cycle = ",i5,2x,"grind = ",1pe12.5,
     &        " ms","  iter = ",i3)
  150 format (1x,29hdelt less than dtend on cycle,i6,1x,3ht =,1pe15.7)
  160 format (//1x,25htoo many pressit failures,3x,5hncyc=,i7,1x,2ht=,1p
     1         e14.6//)
  170 format (//1x,24htoo many vofadv failures,3x,5hncyc=,i7,1x,2ht=,1pe
     1         14.6//)
  180 format (/,18x,"* * * * Normal Termination * * * *")
  190 format (/,25x,"RIPPLE Timing Statistics",/,
     &        25x,"(floating-point seconds)",/,
     &        "CPU: ",1pe12.5,2x,"System Call: ",1pe12.5,2x,
     &        "I/O + Memory: ",1pe12.5)
      end



C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C
C		SUBROUTINES
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


      subroutine deltadj
c
c ======================================================================
c
c   Purpose -
c     time step adjustment
c
c   DELTADJ is called by -
c
c          name      name      name      name      name      name
c        --------  --------  --------  --------  --------  --------
c          NEWCYC   SRFEMIN
c
c
c   DELTADJ calls the following subroutines and functions -
c
c          name      source    name      source    name      source
c        --------  --------  --------  --------  --------  --------
c            none
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
c
      data itmin /20/, itmost /200/
      data itcjr / 100 /
c
c <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c
      deltn=delt
      itc=0
      jtc=0
      if (flgc.lt.0.5) go to 20
      t=t-delt
      ncyc=ncyc-1
      delt=0.8*delt
c
      do 11 i=1,imax
        do 10 j=1,jmax
          ij=(j-1)*imax+i
          p(ij)=pn(ij)
          f(ij)=fn(ij)
          u(ij)=0.0
          v(ij)=0.0
   10   continue
   11 continue
c
      flgc=0.0
      nflgc=nflgc+1
   20 continue
      if (autot.eq.0.0) go to 40
      dumx=em10
      dvmx=em10
c
      ijx=0
      ijy=0
      do 31 i=1,im1
        do 30 j=2,jm1
          ij=(j-1)*imax+i
          udm=abs(un(ij))/(xi(i+1)-xi(i))
          vdm=abs(vn(ij))/(yj(j+1)-yj(j))
          dumx=dmax1(dumx,udm)
          dvmx=dmax1(dvmx,vdm)
          if (dumx.eq.udm) ijx = ij
          if (dvmx.eq.vdm) ijy = ij
   30   continue
   31 continue
c
      jx=ijx/imax + 1
      ix=ijx-(jx-1)*imax
      jy=ijy/imax + 1
      iy=ijy-(jy-1)*imax
c
      dtmp=1.05d0
      if(iter.lt.itmin) dtmp=1.10d0
      if(iter.gt.itmost.and.liter.gt.itmost) dtmp=0.99d0
      delto=delt*dtmp
      delto=dmin1(delto,dtmax)
c
      delt=dmin1(delto,con/dumx,con/dvmx,dtvis,dtsft)
      if (delt.eq.con/dumx) then
        itc=ix
        jtc=jx
      endif
      if (delt.eq.con/dvmx) then
        itc=iy
        jtc=jy
      endif
      if (delt.eq.dtsft) then
        itc=isft
        jtc=jsft
      endif
      if (delt.eq.dtvis) then
        itc=ivis
        jtc=jvis
      endif
   40 if (delt.eq.deltn) go to 60
c
   60 continue
      liter=iter
c.... time step growing
      if ((delt.eq.delto).and.(dtmp.gt.1.0d0)) idt="g"
c.... time step decaying
      if ((delt.eq.delto).and.(dtmp.lt.1.0d0)) idt="d"
c.... time step unchanged
      if ((delt.eq.deltn).or.(dtmp.eq.1.0d0)) idt="f"
c.... x Courant limit
      if (delt.eq.con/dumx) idt="cx"
c.... y Courant limit
      if (delt.eq.con/dvmx) idt="cy"
c.... viscous limit
      if (delt.eq.dtvis) idt="v"
c.... surface tension limit
      if (delt.eq.dtsft) idt="s"
c.... maximum allowed time step
      if (delt.eq.dtmax) idt="m"
c
      return
      end


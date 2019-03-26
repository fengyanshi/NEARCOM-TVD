C       MASTER PROGRAM
C       master.f is used to couple ADCIRC module and WAVE module.
C       It is modified based on the master program in NearCoM model package
C       The interpolation scheme used between structured and unstructured grids
C       is the same as the structured-structured except 1D arrays are used
C       for unstructured grid.
C       Fengyan Shi  05/06/2006
C
C       NOTE: in minput.dat, Nx_Circ is always set to "1" for unstructured grid 
C
C
c ------------------------------------------------------------------

        program master
        USE SIZES
        USE PASS 
        USE INTERP_COEF
        implicit none
        integer i,j,interp_coefi
        integer master_steps
        real(SZ) Time_Circ, Time_Wave, Time_Sedi, Time_Output
        data Time_Circ /0./, Time_Wave /0./, Time_Sedi /0./, 
     &       Time_Output /0./

! test

!        Master_Start=1
!        do istep =1,2
!        call WaveModule
!        Master_Start=0
!        enddo
!        stop


c --- read file and allocate pass and interp_coef

        call readfile
       
c --- calculate interpolation coefficients
       write(*,*)'Master currently read interp_coefficient.dat,'
       write(*,*)'if you want to re-calculate interp, '
       write(*,*)'modify interp_coef.dat'
       open(10,file='interp_coef.dat')
       read(10,*) interp_coefi
       close(10)
       if(interp_coefi.eq.1)then
       open(10,file='interp_coefficient.dat')
       read(10,*) Sc_01,S1_01, S2_01,S3_01, Sc_02,S1_02, S2_02,S3_02,
     &       Sc_03,S1_03, S2_03,S3_03, Sc_12,S1_12, S2_12,S3_12,
     &       Sc_13,S1_13, S2_13,S3_13, Sc_21,S1_21, S2_21,S3_21,
     &       Sc_23,S1_23, S2_23,S3_23, Sc_31,S1_31, S2_31,S3_31,
     &       Sc_32,S1_32, S2_32,S3_32,
     &       nx1_01,ny1_01,nx2_01,ny2_01,nx3_01,ny3_01,
     &       nx1_02,ny1_02,nx2_02,ny2_02,nx3_02,ny3_02,
     &       nx1_03,ny1_03,nx2_03,ny2_03,nx3_03,ny3_03,
     &       nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &       nx1_13,ny1_13,nx2_13,ny2_13,nx3_13,ny3_13,
     &       nx1_21,ny1_21,nx2_21,ny2_21,nx3_21,ny3_21,
     &       nx1_23,ny1_23,nx2_23,ny2_23,nx3_23,ny3_23,
     &       nx1_31,ny1_31,nx2_31,ny2_31,nx3_31,ny3_31,
     &       nx1_32,ny1_32,nx2_32,ny2_32,nx3_32,ny3_32
         close(10)
        else
        call get_interpolation_coef
        open(10,file='interp_coefficient.dat')
       write(10,*) Sc_01,S1_01, S2_01,S3_01, Sc_02,S1_02, S2_02,S3_02,
     &       Sc_03,S1_03, S2_03,S3_03, Sc_12,S1_12, S2_12,S3_12,
     &       Sc_13,S1_13, S2_13,S3_13, Sc_21,S1_21, S2_21,S3_21,
     &       Sc_23,S1_23, S2_23,S3_23, Sc_31,S1_31, S2_31,S3_31,
     &       Sc_32,S1_32, S2_32,S3_32,
     &       nx1_01,ny1_01,nx2_01,ny2_01,nx3_01,ny3_01,
     &       nx1_02,ny1_02,nx2_02,ny2_02,nx3_02,ny3_02,
     &       nx1_03,ny1_03,nx2_03,ny2_03,nx3_03,ny3_03,
     &       nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &       nx1_13,ny1_13,nx2_13,ny2_13,nx3_13,ny3_13,
     &       nx1_21,ny1_21,nx2_21,ny2_21,nx3_21,ny3_21,
     &       nx1_23,ny1_23,nx2_23,ny2_23,nx3_23,ny3_23,
     &       nx1_31,ny1_31,nx2_31,ny2_31,nx3_31,ny3_31,
     &       nx1_32,ny1_32,nx2_32,ny2_32,nx3_32,ny3_32
          close(10)
          endif

c --- depth interplation/extrapolation

        call interp_depth
  
c ---   module Initialization and first call
c       the first call is for hot start

         Master_Start = 1 

         call MasterInit() 

c ---    wave
         if(N_Interval_CallWave.ge.0)call WaveModule() 

            Waveupdat_for_Circ = .true.
            Waveupdat_for_Sedi = .true.

         if(Waveupdat_for_Circ) then
            call interp_wave_circ
         endif

c ---    circulation
         if(N_Interval_CallCirc.ge.0)call CircModule()

           Circupdat_for_Wave = .true.
           Circupdat_for_Sedi = .true.
c ---    sediment
         if(N_Interval_CallSedi.ge.0)call SediModule() 

c ---    set initialization switch
         Master_Start=-1

c ---   get Master_steps and ...
         
        Master_steps = Total_Time / Master_dt

        if(N_Interval_CallCirc.gt.0)then 
          nCirc = N_Interval_CallCirc
        else
          nCirc= Master_steps+1
        endif

        if(N_Interval_CallWave.gt.0)then 
          nWave = N_Interval_CallWave
        else
          nWave= Master_steps+1
        endif

        if(N_Interval_CallSedi.gt.0)then 
          nSedi = N_Interval_CallSedi
        else
          nSedi= Master_steps+1
        endif

        nOut = N_Interval_Output

c ---   Do timestepping   
        
        do 100 istep = 1, Master_steps-1

        Time_Master= istep*Master_dt

        write(*,*) 'Master Time = ',Time_Master, 's'
 
c --- call wave module

         if (N_Interval_CallWave.ge.0.and.mod(istep,nWave).eq.0) then

            if(Circupdat_for_Wave)then
               call interp_circ_wave
            endif

            if(Sediupdat_for_Wave)then
               call interp_sedi_wave
            endif

            call WaveModule

            Circupdat_for_Wave = .false.
            Sediupdat_for_Wave =.false.
            Waveupdat_for_Circ = .true.
            Waveupdat_for_Sedi = .true.

         end if 

c ---  call circulation module

         if (N_Interval_CallCirc.ge.0.and.mod(istep,nCirc).eq.0) then

         if(Waveupdat_for_Circ) then
            call interp_wave_circ
         endif

         if(Sediupdat_for_Circ) then
            call interp_sedi_circ
         endif

           call CircModule()

           Waveupdat_for_Circ = .false.
           Sediupdat_for_Circ = .false.
           Circupdat_for_Wave = .true.
           Circupdat_for_Sedi = .true.

         end if 

c ---  call sediment module

         if(N_Interval_CallSedi.ge.0
     & .and.mod(istep,nSedi).eq.0.and.istep.ge.N_Delay_CallSedi)then 

         if(Waveupdat_for_Sedi) then
            call interp_wave_sedi
         endif

         if(Circupdat_for_Sedi)then
            call interp_circ_Sedi
         endif

         call SediModule()

         Waveupdat_for_Sedi = .false.
         Circupdat_for_Sedi = .false.
         Sediupdat_for_Wave = .true.
         Sediupdat_for_Circ = .true.

         end if 

c --- output
           call interp_circ_wave

         if (N_Interval_Output.gt.0.and.mod(istep-1,nOut).eq.0) then 
           call Mexport() 
         end if 

100      continue



c   ---  program end

        print*,'Program end'

       end 
c ------------------------------------------------------------------
        subroutine readfile
        USE PASS
        USE INTERP_COEF
        implicit none
        integer i,j

        namelist /f_names/ f_depth,f_xymast,f_xywave,f_xycirc,
     &                 f_xysedi,f_name6,
     &                 f_name7,
     &                 f_name8,f_name9,f_name10,f_name11,f_name12,
     &                 f_name13,f_name14,
     &                 f_name15,f_name16

     &           /gridin/ Nx_Mast, Ny_Mast, Nx_Circ,
     &                  Ny_Circ, Nx_Wave, Ny_Wave, Nx_Sedi, Ny_Sedi,
     &                  Grid_Mast_Wave_Same, Grid_Mast_Circ_Same,
     &                  Grid_Mast_Sedi_Same,
     &                  Wave_Staggered, Circ_Staggered,Sedi_Staggered,
     &                  Wave_Stag_huv,  Circ_Stag_huv, Sedi_Stag_huv,
     &                  Wave_structured, Circ_Structured, 
     &                  Sedi_Structured,
     &                  Grid_Extrapolation

     &           /interaction/
     &                  Wave_Curr_Interact,
     &                  Wave_Bed_Interact,
     &                  Curr_Bed_Interact  

     &           /passvariables/
     &               Wave_To_Circ_Height,
     &               Wave_To_Circ_Angle,
     &               Wave_To_Circ_WaveNum,
     &               Wave_To_Circ_C,
     &               Wave_To_Circ_Radiation,
     &               Wave_To_Circ_Rad_Surf,
     &               Wave_To_Circ_Rad_Body,
     &               Wave_To_Circ_Forcing,
     &               Wave_To_Circ_MassFlux,
     &               Wave_To_Circ_Dissipation,
     &               Wave_To_Circ_BottomUV,
     &               Wave_To_Circ_Brkindex,
     &               Circ_To_Wave_UV,
     &               Circ_To_Wave_eta,
     &               Wave_To_Sedi_Height,
     &               Wave_To_Sedi_Angle,
     &               Wave_To_Sedi_BottomUV,
     &               Circ_To_Sedi_SedFlux,
     &               Circ_To_Sedi_UV,
     &               Circ_To_Sedi_UVb,
     &               Circ_To_Sedi_eta,
     &               Circ_To_Sedi_UV3D,
     &               Circ_To_Sedi_fw,
     &               Circ_To_Sedi_vt,
     &               Circ_To_Sedi_UVquasi3D,
     &               Sedi_To_Wave_Depth,
     &               Sedi_To_Circ_Depth,
     &               Circ_POM, Circ_SC  

     &           /vectorrotate/
     &                Circ_Rotate_Angle, Wave_Rotate_Angle,
     &                Sedi_Rotate_Angle   

     &           /timein/ Total_Time,Master_dt, N_Interval_CallWave,
     &                  N_Interval_CallCirc,N_Interval_CallSedi,
     &                  N_Delay_CallSedi,WaveSpectral_Data_Interval,
     &                  N_Interval_Output

c        include 'ini_logical.f'

        Grid_Wave_Circ_Same = .false.
        Grid_Wave_Sedi_Same = .false.
        Grid_Circ_Sedi_Same = .false.
        Wave_Staggered = .false.
        Circ_Staggered = .false.
        Sedi_Staggered = .false.
        Wave_Structured = .false.
        Circ_Structured = .false.
        Sedi_Structured = .false.
        Grid_Extrapolation = .false.
        Wave_Curr_Interact = .false.
        Wave_Bed_Interact = .false.
        Curr_Bed_Interact = .false.
        Wave_To_Circ_Height = .false.
        Wave_To_Circ_Angle = .false.
        Wave_To_Circ_WaveNum = .false.
        Wave_To_Circ_C = .false.
        Wave_To_Circ_Cg = .false.
        Wave_To_Circ_Radiation = .false.
        Wave_To_Circ_Rad_Surf = .false.
        Wave_To_Circ_Rad_Body = .false.
        Wave_To_Circ_Forcing = .false.
        Wave_To_Circ_MassFlux = .false.
        Wave_To_Circ_Dissipation = .false.
        Wave_To_Circ_BottomUV = .false.
        Wave_To_Circ_Brkindex = .false.
        Circ_To_Wave_UV = .false.
        Circ_To_Wave_eta = .false.
        Wave_To_Sedi_Height = .false.
        Wave_To_Sedi_Angle = .false.
        Wave_To_Sedi_BottomUV = .false.
        Circ_To_Sedi_UV = .false.
        Circ_To_Sedi_UVb = .false.
        Circ_To_Sedi_eta = .false.
        Circ_To_Sedi_UV3D = .false.
        Circ_To_Sedi_fw = .false.
        Circ_To_Sedi_vt = .false.
        Circ_To_Sedi_UVquasi3D = .false.
        Sedi_To_Wave_Depth = .false.
        Sedi_To_Circ_Depth = .false.
        Circ_POM = .false.
        Circ_SC = .false.

c -- initialize unstructured grid

        Nele_Circ=1
        Nele_Wave=1
        Nele_Sedi=1
        Nele_Mast=1

        open(1,file='minput.dat')
        
        read(1,nml=f_names)
        read(1,nml=gridin)
        read(1,nml=interaction)
        read(1,nml=passvariables)
        read(1,nml=vectorrotate)
        read(1,nml=timein)
        close(1)        

C -- allocate pass and interp_coef

        call ALLOCATE_PASS
        call ALLOCATE_INTERP_COEF

C minimal consistency check
        if(Circ_POM .and. Circ_SC)then
          write(0,*)' Specify only one circulation module'
          stop
        endif
        if(Circ_POM)then
          Wave_To_Circ_Height = .true.
          Wave_To_Circ_Angle = .true.
          Wave_To_Circ_WaveNum = .true.
          Wave_To_Circ_C = .true.
          Wave_To_Circ_Cg = .true.
          Wave_To_Circ_MassFlux = .true.
          Wave_To_Circ_Dissipation = .true.
          Wave_To_Circ_BottomUV = .true.
          Circ_Staggered = .true.
C check
          Circ_Stag_huv(1)=1
          Circ_Stag_huv(2)=2
          Circ_Stag_huv(3)=3
        endif
c put definitions for other modules here so all pass flags are set
        if(Circ_POM)then
          if(Circ_to_Sedi_vt)then
            write(0,*)'vt not implimented in POM'
            stop
          endif
          if(Circ_to_Sedi_fw)then
            write(0,*)'fw not implimented in POM'
            stop
          endif
          if(Circ_To_Sedi_UVquasi3D)then
            write(0,*)'use UV3D with POM'
            stop
          endif
        endif

c --- read initial depth

c -- for adcirc

        if(f_depth.eq.'adcirc')then
         print*,'read depth from adcirc grid file: fort.14'

        open(14,file='fort.14')
        read(14,*)
        read(14,*)Nele_Circ,Ny_Circ

        call ALLOCATE_ELEMENT

        Nx_Circ=1
        do j=1,Ny_Circ
         read(14,*)i,X_Circ(Nx_Circ,j),Y_Circ(Nx_Circ,j)
     &              ,Depth_Circ(Nx_Circ,j)
        enddo

        do j=1,Nele_Circ
         read(14,*)i,i,Nele1_Circ(j),Nele2_Circ(j),Nele3_Circ(j)
        enddo

        close(14)

        else

        call ALLOCATE_ELEMENT

        open(1,file=f_depth)
        do j=1,Ny_Mast
          read(1,*)(Depth_Mast(i,j),i=1,Nx_Mast)
        enddo
        close(1)

        endif

c -- end adcirc

c --- read xy of master grid

        open(1,file=f_xymast)
        do j=1,Ny_Mast
          read(1,*)(X_Mast(i,j),i=1,Nx_Mast)
        enddo
        do j=1,Ny_Mast
          read(1,*)(Y_Mast(i,j),i=1,Nx_Mast)
        enddo
        close(1)

c --- read xy of wave module

        if(f_xywave.ne.' ')then
        open(1,file=f_xywave)
        do j=1,Ny_Wave
          read(1,*)(X_Wave(i,j),i=1,Nx_Wave)
        enddo
        do j=1,Ny_Wave
          read(1,*)(Y_Wave(i,j),i=1,Nx_Wave)
        enddo
        close(1)
        else
        do j=1,Ny_Wave
        do i=1,Nx_Wave
          X_Wave(i,j)=X_Mast(i,j)
          Y_Wave(i,j)=Y_Mast(i,j)
        enddo
        enddo
        Grid_Mast_Wave_Same = .true.
        endif

c --- read xy of circulation module
       
c -- for adcirc

       if(f_depth.eq.'adcirc')then
         write(*,*) 'because of adcirc, skip reading XY_Circ'        
       else

        if(f_xycirc.ne.' ')then
        open(1,file=f_xycirc)
        do j=1,Ny_Circ
          read(1,*)(X_Circ(i,j),i=1,Nx_Circ)
        enddo
        do j=1,Ny_Circ
          read(1,*)(Y_Circ(i,j),i=1,Nx_Circ)
        enddo
        close(1)
        else
        do j=1,Ny_Circ
        do i=1,Nx_Circ
          X_Circ(i,j)=X_Mast(i,j)
          Y_Circ(i,j)=Y_Mast(i,j)
        enddo
        enddo
        Grid_Mast_Circ_Same = .true.
        endif
        endif

c --- read xy of sediment module

        if(f_xysedi.ne.' ')then
        open(1,file=f_xysedi)
        do j=1,Ny_Sedi
          read(1,*)(X_Sedi(i,j),i=1,Nx_Sedi)
        enddo
        do j=1,Ny_Sedi
          read(1,*)(Y_Sedi(i,j),i=1,Nx_Sedi)
        enddo
        close(1)
        else
        do j=1,Ny_Sedi
        do i=1,Nx_Sedi
          X_Sedi(i,j)=X_Mast(i,j)
          Y_Sedi(i,j)=Y_Mast(i,j)
        enddo
        enddo
        Grid_Mast_Sedi_Same = .true.
        endif

c --- grid relations between wave-circ, wave-sedi, and circ-sedi

        Grid_Wave_Circ_Same = .false.
        Grid_Wave_Sedi_Same = .false.
        Grid_Circ_Sedi_Same = .false.

        if(Grid_Mast_Wave_Same.and.Grid_Mast_Circ_Same)
     &     Grid_Wave_Circ_Same = .true.

        if(Grid_Mast_Wave_Same.and.Grid_Mast_Sedi_Same)
     &     Grid_Wave_Sedi_Same = .true.

        if(Grid_Mast_Circ_Same.and.Grid_Mast_Sedi_Same)
     &     Grid_Circ_Sedi_Same = .true.

100      format(800f16.8)

        return
        end
c --------------------------------------------------------------

        subroutine get_interpolation_coef
        USE PASS
        USE INTERP_COEF
        implicit none
        integer i,j

       if(f_depth.eq.'adcirc')then

          call interpolation
     .    (Nx_Circ,Ny_Circ,X_Circ,Y_Circ,
     .     Nele_Circ,Nele1_Circ,Nele2_Circ,Nele3_Circ,
     .     Nx_Mast,Ny_Mast,X_Mast,Y_Mast,Sc_20,S1_20,S2_20,S3_20,
     .     nx1_20,ny1_20,nx2_20,ny2_20,nx3_20,ny3_20)

          write(*,*) 'get interp_coef circ-mast because of using adcirc'

       endif

c ---   Mast-Wave
        if(Grid_Mast_Wave_Same) then
          call interpsame(Nx_Wave,Ny_Wave,Sc_01,S1_01,S2_01,S3_01,
     .          nx1_01,ny1_01,nx2_01,ny2_01,nx3_01,ny3_01)
        else
         write(*,*)'Grid_Mast & Grid_Wave are different, calc coef...'
          call interpolation
     .    (Nx_Mast,Ny_Mast,X_Mast,Y_Mast,
     .     Nele_Mast,Nele1_Mast,Nele2_Mast,Nele3_Mast,
     .     Nx_Wave,Ny_Wave,X_Wave,Y_Wave,Sc_01,S1_01,S2_01,S3_01,
     .     nx1_01,ny1_01,nx2_01,ny2_01,nx3_01,ny3_01)
        endif

c ---   Mast-Circ
        if(Grid_Mast_Circ_Same) then
          call interpsame(Nx_Circ,Ny_Circ,Sc_02,S1_02,S2_02,S3_02,
     .          nx1_02,ny1_02,nx2_02,ny2_02,nx3_02,ny3_02)
        else
         write(*,*)'Grid_Mast & Grid_Circ are different, calc coef...'

          call interpolation
     .    (Nx_Mast,Ny_Mast,X_Mast,Y_Mast,
     .     Nele_Mast,Nele1_Mast,Nele2_Mast,Nele3_Mast,
     .     Nx_Circ,Ny_Circ,X_Circ,Y_Circ,Sc_02,S1_02,S2_02,S3_02,
     .     nx1_02,ny1_02,nx2_02,ny2_02,nx3_02,ny3_02)       
        endif

c ---   Mast-Sedi
        if(Grid_Mast_Sedi_Same) then
          call interpsame(Nx_Sedi,Ny_Sedi,Sc_03,S1_03,S2_03,S3_03,
     .          nx1_03,ny1_03,nx2_03,ny2_03,nx3_03,ny3_03)
        else
         write(*,*)'Grid_Mast & Grid_Sedi are different, calc coef...'
          call interpolation
     .    (Nx_Mast,Ny_Mast,X_Mast,Y_Mast,
     .     Nele_Mast,Nele1_Mast,Nele2_Mast,Nele3_Mast,
     .     Nx_Sedi,Ny_Sedi,X_Sedi,Y_Sedi,Sc_03,S1_03,S2_03,S3_03,
     .     nx1_03,ny1_03,nx2_03,ny2_03,nx3_03,ny3_03)
        endif

c ---   Circ-Wave
        if(Grid_Wave_Circ_Same) then
          call interpsame(Nx_Wave,Ny_Wave,Sc_21,S1_21,S2_21,S3_21,
     .          nx1_21,ny1_21,nx2_21,ny2_21,nx3_21,ny3_21)
        else
         write(*,*)'Grid_Circ & Grid_Wave are different, calc coef...'
          call interpolation
     .    (Nx_Circ,Ny_Circ,X_Circ,Y_Circ,
     .     Nele_Circ,Nele1_Circ,Nele2_Circ,Nele3_Circ,
     .     Nx_Wave,Ny_Wave,X_Wave,Y_Wave,Sc_21,S1_21,S2_21,S3_21,
     .     nx1_21,ny1_21,nx2_21,ny2_21,nx3_21,ny3_21)
        endif

c ---   Sedi-Wave
        if(Grid_Wave_Sedi_Same) then
          call interpsame(Nx_Wave,Ny_Wave,Sc_31,S1_31,S2_31,S3_31,
     .          nx1_31,ny1_31,nx2_31,ny2_31,nx3_31,ny3_31)
        else
         write(*,*)'Grid_Sedi & Grid_Wave are different, calc coef...'
          call interpolation
     .    (Nx_Sedi,Ny_Sedi,X_Sedi,Y_Sedi,
     .     Nele_Sedi,Nele1_Sedi,Nele2_Sedi,Nele3_Sedi,
     .     Nx_Wave,Ny_Wave,X_Wave,Y_Wave,Sc_31,S1_31,S2_31,S3_31,
     .     nx1_31,ny1_31,nx2_31,ny2_31,nx3_31,ny3_31)
        endif

c ---   Wave-Circ
        if(Grid_Wave_Circ_Same) then
          call interpsame(Nx_Circ,Ny_Circ,Sc_12,S1_12,S2_12,S3_12,
     .          nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12)
        else
         write(*,*)'Grid_Wave & Grid_Circ are different, calc coef...'
          call interpolation
     .    (Nx_Wave,Ny_Wave,X_Wave,Y_Wave,
     .     Nele_Wave,Nele1_Wave,Nele2_Wave,Nele3_Wave,
     .     Nx_Circ,Ny_Circ,X_Circ,Y_Circ,Sc_12,S1_12,S2_12,S3_12,
     .     nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12)
        endif

c ---   Sedi-Circ
        if(Grid_Circ_Sedi_Same) then
          call interpsame(Nx_Circ,Ny_Circ,Sc_32,S1_32,S2_32,S3_32,
     .          nx1_32,ny1_32,nx2_32,ny2_32,nx3_32,ny3_32)
        else
         write(*,*)'Grid_Sedi & Grid_Circ are different, calc coef...'
          call interpolation
     .    (Nx_Sedi,Ny_Sedi,X_Sedi,Y_Sedi,
     .     Nele_Sedi,Nele1_Sedi,Nele2_Sedi,Nele3_Sedi,
     .     Nx_Circ,Ny_Circ,X_Circ,Y_Circ,Sc_32,S1_32,S2_32,S3_32,
     .     nx1_32,ny1_32,nx2_32,ny2_32,nx3_32,ny3_32)
        endif

c ---   Wave-Sedi
        if(Grid_Wave_Sedi_Same) then
          call interpsame(Nx_Sedi,Ny_Sedi,Sc_13,S1_13,S2_13,S3_13,
     .          nx1_13,ny1_13,nx2_13,ny2_13,nx3_13,ny3_13)
        else
         write(*,*)'Grid_Wave & Grid_Sedi are different, calc coef...'
          call interpolation
     .    (Nx_Wave,Ny_Wave,X_Wave,Y_Wave,
     .     Nele_Wave,Nele1_Wave,Nele2_Wave,Nele3_Wave,
     .     Nx_Sedi,Ny_Sedi,X_Sedi,Y_Sedi,Sc_13,S1_13,S2_13,S3_13,
     .     nx1_13,ny1_13,nx2_13,ny2_13,nx3_13,ny3_13)
        endif

c ---   Circ-Sedi
        if(Grid_Circ_Sedi_Same) then
          call interpsame(Nx_Sedi,Ny_Sedi,Sc_23,S1_23,S2_23,S3_23,
     .          nx1_23,ny1_23,nx2_23,ny2_23,nx3_23,ny3_23)
        else
         write(*,*)'Grid_Circ & Grid_Sedi are different, calc coef...'
          call interpolation
     .    (Nx_Circ,Ny_Circ,X_Circ,Y_Circ,
     .     Nele_Circ,Nele1_Circ,Nele2_Circ,Nele3_Circ,
     .     Nx_Sedi,Ny_Sedi,X_Sedi,Y_Sedi,Sc_23,S1_23,S2_23,S3_23,
     .     nx1_23,ny1_23,nx2_23,ny2_23,nx3_23,ny3_23)
        endif

100     continue

        return
        end
c --------------------------------------------------------------

        subroutine interp_depth
        USE PASS
        USE INTERP_COEF
        implicit none
        integer i,j

       if(f_depth.eq.'adcirc')then

        call grid1_to_grid2(Nx_Circ,Ny_Circ,Nx_Wave,Ny_Wave,
     &                      Sc_21,S1_21,S2_21,S3_21,
     &                      nx1_21,ny1_21,nx2_21,ny2_21,nx3_21,ny3_21,
     &                      Depth_Circ,Depth_Wave,
     &                      .false.,0,Wave_Staggered,Wave_Stag_huv(1))

        call grid1_to_grid2(Nx_Circ,Ny_Circ,Nx_Mast,Ny_Mast,
     &                      Sc_20,S1_20,S2_20,S3_20,
     &                      nx1_20,ny1_20,nx2_20,ny2_20,nx3_20,ny3_20,
     &                      Depth_Circ,Depth_Mast,
     &                      .false.,0,Circ_Staggered,Circ_Stag_huv(1))

        call grid1_to_grid2(Nx_Circ,Ny_Circ,Nx_Sedi,Ny_Sedi,
     &                      Sc_23,S1_23,S2_23,S3_23,
     &                      nx1_23,ny1_23,nx2_23,ny2_23,nx3_23,ny3_23,
     &                      Depth_Circ,Depth_Sedi,
     &                      .false.,0,Sedi_Staggered,Sedi_Stag_huv(1))


        else

        call grid1_to_grid2(Nx_Mast,Ny_Mast,Nx_Wave,Ny_Wave,
     &                      Sc_01,S1_01,S2_01,S3_01,
     &                      nx1_01,ny1_01,nx2_01,ny2_01,nx3_01,ny3_01,
     &                      Depth_Mast,Depth_Wave,
     &                      .false.,0,Wave_Staggered,Wave_Stag_huv(1))

        call grid1_to_grid2(Nx_Mast,Ny_Mast,Nx_Circ,Ny_Circ,
     &                      Sc_02,S1_02,S2_02,S3_02,
     &                      nx1_02,ny1_02,nx2_02,ny2_02,nx3_02,ny3_02,
     &                      Depth_Mast,Depth_Circ,
     &                      .false.,0,Circ_Staggered,Circ_Stag_huv(1))

        call grid1_to_grid2(Nx_Mast,Ny_Mast,Nx_Sedi,Ny_Sedi,
     &                      Sc_03,S1_03,S2_03,S3_03,
     &                      nx1_03,ny1_03,nx2_03,ny2_03,nx3_03,ny3_03,
     &                      Depth_Mast,Depth_Sedi,
     &                      .false.,0,Sedi_Staggered,Sedi_Stag_huv(1))

        endif


c --- test interpolatoin
c        call output(Nx_Mast,Ny_Mast,1,Depth_Mast)
c        call output(Nx_Circ,Ny_Circ,2,Depth_Circ)

        return
        end
c --------------------------------------------------------------

        subroutine interp_circ_wave
        USE SIZES
        USE PASS
        USE INTERP_COEF
        implicit none
        real(SZ) Tmp1, Tmp2, tht
        integer i,j

        if (Wave_Curr_Interact) then

          if(Circ_To_Wave_UV) then


        call grid1_to_grid2(Nx_Circ,Ny_Circ,Nx_Wave,Ny_Wave,
     &                      Sc_21,S1_21,S2_21,S3_21,
     &                      nx1_21,ny1_21,nx2_21,ny2_21,nx3_21,ny3_21,
     &                      Pass_U,Intp_U_Wave,
     &                      Circ_Staggered,Circ_Stag_huv(2),
     &                      Wave_Staggered,Wave_Stag_huv(2))

        call grid1_to_grid2(Nx_Circ,Ny_Circ,Nx_Wave,Ny_Wave,
     &                      Sc_21,S1_21,S2_21,S3_21,
     &                      nx1_21,ny1_21,nx2_21,ny2_21,nx3_21,ny3_21,
     &                      Pass_V,Intp_V_Wave,
     &                      Circ_Staggered,Circ_Stag_huv(3),
     &                      Wave_Staggered,Wave_Stag_huv(3))


             if((Circ_Rotate_Angle-Wave_Rotate_Angle).ne.0)then
               do j=1,Ny_Wave
               do i=1,Nx_Wave
                 Tmp1=Intp_U_Wave(i,j)
                 Tmp2=Intp_V_Wave(i,j)
                 tht=(Circ_Rotate_Angle-Wave_Rotate_Angle)*3.14159/180.
                  Intp_U_Wave(i,j)=Tmp1*cos(tht)- Tmp2*sin(tht)
                  Intp_V_Wave(i,j)=Tmp1*sin(tht)+ Tmp2*cos(tht)
               enddo
               enddo
        
             endif
          endif

          endif

          if(Circ_To_Wave_eta)then

        call grid1_to_grid2(Nx_Circ,Ny_Circ,Nx_Wave,Ny_Wave,
     &                      Sc_21,S1_21,S2_21,S3_21,
     &                      nx1_21,ny1_21,nx2_21,ny2_21,nx3_21,ny3_21,
     &                      Pass_eta,Intp_eta_Wave,
     &                      Circ_Staggered,Circ_Stag_huv(1),
     &                      Wave_Staggered,Wave_Stag_huv(1))    
         
          endif


        return
        end     
c --------------------------------------------------------------

        subroutine interp_sedi_wave
        USE SIZES
        USE PASS
        USE INTERP_COEF
        implicit none
        real(SZ) Tmp1,Tmp2,tht
        integer i,j


        if (Wave_Bed_Interact) then
          if(Sedi_To_Wave_Depth)then
        call grid1_to_grid2(Nx_Sedi,Ny_Sedi,Nx_Wave,Ny_Wave,
     &                      Sc_31,S1_31,S2_31,S3_31,
     &                      nx1_31,ny1_31,nx2_31,ny2_31,nx3_31,ny3_31,
     &                      Pass_Dupdated,Depth_Wave,
     &                      Sedi_Staggered,Sedi_Stag_huv(1),
     &                      Wave_Staggered,Sedi_Stag_huv(1))
          endif
        endif

        return
        end     
c --------------------------------------------------------------

        subroutine interp_wave_circ
        USE SIZES
        USE PASS
        USE INTERP_COEF
        implicit none
        real(SZ) Tmp1, Tmp2, tht
        integer i,j

c --- wave height

        if(Wave_To_Circ_Height)then

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Circ,Ny_Circ,
     &                      Sc_12,S1_12,S2_12,S3_12,
     &                      nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &                      Pass_Height,Intp_Height_Circ,
     &                      Wave_Staggered,0,
     &                      Circ_Staggered,0)

        endif

c --- wave angle

        if(Wave_To_Circ_Angle)then

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Circ,Ny_Circ,
     &                      Sc_12,S1_12,S2_12,S3_12,
     &                      nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &                      Pass_Theta,Intp_Theta_Circ,
     &                      Wave_Staggered,0,
     &                      Circ_Staggered,0)

             if((Wave_Rotate_Angle-Circ_Rotate_Angle).ne.0)then
               do j=1,Ny_Circ
               do i=1,Nx_Circ
                 tht=Wave_Rotate_Angle-Circ_Rotate_Angle
                 Intp_Theta_Circ(i,j)=Intp_Theta_Circ(i,j)+tht
               enddo
               enddo       
             endif

        endif


c --- wave number

        if(Wave_To_Circ_WaveNum)then

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Circ,Ny_Circ,
     &                      Sc_12,S1_12,S2_12,S3_12,
     &                      nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &                      Pass_WaveNum,Intp_WaveNum_Circ,
     &                      Wave_Staggered,0,
     &                      Circ_Staggered,0)

        endif

c --- wave C

        if(Wave_To_Circ_C)then

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Circ,Ny_Circ,
     &                      Sc_12,S1_12,S2_12,S3_12,
     &                      nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &                      Pass_C,Intp_C_Circ,
     &                      Wave_Staggered,0,
     &                      Circ_Staggered,0)

        endif

c --- radiation stress (depth-average form)

        if(Wave_To_Circ_Radiation)then

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Circ,Ny_Circ,
     &                      Sc_12,S1_12,S2_12,S3_12,
     &                      nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &                      Pass_Sxx,Intp_Sxx_Circ,
     &                      Wave_Staggered,0,
     &                      Circ_Staggered,0)

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Circ,Ny_Circ,
     &                      Sc_12,S1_12,S2_12,S3_12,
     &                      nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &                      Pass_Sxy,Intp_Sxy_Circ,
     &                      Wave_Staggered,0,Circ_Staggered,0) 

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Circ,Ny_Circ,
     &                      Sc_12,S1_12,S2_12,S3_12,
     &                      nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &                      Pass_Syy,Intp_Syy_Circ,
     &                      Wave_Staggered,0,Circ_Staggered,0) 

        endif

c --- radiation stress - surface

        if(Wave_To_Circ_Rad_Surf)then

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Circ,Ny_Circ,
     &                      Sc_12,S1_12,S2_12,S3_12,
     &                      nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &                      Pass_Sxx_surf,Intp_Sxx_surf,
     &                      Wave_Staggered,0,Circ_Staggered,0)

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Circ,Ny_Circ,
     &                      Sc_12,S1_12,S2_12,S3_12,
     &                      nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &                      Pass_Sxy_surf,Intp_Sxy_surf,
     &                      Wave_Staggered,0,Circ_Staggered,0) 

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Circ,Ny_Circ,
     &                      Sc_12,S1_12,S2_12,S3_12,
     &                      nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &                      Pass_Syy_surf,Intp_Syy_surf,
     &                      Wave_Staggered,0,Circ_Staggered,0) 

        endif

c --- radiation stress - body

        if(Wave_To_Circ_Rad_Body)then

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Circ,Ny_Circ,
     &                      Sc_12,S1_12,S2_12,S3_12,
     &                      nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &                      Pass_Sxx_body,Intp_Sxx_body,
     &                      Wave_Staggered,0,Circ_Staggered,0)

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Circ,Ny_Circ,
     &                      Sc_12,S1_12,S2_12,S3_12,
     &                      nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &                      Pass_Sxy_body,Intp_Sxy_body,
     &                      Wave_Staggered,0,Circ_Staggered,0) 

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Circ,Ny_Circ,
     &                      Sc_12,S1_12,S2_12,S3_12,
     &                      nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &                      Pass_Syy_body,Intp_Syy_body,
     &                      Wave_Staggered,0,Circ_Staggered,0) 

        endif


c --- short wave forcing

        if(Wave_To_Circ_Forcing)then

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Circ,Ny_Circ,
     &                      Sc_12,S1_12,S2_12,S3_12,
     &                      nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &                      Pass_Wave_Fx,Intp_Fx_Circ,
     &                      Wave_Staggered,Wave_Stag_huv(2),
     &                      Circ_Staggered,Circ_Stag_huv(2)) 

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Circ,Ny_Circ,
     &                      Sc_12,S1_12,S2_12,S3_12,
     &                      nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &                      Pass_Wave_Fy,Intp_Fy_Circ,
     &                      Wave_Staggered,Wave_Stag_huv(3),
     &                      Circ_Staggered,Circ_Stag_huv(3)) 

             if((Wave_Rotate_Angle-Circ_Rotate_Angle).ne.0)then
               do j=1,Ny_Circ
               do i=1,Nx_Circ
                 Tmp1=Intp_Fx_Circ(i,j)
                 Tmp2=Intp_Fy_Circ(i,j)
                 tht=(Wave_Rotate_Angle-CIrc_Rotate_Angle)*3.14159/180.
                  Intp_Fx_Circ(i,j)=Tmp1*cos(tht)- Tmp2*sin(tht)
                  Intp_Fy_Circ(i,j)=Tmp1*sin(tht)+ Tmp2*cos(tht)
               enddo
               enddo       
             endif

        endif

c ---  mass flux

        if(Wave_To_Circ_MassFlux) then

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Circ,Ny_Circ,
     &                      Sc_12,S1_12,S2_12,S3_12,
     &                      nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &                      Pass_MassFluxU,Intp_MassFluxU_Circ,
     &                      Wave_Staggered,Wave_Stag_huv(2),
     &                      Circ_Staggered,Circ_Stag_huv(2)) 

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Circ,Ny_Circ,
     &                      Sc_12,S1_12,S2_12,S3_12,
     &                      nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &                      Pass_MassFluxV,Intp_MassFluxV_Circ,
     &                      Wave_Staggered,Wave_Stag_huv(3),
     &                      Circ_Staggered,Circ_Stag_huv(3)) 

             if((Wave_Rotate_Angle-Circ_Rotate_Angle).ne.0)then
               do j=1,Ny_Circ
               do i=1,Nx_Circ
                 Tmp1=Intp_MassFluxU_Circ(i,j)
                 Tmp2=Intp_MassFluxV_Circ(i,j)
                 tht=(Wave_Rotate_Angle-CIrc_Rotate_Angle)*3.14159/180.
                  Intp_MassFluxU_Circ(i,j)=Tmp1*cos(tht)- Tmp2*sin(tht)
                  Intp_MassFluxV_Circ(i,j)=Tmp1*sin(tht)+ Tmp2*cos(tht)
               enddo
               enddo       
             endif

        endif

c --- wave dissipation

        if(Wave_To_Circ_Dissipation) then

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Circ,Ny_Circ,
     &                      Sc_12,S1_12,S2_12,S3_12,
     &                      nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &                      Pass_Diss,Intp_Diss_Circ,
     &                      Wave_Staggered,0,Circ_Staggered,0) 

        endif

c --- wave bottom velocity

        if(Wave_To_Circ_BottomUV) then

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Circ,Ny_Circ,
     &                      Sc_12,S1_12,S2_12,S3_12,
     &                      nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &                      Pass_ubott,Intp_ubott_Circ,
     &                      Wave_Staggered,0,Circ_Staggered,0) 

        endif

c ---  break index

        if(Wave_To_Circ_Dissipation) then

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Circ,Ny_Circ,
     &                      Sc_12,S1_12,S2_12,S3_12,
     &                      nx1_12,ny1_12,nx2_12,ny2_12,nx3_12,ny3_12,
     &                      Pass_ibrk,Intp_ibrk_Circ,
     &                      Wave_Staggered,0,Circ_Staggered,0) 

        endif

        return
        end     
c --------------------------------------------------------------

        subroutine interp_sedi_circ
        USE SIZES
        USE PASS
        USE INTERP_COEF
        implicit none
        real(SZ) Tmp1,Tmp2,tht
        integer i,j


        if (Curr_Bed_Interact) then

        if(Sedi_To_Circ_Depth)then

        call grid1_to_grid2(Nx_Sedi,Ny_Sedi,Nx_Circ,Ny_Circ,
     &                      Sc_32,S1_32,S2_32,S3_32,
     &                      nx1_32,ny1_32,nx2_32,ny2_32,nx3_32,ny3_32,
     &                      Pass_Dupdated,Depth_Circ,
     &                      Sedi_Staggered,Sedi_Stag_huv(1),
     &                      Circ_Staggered,Circ_Stag_huv(1)) 
        endif

        endif

        return
        end     
c --------------------------------------------------------------

        subroutine interp_wave_sedi
        USE SIZES
        USE PASS
        USE INTERP_COEF
        implicit none
        real(SZ) Tmp1,Tmp2,tht
        integer i,j

        if(Wave_To_Sedi_Height)then

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Sedi,Ny_Sedi,
     &                      Sc_13,S1_13,S2_13,S3_13,
     &                      nx1_13,ny1_13,nx2_13,ny2_13,nx3_13,ny3_13,
     &                      Pass_Height,Intp_Height_Sedi,
     &                      Wave_Staggered,0,Sedi_Staggered,0) 

        endif

        if(Wave_To_Sedi_Angle)then

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Sedi,Ny_Sedi,
     &                      Sc_13,S1_13,S2_13,S3_13,
     &                      nx1_13,ny1_13,nx2_13,ny2_13,nx3_13,ny3_13,
     &                      Pass_Theta,Intp_Theta_Sedi,
     &                      Wave_Staggered,0,Sedi_Staggered,0) 

             if((Wave_Rotate_Angle-Sedi_Rotate_Angle).ne.0)then
               do j=1,Ny_Sedi
               do i=1,Nx_Sedi
                 tht=Wave_Rotate_Angle-Sedi_Rotate_Angle
                 Intp_Theta_Sedi(i,j)=Intp_Theta_Sedi(i,j)+tht
               enddo
               enddo       
             endif

        endif

        if(Wave_To_Sedi_BottomUV)then

        call grid1_to_grid2(Nx_Wave,Ny_Wave,Nx_Sedi,Ny_Sedi,
     &                      Sc_13,S1_13,S2_13,S3_13,
     &                      nx1_13,ny1_13,nx2_13,ny2_13,nx3_13,ny3_13,
     &                      Pass_ubott,Intp_ubott_Sedi,
     &                      Wave_Staggered,0,Sedi_Staggered,0) 


        endif

        return
        end     
c --------------------------------------------------------------

        subroutine interp_circ_sedi
        USE SIZES
        USE PASS
        USE INTERP_COEF
        implicit none
        real(SZ) Tmp1,Tmp2,tht
        integer i,j


c --- Sediment Flux
        if(Circ_To_Sedi_SedFlux)then

        call grid1_to_grid2(Nx_Circ,Ny_Circ,Nx_Sedi,Ny_Sedi,
     &                      Sc_23,S1_23,S2_23,S3_23,
     &                      nx1_23,ny1_23,nx2_23,ny2_23,nx3_23,ny3_23,
     &                      Pass_SedFluxCum_x,Intp_SedFluxCum_x,
     &                      Circ_Staggered,Circ_Stag_huv(1),
     &                      Sedi_Staggered,Sedi_Stag_huv(1)) 

        call grid1_to_grid2(Nx_Circ,Ny_Circ,Nx_Sedi,Ny_Sedi,
     &                      Sc_23,S1_23,S2_23,S3_23,
     &                      nx1_23,ny1_23,nx2_23,ny2_23,nx3_23,ny3_23,
     &                      Pass_SedFluxCum_y,Intp_SedFluxCum_y,
     &                      Circ_Staggered,Circ_Stag_huv(1),
     &                      Sedi_Staggered,Sedi_Stag_huv(1)) 


             if((Circ_Rotate_Angle-Sedi_Rotate_Angle).ne.0)then
               do j=1,Ny_Sedi
               do i=1,Nx_Sedi
                 Tmp1=Intp_SedFluxCum_x(i,j)
                 Tmp2=Intp_SedFluxCum_y(i,j)
                 tht=(Circ_Rotate_Angle-Sedi_Rotate_Angle)*3.14159/180.
                  Intp_SedFluxCum_x(i,j)=Tmp1*cos(tht)- Tmp2*sin(tht)
                  Intp_SedFluxCum_y(i,j)=Tmp1*sin(tht)+ Tmp2*cos(tht)
               enddo
               enddo       
             endif
        endif

c --- UV
        if(Circ_To_Sedi_UV)then

        call grid1_to_grid2(Nx_Circ,Ny_Circ,Nx_Sedi,Ny_Sedi,
     &                      Sc_23,S1_23,S2_23,S3_23,
     &                      nx1_23,ny1_23,nx2_23,ny2_23,nx3_23,ny3_23,
     &                      Pass_U,Intp_U_Sedi,
     &                      Circ_Staggered,Circ_Stag_huv(2),
     &                      Sedi_Staggered,Sedi_Stag_huv(2)) 


        call grid1_to_grid2(Nx_Circ,Ny_Circ,Nx_Sedi,Ny_Sedi,
     &                      Sc_23,S1_23,S2_23,S3_23,
     &                      nx1_23,ny1_23,nx2_23,ny2_23,nx3_23,ny3_23,
     &                      Pass_V,Intp_V_Sedi,
     &                      Circ_Staggered,Circ_Stag_huv(3),
     &                      Sedi_Staggered,Sedi_Stag_huv(3)) 

             if((Circ_Rotate_Angle-Sedi_Rotate_Angle).ne.0)then
               do j=1,Ny_Sedi
               do i=1,Nx_Sedi
                 Tmp1=Intp_U_Sedi(i,j)
                 Tmp2=Intp_V_Sedi(i,j)
                 tht=(Circ_Rotate_Angle-Sedi_Rotate_Angle)*3.14159/180.
                  Intp_U_Sedi(i,j)=Tmp1*cos(tht)- Tmp2*sin(tht)
                  Intp_V_Sedi(i,j)=Tmp1*sin(tht)+ Tmp2*cos(tht)
               enddo
               enddo       
             endif
        endif

c -- Ub Vb
        if(Circ_To_Sedi_UVb)then

        call grid1_to_grid2(Nx_Circ,Ny_Circ,Nx_Sedi,Ny_Sedi,
     &                      Sc_23,S1_23,S2_23,S3_23,
     &                      nx1_23,ny1_23,nx2_23,ny2_23,nx3_23,ny3_23,
     &                      Pass_Ub,Intp_Ub_Sedi,
     &                      Circ_Staggered,Circ_Stag_huv(2),
     &                      Sedi_Staggered,Sedi_Stag_huv(2)) 


        call grid1_to_grid2(Nx_Circ,Ny_Circ,Nx_Sedi,Ny_Sedi,
     &                      Sc_23,S1_23,S2_23,S3_23,
     &                      nx1_23,ny1_23,nx2_23,ny2_23,nx3_23,ny3_23,
     &                      Pass_Vb,Intp_Vb_Sedi,
     &                      Circ_Staggered,Circ_Stag_huv(3),
     &                      Sedi_Staggered,Sedi_Stag_huv(3)) 

             if((Circ_Rotate_Angle-Sedi_Rotate_Angle).ne.0)then
               do j=1,Ny_Sedi
               do i=1,Nx_Sedi
                 Tmp1=Intp_Ub_Sedi(i,j)
                 Tmp2=Intp_Vb_Sedi(i,j)
                 tht=(Circ_Rotate_Angle-Sedi_Rotate_Angle)*3.14159/180.
                  Intp_Ub_Sedi(i,j)=Tmp1*cos(tht)- Tmp2*sin(tht)
                  Intp_Vb_Sedi(i,j)=Tmp1*sin(tht)+ Tmp2*cos(tht)
               enddo
               enddo       
             endif
        endif

c --- eta
        
        if(Circ_To_Sedi_eta)then

        call grid1_to_grid2(Nx_Circ,Ny_Circ,Nx_Sedi,Ny_Sedi,
     &                      Sc_23,S1_23,S2_23,S3_23,
     &                      nx1_23,ny1_23,nx2_23,ny2_23,nx3_23,ny3_23,
     &                      Pass_eta,Intp_eta_Sedi,
     &                      Circ_Staggered,Circ_Stag_huv(1),
     &                      Sedi_Staggered,Sedi_Stag_huv(1)) 

        endif

c ---  fw

        if(Circ_To_Sedi_fw)then

        call grid1_to_grid2(Nx_Circ,Ny_Circ,Nx_Sedi,Ny_Sedi,
     &                      Sc_23,S1_23,S2_23,S3_23,
     &                      nx1_23,ny1_23,nx2_23,ny2_23,nx3_23,ny3_23,
     &                      Pass_fw,Intp_fw_Sedi,
     &                      Circ_Staggered,0,Sedi_Staggered,0) 

        endif

c --- vt

        if(Circ_To_Sedi_vt)then

        call grid1_to_grid2(Nx_Circ,Ny_Circ,Nx_Sedi,Ny_Sedi,
     &                      Sc_23,S1_23,S2_23,S3_23,
     &                      nx1_23,ny1_23,nx2_23,ny2_23,nx3_23,ny3_23,
     &                      Pass_vt,Intp_vt_Sedi,
     &                      Circ_Staggered,0,Sedi_Staggered,0) 

        endif


        if(Circ_To_Sedi_UV3D)then
        print*,'not done yet'
        
        endif


        if(Circ_To_Sedi_UVquasi3D)then
        print*,'not done yet'
        
        endif



        return
        end     


c ----------------------------------------------------------------

        subroutine interpsame
     .  (m_grid2,n_grid2,Sc,S1,S2,S3,
     .   nx1,ny1,nx2,ny2,nx3,ny3)
        USE SIZES
        USE PASS
        USE INTERP_COEF
        implicit none
        integer i,j
! ---   m(x direction) and n(y direction)
        integer m_grid2,n_grid2


! ---   (i,j) of three points surrounded 

! ---  areas of the four triangles, area will be negative
!      if an order is clockwise
!       Sc -- triangle 1,2,3
!       S1 -- triangle 2,3,c
!       S2 -- triangle 3,1,c
!       S3 -- triangle 1,2,c

       integer  nx1(m_grid2,n_grid2)
     .         ,ny1(m_grid2,n_grid2)
     .         ,nx2(m_grid2,n_grid2)
     .         ,ny2(m_grid2,n_grid2)
     .         ,nx3(m_grid2,n_grid2)
     .         ,ny3(m_grid2,n_grid2)  

      real(SZ) Sc(m_grid2,n_grid2)
     .      ,S1(m_grid2,n_grid2)
     .      ,S2(m_grid2,n_grid2)
     .      ,S3(m_grid2,n_grid2) 


        do j=1,n_grid2
        do i=1,m_grid2
          Sc(i,j)=1.
          S1(i,j)=1.
          S2(i,j)=0.
          S3(i,j)=0.
          nx1(i,j)=i
          ny1(i,j)=j
          nx2(i,j)=1
          ny2(i,j)=1
          nx3(i,j)=1
          ny3(i,j)=1
        enddo
       enddo

        return
        end


! -----------------------------------------------------------------
        subroutine interpolation
     .  (m_grid1,n_grid1,x_grid1,y_grid1,
     .   n_ele,n1_ele,n2_ele,n3_ele,
     .   m_grid2,n_grid2,x_grid2,y_grid2,Sc,S1,S2,S3,
     .   nx1,ny1,nx2,ny2,nx3,ny3)

        USE SIZES
        USE PASS
        implicit none
        integer i,j
! ---   m(x direction) and n(y direction)
        integer m_grid1,n_grid1,m_grid2,n_grid2


! ---   types, interpolation -- 0, extrapolation -- 1
        integer ntype(m_grid2,n_grid2)

! ---  areas of the four triangles, area will be negative
!      if an order is clockwise
!       Sc -- triangle 1,2,3
!       S1 -- triangle 2,3,c
!       S2 -- triangle 3,1,c
!       S3 -- triangle 1,2,c

! ---   x, y and variables of grid1 and grid2

        real(SZ) x1,y1,x2,y2,x3,y3,area1,area2,area3,area,dist,
     .       dist_init

        integer ii,jj,nx_near,ny_near

       integer nx1(m_grid2,n_grid2)
     .         ,ny1(m_grid2,n_grid2)
     .         ,nx2(m_grid2,n_grid2)
     .         ,ny2(m_grid2,n_grid2)
     .         ,nx3(m_grid2,n_grid2)
     .         ,ny3(m_grid2,n_grid2)
       integer n_ele, n1_ele(n_ele),n2_ele(n_ele),n3_ele(n_ele)  

      real(SZ) Sc(m_grid2,n_grid2)
     .      ,S1(m_grid2,n_grid2)
     .      ,S2(m_grid2,n_grid2)
     .      ,S3(m_grid2,n_grid2)


      real(SZ) x_grid1(m_grid1,n_grid1),y_grid1(m_grid1,n_grid1)
     .      ,var_grid1(m_grid1,n_grid1)
     .      ,x_grid2(m_grid2,n_grid2),y_grid2(m_grid2,n_grid2)
     .      ,var_grid2(m_grid2,n_grid2)



! ---  control parameter, the initial -- 0
        integer Iconv
        data Iconv /0/

! --- for 1-D case

        if(n_grid2.eq.1.and.n_grid1.eq.1)then
          
          j=1
          do i=1,m_grid2
            x1=x_grid2(i,j)
            jj=1
            do ii=1,m_grid1-1
            x2=x_grid1(ii,jj)
            x3=x_grid1(ii+1,jj)
            area1=x1-x2
            area2=x3-x1
            if(area1.ge.0.and.area2.ge.0)then
            nx1(i,j)=ii+1
            nx2(i,j)=ii
            S1(i,j)=area1
            S2(i,j)=area2
            Sc(i,j)=area1+area2
            ntype(i,j)=0
            goto 1100
            endif
            enddo ! ii

            ntype(i,j)=1

1100        continue
            
          enddo ! i

          j=1
          do i=1,m_grid2
          if(ntype(i,j).eq.1)then
            jj=1
            x1=x_grid1(1,jj)
            x2=x_grid1(m_grid1,jj)
            x3=x_grid2(i,j)
            if(abs(x3-x1).lt.abs(x2-x3))then
              nx1(i,j)=2
              nx2(i,j)=1
              S1(i,j)=x3-x1
              S2(i,j)=x_grid1(2,jj)-x3
              Sc(i,j)=x_grid1(2,jj)-x1
            else
              nx1(i,j)=m_grid1-1
              nx2(i,j)=m_grid1
              S1(i,j)=x2-x3
              S2(i,j)=x3-x_grid1(m_grid1-1,jj)
              S3(i,j)=x2-x_grid1(m_grid1-1,jj)
            endif

           endif ! ntype=1

          enddo ! i

        endif ! end of 1D



c -- unstructured grid
        if(m_grid1.eq.1) then
! --- find the triangle includes the points in grid2

        do j=1,n_grid2
        do i=1,m_grid2
          x1=x_grid2(i,j)
          y1=y_grid2(i,j)

          ii=1
          do jj=1,n_ele

            x2=x_grid1(ii,n2_ele(jj))
            y2=y_grid1(ii,n2_ele(jj))
            x3=x_grid1(ii,n3_ele(jj))
            y3=y_grid1(ii,n3_ele(jj))
            area1=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            x2=x_grid1(ii,n3_ele(jj))
            y2=y_grid1(ii,n3_ele(jj))
            x3=x_grid1(ii,n1_ele(jj))
            y3=y_grid1(ii,n1_ele(jj))
            area2=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)


            x2=x_grid1(ii,n1_ele(jj))
            y2=y_grid1(ii,n1_ele(jj))
            x3=x_grid1(ii,n2_ele(jj))
            y3=y_grid1(ii,n2_ele(jj))
            area3=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            if(area1.ge.0.and.area2.ge.0.and.area3.ge.0)then
              ntype(i,j)=0
              nx1(i,j)=ii
              ny1(i,j)=n1_ele(jj)
              nx2(i,j)=ii
              ny2(i,j)=n2_ele(jj)
              nx3(i,j)=ii
              ny3(i,j)=n3_ele(jj)
              S1(i,j)=area1
              S2(i,j)=area2
              S3(i,j)=area3

              x1=x_grid1(nx1(i,j),ny1(i,j))
              y1=y_grid1(nx1(i,j),ny1(i,j))
              x2=x_grid1(nx2(i,j),ny2(i,j))
              y2=y_grid1(nx2(i,j),ny2(i,j))
              x3=x_grid1(nx3(i,j),ny3(i,j))
              y3=y_grid1(nx3(i,j),ny3(i,j))
              Sc(i,j)=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

              goto 110
            endif

            ntype(i,j)=1

          enddo

110       continue

        enddo
        enddo


! --- find the nearest point in grid1 for grid2-points with ntype=1
!     these points will be used for extrapolation

        do j=1,n_grid2
        do i=1,m_grid2

          if (ntype(i,j).eq.1)then

!           -- find the nearest point

            x1=x_grid2(i,j)
            y1=y_grid2(i,j)
            x2=1./3.*(x_grid1(1,n1_ele(1))+x_grid1(1,n2_ele(1))
     &          +x_grid1(1,n3_ele(1)))
            y2=1./3.*(y_grid1(1,n1_ele(1))+y_grid1(1,n2_ele(1))
     &          +y_grid1(1,n3_ele(1)))
            dist_init=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)
            nx_near=1
            ny_near=1

            ii=1
            do jj=2,n_ele
            x2=1./3.*(x_grid1(ii,n1_ele(jj))+x_grid1(ii,n2_ele(jj))
     &          +x_grid1(ii,n3_ele(jj)))
            y2=1./3.*(y_grid1(ii,n1_ele(jj))+y_grid1(ii,n2_ele(jj))
     &          +y_grid1(ii,n3_ele(jj)))
              dist=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)
              if(dist.lt.dist_init)then
                dist_init=dist
                nx_near=ii
                ny_near=jj
              endif
            enddo

              nx1(i,j)=ii
              ny1(i,j)=n1_ele(ny_near)
              nx2(i,j)=ii
              ny2(i,j)=n2_ele(ny_near)
              nx3(i,j)=ii
              ny3(i,j)=n3_ele(ny_near)

c -- no exterpolation for unstructured grid

            if (Grid_Extrapolation) then
              x2=x_grid1(nx2(i,j),ny2(i,j))
              y2=y_grid1(nx2(i,j),ny2(i,j))
              x3=x_grid1(nx3(i,j),ny3(i,j))
              y3=y_grid1(nx3(i,j),ny3(i,j))
              S1(i,j)=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

              x2=x_grid1(nx3(i,j),ny3(i,j))
              y2=y_grid1(nx3(i,j),ny3(i,j))
              x3=x_grid1(nx1(i,j),ny1(i,j))
              y3=y_grid1(nx1(i,j),ny1(i,j))
              S2(i,j)=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

              x2=x_grid1(nx1(i,j),ny1(i,j))
              y2=y_grid1(nx1(i,j),ny1(i,j))
              x3=x_grid1(nx2(i,j),ny2(i,j))
              y3=y_grid1(nx2(i,j),ny2(i,j))
              S3(i,j)=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

              x1=x_grid1(nx1(i,j),ny1(i,j))
              y1=y_grid1(nx1(i,j),ny1(i,j))
              x2=x_grid1(nx2(i,j),ny2(i,j))
              y2=y_grid1(nx2(i,j),ny2(i,j))
              x3=x_grid1(nx3(i,j),ny3(i,j))
              y3=y_grid1(nx3(i,j),ny3(i,j))
              Sc(i,j)=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            else

              S1(i,j)=1.
              S2(i,j)=0.
              S3(i,j)=0.
              Sc(i,j)=1.

            endif


          endif
        enddo
        enddo

        endif ! end of unstructured grid


c -- 2D structured grid
        if(m_grid1.ne.1.and.n_grid1.ne.1)then
! --- find the triangle includes the points in grid2

        do j=1,n_grid2
        do i=1,m_grid2
          x1=x_grid2(i,j)
          y1=y_grid2(i,j)
          do jj=1,n_grid1-1
          do ii=1,m_grid1-1

            x2=x_grid1(ii+1,jj)
            y2=y_grid1(ii+1,jj)
            x3=x_grid1(ii,jj+1)
            y3=y_grid1(ii,jj+1)
            area1=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            x2=x_grid1(ii,jj+1)
            y2=y_grid1(ii,jj+1)
            x3=x_grid1(ii,jj)
            y3=y_grid1(ii,jj)
            area2=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            x2=x_grid1(ii,jj)
            y2=y_grid1(ii,jj)
            x3=x_grid1(ii+1,jj)
            y3=y_grid1(ii+1,jj)
            area3=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            if(area1.ge.0.and.area2.ge.0.and.area3.ge.0)then
              ntype(i,j)=0
              nx1(i,j)=ii
              ny1(i,j)=jj
              nx2(i,j)=ii+1
              ny2(i,j)=jj
              nx3(i,j)=ii
              ny3(i,j)=jj+1
              S1(i,j)=area1
              S2(i,j)=area2
              S3(i,j)=area3
               
              x1=x_grid1(nx1(i,j),ny1(i,j))
              y1=y_grid1(nx1(i,j),ny1(i,j))
              x2=x_grid1(nx2(i,j),ny2(i,j))
              y2=y_grid1(nx2(i,j),ny2(i,j))
              x3=x_grid1(nx3(i,j),ny3(i,j))
              y3=y_grid1(nx3(i,j),ny3(i,j))
              Sc(i,j)=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

              goto 120
            endif

            x2=x_grid1(ii+1,jj)
            y2=y_grid1(ii+1,jj)
            x3=x_grid1(ii+1,jj+1)
            y3=y_grid1(ii+1,jj+1)
            area1=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            x2=x_grid1(ii+1,jj+1)
            y2=y_grid1(ii+1,jj+1)
            x3=x_grid1(ii,jj+1)
            y3=y_grid1(ii,jj+1)
            area2=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            x2=x_grid1(ii,jj+1)
            y2=y_grid1(ii,jj+1)
            x3=x_grid1(ii+1,jj)
            y3=y_grid1(ii+1,jj)
            area3=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            if(area1.ge.0.and.area2.ge.0.and.area3.ge.0)then
              ntype(i,j)=0
              nx1(i,j)=ii
              ny1(i,j)=jj+1
              nx2(i,j)=ii+1
              ny2(i,j)=jj
              nx3(i,j)=ii+1
              ny3(i,j)=jj+1
              S1(i,j)=area1
              S2(i,j)=area2
              S3(i,j)=area3
               
              x1=x_grid1(nx1(i,j),ny1(i,j))
              y1=y_grid1(nx1(i,j),ny1(i,j))
              x2=x_grid1(nx2(i,j),ny2(i,j))
              y2=y_grid1(nx2(i,j),ny2(i,j))
              x3=x_grid1(nx3(i,j),ny3(i,j))
              y3=y_grid1(nx3(i,j),ny3(i,j))
              Sc(i,j)=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

              goto 120
            endif

            ntype(i,j)=1

          enddo
          enddo

120       continue

        enddo
        enddo


! --- find the nearest point in grid1 for grid2-points with ntype=1
!     these points will be used for extrapolation

        do j=1,n_grid2
        do i=1,m_grid2

          if (ntype(i,j).eq.1)then

!           -- find the nearest point

            x1=x_grid2(i,j)
            y1=y_grid2(i,j)
            x2=x_grid1(1,1)
            y2=y_grid1(1,1)         
            dist_init=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)
            nx_near=1
            ny_near=1

            do jj=2,n_grid1-1
            do ii=2,m_grid1-1       
              x2=x_grid1(ii,jj)
              y2=y_grid1(ii,jj)
              dist=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)
              if(dist.lt.dist_init)then
                dist_init=dist
                nx_near=ii
                ny_near=jj
              endif
            enddo
            enddo

!           -- calculate four areas -- S1, S2, S3, Sc
!              choose the nearest triangle by using the sign of the area

            x2=x_grid1(nx_near+1,ny_near)
            y2=y_grid1(nx_near+1,ny_near)
            x3=x_grid1(nx_near,ny_near+1)
            y3=y_grid1(nx_near,ny_near+1)
            area=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            if(area.ge.0)then

              nx1(i,j)=nx_near
              ny1(i,j)=ny_near
              nx2(i,j)=nx_near+1
              ny2(i,j)=ny_near
              nx3(i,j)=nx_near
              ny3(i,j)=ny_near+1

            else

              nx1(i,j)=nx_near
              ny1(i,j)=ny_near+1
              nx2(i,j)=nx_near+1
              ny2(i,j)=ny_near
              nx3(i,j)=nx_near+1
              ny3(i,j)=ny_near+1

            endif
           
! --- if no extrapolation is allowed , evaluated variable will equal 
!     to the variable at neast grid point

            if (Grid_Extrapolation) then
              x2=x_grid1(nx2(i,j),ny2(i,j))
              y2=y_grid1(nx2(i,j),ny2(i,j))
              x3=x_grid1(nx3(i,j),ny3(i,j))
              y3=y_grid1(nx3(i,j),ny3(i,j))
              S1(i,j)=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

              x2=x_grid1(nx3(i,j),ny3(i,j))
              y2=y_grid1(nx3(i,j),ny3(i,j))
              x3=x_grid1(nx1(i,j),ny1(i,j))
              y3=y_grid1(nx1(i,j),ny1(i,j))
              S2(i,j)=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

              x2=x_grid1(nx1(i,j),ny1(i,j))
              y2=y_grid1(nx1(i,j),ny1(i,j))
              x3=x_grid1(nx2(i,j),ny2(i,j))
              y3=y_grid1(nx2(i,j),ny2(i,j))
              S3(i,j)=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

              x1=x_grid1(nx1(i,j),ny1(i,j))
              y1=y_grid1(nx1(i,j),ny1(i,j))
              x2=x_grid1(nx2(i,j),ny2(i,j))
              y2=y_grid1(nx2(i,j),ny2(i,j))
              x3=x_grid1(nx3(i,j),ny3(i,j))
              y3=y_grid1(nx3(i,j),ny3(i,j))
              Sc(i,j)=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

            else

              S1(i,j)=1.
              S2(i,j)=0.
              S3(i,j)=0.
              Sc(i,j)=1.

            endif
            
          endif
        enddo
        enddo

        endif ! for n_grid ?= 1

        return
        end


! -----------------------------------------------------------------
        subroutine interpolation_nonstruc
     .  (m_grid1,n_grid1,x_grid1,y_grid1,
     .   m_grid2,n_grid2,x_grid2,y_grid2,Sc,S1,S2,S3,
     .   nx1,ny1,nx2,ny2,nx3,ny3)
        USE SIZES
        USE PASS
        implicit none
        integer i,j
! ---   m(x direction) and n(y direction)
        integer m_grid1,n_grid1,m_grid2,n_grid2

! ---   (i,j) of three points surrounded 

! ---  areas of the four triangles, area will be negative
!      if an order is clockwise
!       Sc -- triangle 1,2,3
!       S1 -- triangle 2,3,c
!       S2 -- triangle 3,1,c
!       S3 -- triangle 1,2,c

! ---   x, y and variables of grid1 and grid2

        real(SZ) x1,y1,x2,y2,x3,y3,area1,area2,area3,area,dist,
     .       dist_init,dist_1,dist_2,dist_3

        integer ii,jj,nx_near,ny_near,nx_near_1,ny_near_1,
     .     nx_near_2,ny_near_2,nx_near_3,ny_near_3

! ---  control parameter, the initial -- 0
        Integer Iconv
        data Iconv /0/

       integer nx1(m_grid2,n_grid2)
     .         ,ny1(m_grid2,n_grid2)
     .         ,nx2(m_grid2,n_grid2)
     .         ,ny2(m_grid2,n_grid2)
     .         ,nx3(m_grid2,n_grid2)
     .         ,ny3(m_grid2,n_grid2) 

      real(SZ) Sc(m_grid2,n_grid2)
     .      ,S1(m_grid2,n_grid2)
     .      ,S2(m_grid2,n_grid2)
     .      ,S3(m_grid2,n_grid2) 


      real(SZ) x_grid1(m_grid1,n_grid1),y_grid1(m_grid1,n_grid1)
     .      ,var_grid1(m_grid1,n_grid1)
     .      ,x_grid2(m_grid2,n_grid2),y_grid2(m_grid2,n_grid2)
     .      ,var_grid2(m_grid2,n_grid2)


! --- find the nearest three points on grid1 
!     these points will be used for interpolation/extrapolation

        do j=1,n_grid2
        do i=1,m_grid2

!  --- find the farest point first

            x1=x_grid2(i,j)
            y1=y_grid2(i,j)
            x2=x_grid1(1,1)
            y2=y_grid1(1,1)         
            dist_1=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)
            nx_near_1=1
            ny_near_1=1

            do jj=1,n_grid1
            do ii=1,m_grid1
              x2=x_grid1(ii,jj)
              y2=y_grid1(ii,jj)
              dist=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)
              if(dist.gt.dist_1)then
                dist_1=dist
                nx_near_1=ii
                ny_near_1=jj
              endif
            enddo
            enddo

            nx_near_2=nx_near_1
            ny_near_2=ny_near_1
            nx_near_3=nx_near_1
            ny_near_3=ny_near_1
            dist_2=dist_1
            dist_3=dist_1

! ---   find nearest three points

            do jj=1,n_grid1
            do ii=1,m_grid1
              x2=x_grid1(ii,jj)
              y2=y_grid1(ii,jj)
              dist=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)
              if(dist.lt.dist_1)then
                dist_1=dist
                nx_near_1=ii
                ny_near_1=jj
              elseif(dist.lt.dist_2)then
                dist_2=dist
                nx_near_2=ii
                ny_near_2=jj
              elseif(dist.lt.dist_3)then
                dist_3=dist
                nx_near_3=ii
                ny_near_3=jj
              endif
            enddo
            enddo

!           -- calculate four areas -- S1, S2, S3, Sc


              nx1(i,j)=nx_near_1
              ny1(i,j)=ny_near_1
              nx2(i,j)=nx_near_2
              ny2(i,j)=ny_near_2
              nx3(i,j)=nx_near_3
              ny3(i,j)=ny_near_3

           
              x2=x_grid1(nx2(i,j),ny2(i,j))
              y2=y_grid1(nx2(i,j),ny2(i,j))
              x3=x_grid1(nx3(i,j),ny3(i,j))
              y3=y_grid1(nx3(i,j),ny3(i,j))
              S1(i,j)=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

              x2=x_grid1(nx3(i,j),ny3(i,j))
              y2=y_grid1(nx3(i,j),ny3(i,j))
              x3=x_grid1(nx1(i,j),ny1(i,j))
              y3=y_grid1(nx1(i,j),ny1(i,j))
              S2(i,j)=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

              x2=x_grid1(nx1(i,j),ny1(i,j))
              y2=y_grid1(nx1(i,j),ny1(i,j))
              x3=x_grid1(nx2(i,j),ny2(i,j))
              y3=y_grid1(nx2(i,j),ny2(i,j))
              S3(i,j)=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)

              x1=x_grid1(nx1(i,j),ny1(i,j))
              y1=y_grid1(nx1(i,j),ny1(i,j))
              x2=x_grid1(nx2(i,j),ny2(i,j))
              y2=y_grid1(nx2(i,j),ny2(i,j))
              x3=x_grid1(nx3(i,j),ny3(i,j))
              y3=y_grid1(nx3(i,j),ny3(i,j))
              Sc(i,j)=0.5*(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)
            
        enddo
        enddo

        print*,'two grids are different, calculate interp coef..'

        return
        end


! --------------------------------------------------------------

        subroutine grid1_to_grid2
     .      (m_grid1,n_grid1,m_grid2,n_grid2,
     .       Sc,S1,S2,S3,nx1,ny1,nx2,ny2,nx3,ny3,
     .       var_grid1,var_grid2,grid1_stag,ntype_grid1,
     .       grid2_stag,ntype_grid2)

        USE SIZES
        USE PASS
        implicit none
        integer i,j

        integer ntype_grid1,ntype_grid2

! ---   x, y and variables of grid1 and grid2

        real(SZ) tmpb

! ---   logical parameters for staggered grid
        logical grid1_stag, grid2_stag

! --- others
        integer m_grid1,n_grid1,m_grid2,n_grid2


        real(SZ)  Sc(m_grid2,n_grid2)
     .      ,S1(m_grid2,n_grid2)
     .      ,S2(m_grid2,n_grid2)
     .      ,S3(m_grid2,n_grid2)

        integer nx1(m_grid2,n_grid2),ny1(m_grid2,n_grid2)
     .      ,nx2(m_grid2,n_grid2),ny2(m_grid2,n_grid2)
     .      ,nx3(m_grid2,n_grid2),ny3(m_grid2,n_grid2)

        real(SZ) var_grid1(m_grid1,n_grid1)
     .      ,var_grid2(m_grid2,n_grid2)
     .      ,tmp(m_grid1,n_grid1),tmpa(m_grid1,n_grid1)


        do j=1,n_grid1
        do i=1,m_grid1
          tmp(i,j)=var_grid1(i,j)
        enddo
        enddo


! ---   for staggered grid1

        if (grid1_stag.and.ntype_grid1.ne.0)then

          if(n_grid1.eq.1)then
            write(*,*)'stagered type defined wrong for 1D grid'
            stop
          endif
! ---      ntype=1

          if(ntype_grid1.eq.1)then

          if(n_grid1.eq.1)then
            write(*,*)'stagered type defined wrong for 1D grid'
            stop
          endif

            do j=2,n_grid1-1
            do i=2,m_grid1-1

             tmp(i,j)=0.25*(var_grid1(i-1,j-1)+var_grid1(i,j-1)
     &                      +var_grid1(i,j)+var_grid1(i-1,j))

            enddo
            enddo

            do i=2,m_grid1-1
              tmp(i,1)=2.*var_grid1(i,2)-var_grid1(i,3)
              tmp(i,n_grid1)=2.*var_grid1(i,n_grid1-1)
     &                             -var_grid1(i,n_grid1-2)
            enddo
            do j=1,n_grid1
              tmp(1,j)=2.*var_grid1(2,j)-var_grid1(3,j)
              tmp(m_grid1,j)=2.*var_grid1(m_grid1-1,j)
     &                             -var_grid1(m_grid1-2,j)
            enddo
          endif
! ---      ntype=2
          if(ntype_grid1.eq.2)then

          if(n_grid1.eq.1)then
            write(*,*)'stagered type defined wrong for 1D grid'
            stop
          endif

            do j=2,n_grid1-1
            do i=1,m_grid1
             tmp(i,j)=0.5*(var_grid1(i,j)+var_grid1(i,j-1))
            enddo
            enddo
            
            do i=1,m_grid1
              tmp(i,1)=0.5*(3.*var_grid1(i,1)-var_grid1(i,2))
              tmp(i,n_grid1)=0.5*(3.*var_grid1(i,n_grid1-1)
     &                             -var_grid1(i,n_grid1-2))
            enddo
          endif
! ---      ntype=3
          if(ntype_grid1.eq.3)then
            do j=1,n_grid1
            do i=2,m_grid1-1
             tmp(i,j)=0.5*(var_grid1(i,j)+var_grid1(i-1,j))
            enddo
            enddo
            
            do j=1,n_grid1
              tmp(1,j)=0.5*(3.*var_grid1(1,j)-var_grid1(2,j))
              tmp(m_grid1,j)=0.5*(3.*var_grid1(m_grid1-1,j)
     &                             -var_grid1(m_grid1-2,j))
            enddo
          endif

        endif

! --- interpolation/extrapolation       

        if(n_grid1.eq.1)then

        do j=1,n_grid2
        do i=1,m_grid2
          var_grid2(i,j)=1./Sc(i,1)*(S1(i,1)*tmp(nx1(i,1),1)
     &                 +S2(i,1)*tmp(nx2(i,1),1))
        enddo
        enddo

        elseif(n_grid2.eq.1)then

        do j=1,n_grid1
        do i=1,m_grid1
         tmpa(i,j)=tmp(i,j)
        enddo
        enddo

        do i=1,m_grid1
          tmpb=0.
          do j=1,n_grid1
           tmpb=tmpb+tmpa(i,j)
          enddo
           tmp(i,1)=tmpb/n_grid1
        enddo

        j=1
        do i=1,m_grid2
          var_grid2(i,j)=1./Sc(i,j)*(S1(i,j)*tmp(nx1(i,1),1)
     &                 +S2(i,j)*tmp(nx2(i,1),1))
        enddo

        else

        do j=1,n_grid2
        do i=1,m_grid2

            var_grid2(i,j)=(S1(i,j)*tmp(nx1(i,j),ny1(i,j))
     .                   +S2(i,j)*tmp(nx2(i,j),ny2(i,j))
     .                   +S3(i,j)*tmp(nx3(i,j),ny3(i,j)))
     .                   /Sc(i,j)

        enddo
        enddo

        endif              ! ngrid.eq.1

! ---  for staggered grid2

        if (grid2_stag.and.ntype_grid2.ne.0)then
        do j=1,n_grid2
        do i=1,m_grid2
          tmp(i,j)=var_grid2(i,j)
        enddo
        enddo

! ---  ntype_grid2 = 1
          if(ntype_grid2.eq.1)then
            if(n_grid2.eq.1)then
              write(*,*)'wrong staggered grid type for 1-D case'
              stop
            endif

            do j=1,n_grid2-1
            do i=1,m_grid2-1
              var_grid2(i,j)=0.25*(tmp(i,j)+tmp(i+1,j)
     &                       +tmp(i+1,j+1)+tmp(i,j+1))
            enddo
            enddo
          endif
! ---  ntype_grid2 = 2
          if(ntype_grid2.eq.2)then
            if(n_grid2.eq.1)then
              write(*,*)'wrong staggered grid type for 1-D case'
              stop
            endif


            do j=1,n_grid2-1
            do i=1,m_grid2
              var_grid2(i,j)=0.5*(tmp(i,j)+tmp(i,j+1))
            enddo
            enddo
          endif
! ---  ntype_grid2 = 3
          if(ntype_grid2.eq.3)then
            do j=1,n_grid2
            do i=1,m_grid2-1
              var_grid2(i,j)=0.5*(tmp(i,j)+tmp(i+1,j))
            enddo
            enddo
          endif

        endif

        return
        end

! ---------------------------------------------------------------
        subroutine output(mp,np,num_file,varb)
        USE SIZES
        USE PASS
        implicit none
        integer i,j
        character*2 file_name
        integer nm_first,nm_second,nm_third,nm_fourth,np,mp,
     &          num_file

       real(SZ) varb(mp,np)

        nm_first=mod(num_file/1000,10)
        nm_second=mod(num_file/100,10)
        nm_third=mod(num_file/10,10)
        nm_fourth=mod(num_file,10)

        write(file_name(1:1),'(I1)')nm_third
        write(file_name(2:2),'(I1)')nm_fourth
c        write(file_name(3:3),'(I1)')nm_third
c        write(file_name(4:4),'(I1)')nm_fourth

        open(2,file='data'//file_name//'.dat')
        do j=1,np
        write(2,100)(varb(i,j),i=1,mp)
100     format(801f16.8)
        enddo
        close(2)

        return

        end 
c ------------------------------------------------------------------
        subroutine SediModule_sample() 
        USE PASS
        implicit none         
        integer i,j

        if(Master_Start.eq.1)then

          print*,'Sediment module initialization ...'          
        else
         print*, 'call Sediment module ...' 
        endif

        return
        end 
c ------------------------------------------------------------------

        subroutine  Mexport() 
        USE PASS
        USE INTERP_COEF
        implicit none
        integer i,j
        integer nm_first,nm_second,nm_third,nm_fourth,num_file
        character*4 file_name

        num_file=(istep-1)/N_Interval_Output
        print*,'mexport routine  plot num=', num_file

        nm_first=mod(num_file/1000,10)
        nm_second=mod(num_file/100,10)
        nm_third=mod(num_file/10,10)
        nm_fourth=mod(num_file,10)


        write(file_name(1:1),'(I1)')nm_first
        write(file_name(2:2),'(I1)')nm_second
        write(file_name(3:3),'(I1)')nm_third
        write(file_name(4:4),'(I1)')nm_fourth

c        if(f_name11.ne.' ')then
c        open(2,file=f_name11(1:2)//file_name//'.out')
c          do j=1,Ny_Circ
c            write(2,111)(Pass_Theta(i,j),i=1,Nx_Circ)
c          enddo
c        close(2)
c        endif

        if(f_name11.ne.' ')then
        open(2,file=f_name11(1:2)//file_name//'.out')
          do j=1,Ny_Wave
            write(2,111)(Intp_eta_Wave(i,j),i=1,Nx_Wave)
          enddo
        close(2)
        endif


        if(f_name12.ne.' ')then
        open(2,file=f_name12(1:2)//file_name//'.out')
          do j=1,Ny_Circ
            write(2,111)(Intp_Height_Circ(i,j),i=1,Nx_Circ)
          enddo
        close(2)
        endif

        if(f_name13.ne.' ')then
        open(2,file=f_name13(1:2)//file_name//'.out')
          do j=1,Ny_Circ
            write(2,111)(Intp_Theta_Circ(i,j),i=1,Nx_Circ)
          enddo
        close(2)
        endif

        if(f_name14.ne.' ')then
        open(2,file=f_name14(1:2)//file_name//'.out')
          do j=1,Ny_Circ
            write(2,111)(Pass_U(i,j),i=1,Nx_Circ)
          enddo
        close(2)
        endif

        if(f_name15.ne.' ')then
        open(2,file=f_name15(1:2)//file_name//'.out')
          do j=1,Ny_Circ
            write(2,111)(Pass_V(i,j),i=1,Nx_Circ)
          enddo
        close(2)
        endif

        if(f_name16.ne.' ')then
        open(2,file=f_name16(1:2)//file_name//'.out')
          do j=1,Ny_Circ
            write(2,111)(Pass_eta(i,j),i=1,Nx_Circ)
          enddo
        close(2)
        endif

        if(f_name9.ne.' ')then
        open(2,file=f_name9(1:2)//file_name//'.out')
          do j=1,Ny_Circ
            write(2,111)(Intp_Fx_Circ(i,j),i=1,Nx_Circ)
          enddo
        close(2)
        endif

        if(f_name10.ne.' ')then
        open(2,file=f_name10(1:2)//file_name//'.out')
          do j=1,Ny_Circ
            write(2,111)(Intp_Fy_Circ(i,j),i=1,Nx_Circ)
          enddo
        close(2)
        endif


        if(f_name7.ne.' ')then
        open(2,file=f_name7(1:2)//file_name//'.out')
          do j=1,Ny_Circ
            write(2,111)(Intp_Fx_Circ(i,j),
     &                       i=1,Nx_Circ)
          enddo
        close(2)
        endif

        if(f_name8.ne.' ')then
        open(2,file=f_name8(1:2)//file_name//'.out')
          do j=1,Ny_Circ
            write(2,111)(Intp_Fy_Circ(i,j),
     &                        i=1,Nx_Circ)
          enddo
        close(2)
        endif

c more output for soulsby's formula
       goto 9197
      open(2,file='wx'//file_name//'.out')
                do j=1,Ny_Wave
            write(2,111)(pass_massfluxU(i,j),
     &                     i=1,Nx_Wave)
                enddo
       close(2)

      open(2,file='wy'//file_name//'.out')
         do j=1,Ny_Wave
           write(2,111)(pass_massfluxV(i,j),
     &                        i=1,Nx_Wave)
         enddo
         close(2)
       open(2,file='ub'//file_name//'.out')
          do j=1,Ny_Wave
              write(2,111)(pass_Ub(i,j),
     &                        i=1,Nx_Wave)
          enddo
        close(2)

        open(2,file='vb'//file_name//'.out')
           do j=1,Ny_Wave
              write(2,111)(pass_Vb(i,j),
     &                        i=1,Nx_Wave)
           enddo
        close(2)

      open(2,file='ag'//file_name//'.out')
           do j=1,Ny_Wave
            write(2,111)(pass_theta(i,j),
     &                        i=1,Nx_Wave)
           enddo
       close(2)

       open(2,file='dp'//file_name//'.out')
            do j=1,Ny_Wave
               write(2,111)(depth_circ(i,j),
     &                        i=1,Nx_Wave)
            enddo
       close(2)
       
       open(2,file='advx'//file_name//'.out')
       do j=1,Ny_Circ
       write(2,129)(advx(i,j),
     &                       i=1,Nx_Circ)
       enddo
       close(2)

      open(2,file='advy'//file_name//'.out')
      do j=1,Ny_Circ
      write(2,129)(advy(i,j),
     &                       i=1,Nx_Circ)
      enddo
      close(2)

       open(2,file='prex'//file_name//'.out')
       do j=1,Ny_Circ
       write(2,129)(PREX(i,j),
     &                       i=1,Nx_Circ)
       enddo
       close(2)
       
       open(2,file='prey'//file_name//'.out')
       do j=1,Ny_Circ
       write(2,129)(PREY(i,j),
     &                       i=1,Nx_Circ)
       enddo
       close(2)

      open(2,file='radx'//file_name//'.out')
      do j=1,Ny_Circ
      write(2,129)(RADX(i,j),
     &                       i=1,Nx_Circ)
      enddo
      close(2)

      open(2,file='rady'//file_name//'.out')
      do j=1,Ny_Circ
      write(2,129)(RADY(i,j),
     &                       i=1,Nx_Circ)
      enddo
      close(2)

      open(2,file='frcx'//file_name//'.out')
      do j=1,Ny_Circ
      write(2,129)(FRCX(i,j),
     &                       i=1,Nx_Circ)
      enddo
      close(2)

      open(2,file='frcy'//file_name//'.out')
      do j=1,Ny_Circ
      write(2,129)(FRCY(i,j),
     &                       i=1,Nx_Circ)
      enddo
      close(2)

9197  continue       

111     format(500f16.6)
129     format(500E16.6)


        return 
       end 
 
 
c ---------------------------------------------------------
        subroutine WaveModule_sample()
        USE PASS
        implicit none
        integer i,j

        if(Master_Start.eq.1)then

        print*,'wave module initialization ...'

c          write(*,*)'Do you want to run refdifs? Yes=1'
c         read(*,*)ikey
c         if(ikey.eq.1)then 
c            call refdifs
c         else
c            call load_wave
c          endif  

        else

        print*, 'call wave module ...'
c         call refdifs

        endif

        return
        end

c ---------------------------------------------------------
        subroutine CircModule_sample()
        USE PASS
        implicit none
        integer i,j

        if(Master_Start.eq.1)then

        print*,'circulation module initialization ...'

        else

        print*, 'call circulation module ...'

        endif

        return
        end

c --------------------------------------------------------
        subroutine MasterInit
        USE PASS
        implicit none
        integer i,j
 
        do j=1,Ny_Wave
        do i=1,Nx_Wave
          Pass_Sxx(i,j)=0.
          Pass_Sxy(i,j)=0.
          Pass_Syy(i,j)=0.

          Pass_Sxx_body(i,j)=0.
          Pass_Sxy_body(i,j)=0.
          Pass_Syy_body(i,j)=0.

          Pass_Sxx_surf(i,j)=0.
          Pass_Sxy_surf(i,j)=0.
          Pass_Syy_surf(i,j)=0.

          Pass_Wave_Fx(i,j)=0.
          Pass_Wave_Fy(i,j)=0.

          Pass_MassFluxU(i,j)=0.
          Pass_MassFluxV(i,j)=0.
          Pass_MassFlux(i,j)=0.

          Pass_Diss(i,j)=0.
          Pass_WaveNum(i,j)=0.
          Pass_Theta(i,j)=0.
          Pass_ubott(i,j)=0.
          Pass_Height(i,j)=0.
          Pass_C(i,j)=0.
          Pass_Cg(i,j)=0.
          Intp_U_Wave(i,j)=0.
          Intp_V_Wave(i,j)=0.
          Intp_eta_Wave(i,j)=0.
          Pass_ibrk(i,j)=0
         enddo
         enddo

         do j=1,Ny_Circ
         do i=1,Nx_Circ
          Pass_U(i,j)=0.
          Pass_V(i,j)=0.
          Pass_Ub(i,j)=0.
          Pass_Vb(i,j)=0.
          Pass_eta(i,j)=0.

          Pass_d11(i,j)=0.
          Pass_d12(i,j)=0.
          Pass_e11(i,j)=0.
          Pass_e12(i,j)=0.
          Pass_f11(i,j)=0.
          Pass_f12(i,j)=0.

          Pass_fw(i,j)=0.
          pass_vt(i,j)=0

          Intp_Fx_Circ(i,j)=0.
          Intp_Fy_Circ(i,j)=0.
          Intp_ubott_Circ(i,j)=0.
          Intp_Theta_Circ(i,j)=0.
          Intp_Sxx_Circ(i,j)=0.
          Intp_Sxy_Circ(i,j)=0.
          Intp_Syy_Circ(i,j)=0.
          Intp_Sxx_Surf(i,j)=0.
          Intp_Sxy_Surf(i,j)=0.
          Intp_Syy_Surf(i,j)=0.
          Intp_Sxx_Body(i,j)=0.
          Intp_Sxy_Body(i,j)=0.
          Intp_Syy_Body(i,j)=0.
          Intp_MassFluxU_Circ(i,j)=0.
          Intp_MassFluxV_Circ(i,j)=0.
          Intp_Diss_Circ(i,j)=0.
          Intp_ibrk_Circ(i,j)=0.

          enddo
          enddo

          do j=1,Ny_Sedi
          do i=1,Nx_Sedi

          Pass_Dupdated(i,j)=Depth_Sedi(i,j)
          Intp_U_Sedi(i,j)=0.
          Intp_V_Sedi(i,j)=0.
          Intp_Ub_Sedi(i,j)=0.
          Intp_Vb_Sedi(i,j)=0.
          Intp_ubott_Sedi(i,j)=0.
          Intp_eta_Sedi(i,j)=0.
          Intp_fw_Sedi(i,j)=0.
          Intp_vt_Sedi(i,j)=0.
          Intp_Theta_Sedi(i,j)=0.
          Intp_Height_Sedi(i,j)=0.
          Intp_ibrk_Sedi(i,j)=0.

        enddo
        enddo
   
         Pass_period = 1.
         Time_Wave_Counter = WaveSpectral_Data_Interval

        return

        end


C ------------------------------------------------------------------
c   Common blocks for interface of Nearshore Community Model
c   It is used in master program, wave module, circulation module,
c   and sediment module.
c      Fyshi 01/21/2002   
c ------------------------------------------------------------------
C
       MODULE PASS
       USE SIZES

       integer Nx_Circ, Ny_Circ, Nx_Wave,Ny_Wave,
     &         Nx_Mast, Ny_Mast, Nx_Sedi, Ny_Sedi

c -- for unstructure grid: define element numbers, point numbers keep same but Nx=1
c     Nele_Circ: element number for circulation module
c     Nele#_Circ: #th (1,2,3) point node number in element(i)

       integer Nele_Circ,Nele_Wave,Nele_Sedi,Nele_Mast

       integer, ALLOCATABLE :: Nele1_Circ(:),Nele2_Circ(:),Nele3_Circ(:)
       integer, ALLOCATABLE :: Nele1_Wave(:),Nele2_Wave(:),Nele3_Wave(:)
       integer, ALLOCATABLE :: Nele1_Sedi(:),Nele2_Sedi(:),Nele3_Sedi(:)
       integer, ALLOCATABLE :: Nele1_Mast(:),Nele2_Mast(:),Nele3_Mast(:)
   
c -- wave module:

       real(SZ),ALLOCATABLE :: Pass_Sxx(:,:), Pass_Sxy(:,:),
     &      Pass_Syy(:,:),
     &      Pass_Sxx_body(:,:),Pass_Sxy_body(:,:),
     &      Pass_Syy_body(:,:),
     &      Pass_Sxx_surf(:,:),
     &      Pass_Sxy_surf(:,:),
     &      Pass_Syy_surf(:,:),
     &      Pass_Wave_Fx(:,:),Pass_Wave_Fy(:,:), 
     &      Pass_MassFluxU(:,:),
     &      Pass_MassFluxV(:,:),
     &      Pass_MassFlux(:,:),
     &      Pass_Diss(:,:),
     &      Pass_WaveNum(:,:), Pass_Theta(:,:),
     &      Pass_ubott(:,:), Pass_Height(:,:),
     &      Pass_Cg(:,:),
     &      Pass_C(:,:),   
     &      Pass_PPer(:,:), 
     &      Intp_U_Wave(:,:), Intp_V_Wave(:,:),
     &      Intp_eta_Wave(:,:),Pass_uw(:,:,:)
 
       real (SZ) Pass_period

       real(SZ), ALLOCATABLE :: Pass_ibrk(:,:)

c -- circulation module:

       real(SZ), ALLOCATABLE :: Pass_U(:,:),Pass_V(:,:),
     &      Pass_Ub(:,:),Pass_Vb(:,:),
     &      Pass_eta(:,:),
     &      Pass_d11(:,:), Pass_d12(:,:),
     &      Pass_e11(:,:), Pass_e12(:,:),
     &      Pass_f11(:,:), Pass_f12(:,:),
     &      Pass_fw(:,:),Pass_vt(:,:),
     &      Intp_Fx_Circ(:,:),Intp_Fy_Circ(:,:),
     &      Intp_ubott_Circ(:,:),
     &      Intp_Theta_Circ(:,:),
     &      Intp_Height_Circ(:,:),
     &      Intp_C_Circ(:,:),
     &      Intp_Cg_Circ(:,:),
     &      Intp_WaveNum_Circ(:,:),
     &      Intp_Diss_Circ(:,:),
     &      Intp_ibrk_Circ(:,:),
     &      Intp_Sxx_Circ(:,:),Intp_Sxy_Circ(:,:),
     &      Intp_Syy_Circ(:,:),
     &      Intp_Sxx_Surf(:,:),Intp_Sxy_Surf(:,:),
     &      Intp_Syy_Surf(:,:),
     &      Intp_Sxx_Body(:,:),Intp_Sxy_Body(:,:),
     &      Intp_Syy_Body(:,:),
     &      Intp_MassFluxU_Circ(:,:),
     &      Intp_MassFluxV_Circ(:,:),
c -- analysis
     &      ADVX(:,:),PREX(:,:),FRCX(:,:),RADX(:,:),
     &      ADVY(:,:),PREY(:,:),FRCY(:,:),RADY(:,:)

c -- sediment module:

      real(SZ),ALLOCATABLE :: Pass_Dupdated(:,:),
     &      Pass_sedflux_x(:,:),
     &      Pass_sedflux_y(:,:),
     &      Pass_sedfluxcum_x(:,:),
     &      Pass_sedfluxcum_y(:,:),
c for ideal case fyshi
     &      Pass_sedqxum(:,:),
     &      Pass_sedqyum(:,:),
     &      Pass_sedqxslp(:,:),
     &      Pass_sedqyslp(:,:),
     &      Pass_sedqxsk(:,:),
     &      Pass_sedqysk(:,:),
c end
     &      Intp_sedfluxcum_x(:,:),
     &      Intp_sedfluxcum_y(:,:),
     &      Intp_U_Sedi(:,:),
     &      Intp_V_Sedi(:,:),
     &      Intp_Ub_Sedi(:,:),
     &      Intp_Vb_Sedi(:,:),
     &      Intp_ubott_Sedi(:,:),
     &      Intp_eta_Sedi(:,:),
     &      Intp_fw_Sedi(:,:),
     &      Intp_vt_Sedi(:,:),
     &      Intp_Theta_Sedi(:,:),
     &      Intp_Height_Sedi(:,:),
     &      Intp_ibrk_Sedi(:,:)

c -- coordinate systems, depth and ...

       real(SZ), ALLOCATABLE :: Depth_Circ(:,:),Depth_Wave(:,:),
     &      X_Wave(:,:),Y_Wave(:,:), 
     &      X_Circ(:,:),
     &      Y_Circ(:,:), 
     &      U_wind_Mast(:,:),V_wind_Mast(:,:),
     &      U_wind_Circ(:,:),V_wind_CIrc(:,:),
     &      U_wind_Wave(:,:),V_wind_Wave(:,:),
     &      Depth_Sedi(:,:),X_Sedi(:,:),
     &      Y_Sedi(:,:),
     &      X_Mast(:,:),Y_Mast(:,:),
     &      Depth_Mast(:,:)

       real(SZ) Pass_tide

c -- vector rotate
       real(SZ) Circ_Rotate_Angle, Wave_Rotate_Angle, 
     &      Sedi_Rotate_Angle       
 
c -- control parameters

       integer Master_Start,nWave,nCirc,nSedi,nOut,istep

       integer Wave_Stag_huv(3), Circ_Stag_huv(3),Sedi_Stag_huv(3)

       real(SZ) N_Interval_CallWave,
     &     N_Interval_CallCirc,N_Interval_CallSedi,
     &     N_Delay_CallSedi,N_Interval_Output,
     &     WaveSpectral_Data_Interval   
       real(SZ) Total_Time, Master_dt,Time_Master,Time_Wave_Counter
       logical Grid_Mast_Wave_Same, Grid_Mast_Circ_Same,
     &         Grid_Mast_Sedi_Same, Grid_Wave_Circ_Same,
     &         Grid_Wave_Sedi_Same, Grid_Circ_Sedi_Same,
     &         Wave_Staggered, Circ_Staggered,Sedi_Staggered,
     &         Wave_Structured, Circ_Structured, Sedi_Structured,
     &         Grid_Extrapolation,
     &         Wave_Curr_Interact,
     &         Wave_Bed_Interact,
     &         Curr_Bed_Interact,
     &         Wave_To_Circ_Height,
     &         Wave_To_Circ_Angle,
     &         Wave_To_Circ_WaveNum,
     &         Wave_To_Circ_C,
     &         Wave_To_Circ_Cg,
     &         Wave_To_Circ_Radiation,
     &         Wave_To_Circ_Rad_Surf,
     &         Wave_To_Circ_Rad_Body,
     &         Wave_To_Circ_Forcing,
     &         Wave_To_Circ_MassFlux,
     &         Wave_To_Circ_Dissipation,
     &         Wave_To_Circ_BottomUV,
     &         Wave_To_Circ_Brkindex,
     &         Circ_To_Wave_UV,
     &         Circ_To_Wave_eta,
     &         Wave_To_Sedi_Height,
     &         Wave_To_Sedi_Angle,
     &         Wave_To_Sedi_BottomUV,
     &         Circ_To_Sedi_SedFlux,
     &         Circ_To_Sedi_UV,
     &         Circ_To_Sedi_UVb,
     &         Circ_To_Sedi_eta,
     &         Circ_To_Sedi_UV3D,
     &         Circ_To_Sedi_fw,
     &         Circ_To_Sedi_vt,
     &         Circ_To_Sedi_UVquasi3D,
     &         Sedi_To_Wave_Depth,
     &         Sedi_To_Circ_Depth,
     &         Circ_POM,
     &         Circ_SC,
     &         Waveupdat_for_Circ,
     &         Waveupdat_for_Sedi,
     &         Circupdat_for_Wave,
     &         Circupdat_for_Sedi,
     &         Sediupdat_for_Wave,
     &         Sediupdat_for_Circ

c -- file names

      character*255 f_depth,f_xymast,f_xywave,f_xycirc,f_xysedi,
     &              f_name6,f_name7,
     &              f_name8,
     &              f_name9,f_name10,f_name11,f_name12,
     &              f_name13,f_name14,
     &              f_name15,f_name16


c -- allocate pass variables
       CONTAINS

       subroutine ALLOCATE_ELEMENT

       ALLOCATE (Nele1_Circ(Nele_Circ),Nele2_Circ(Nele_Circ),
     &           Nele3_Circ(Nele_Circ) )
       ALLOCATE (Nele1_Wave(Nele_Wave),Nele2_Wave(Nele_Wave),
     &           Nele3_Wave(Nele_Wave) )
       ALLOCATE (Nele1_Sedi(Nele_Sedi),Nele2_Sedi(Nele_Sedi),
     &           Nele3_Sedi(Nele_Sedi) )
       ALLOCATE (Nele1_Mast(Nele_Mast),Nele2_Mast(Nele_Mast),
     &           Nele3_Mast(Nele_Mast) )


       return
       end subroutine ALLOCATE_ELEMENT

       subroutine ALLOCATE_PASS


c -- wave module
       ALLOCATE ( Pass_Sxx(Nx_Wave,Ny_Wave), 
     &      Pass_Sxy(Nx_Wave,Ny_Wave),
     &      Pass_Syy(Nx_Wave,Ny_Wave),
     &      Pass_Sxx_body(Nx_Wave,Ny_Wave),
     &      Pass_Sxy_body(Nx_Wave,Ny_Wave),
     &      Pass_Syy_body(Nx_Wave,Ny_Wave),
     &      Pass_Sxx_surf(Nx_Wave,Ny_Wave),
     &      Pass_Sxy_surf(Nx_Wave,Ny_Wave),
     &      Pass_Syy_surf(Nx_Wave,Ny_Wave),
     &      Pass_Wave_Fx(Nx_Wave,Ny_Wave),
     &      Pass_Wave_Fy(Nx_Wave,Ny_Wave),
     &      Pass_MassFluxU(Nx_Wave,Ny_Wave),
     &      Pass_MassFluxV(Nx_Wave,Ny_Wave),
     &      Pass_MassFlux(Nx_Wave,Ny_Wave),
     &      Pass_Diss(Nx_Wave,Ny_Wave),
     &      Pass_WaveNum(Nx_Wave,Ny_Wave), 
     &      Pass_Theta(Nx_Wave,Ny_Wave),
     &      Pass_ubott(Nx_Wave,Ny_Wave), 
     &      Pass_Height(Nx_Wave,Ny_Wave),
     &      Pass_Cg(Nx_Wave,Ny_Wave),
     &      Pass_C(Nx_Wave,Ny_Wave),
     &      Pass_PPER(Nx_Wave,Ny_Wave),
     &      Intp_U_Wave(Nx_Wave,Ny_Wave), Intp_V_Wave(Nx_Wave,Ny_Wave),
     &      Intp_eta_Wave(Nx_Wave,Ny_Wave),
     &      Pass_uw(nx_wave,ny_wave,100),
     &      Pass_ibrk(Nx_Wave,Ny_Wave) )

c -- circulation module:

       ALLOCATE (  Pass_U(Nx_Circ,Ny_Circ),Pass_V(Nx_Circ,Ny_Circ),
     &      Pass_Ub(Nx_Circ,Ny_Circ),Pass_Vb(Nx_Circ,Ny_Circ),
     &      Pass_eta(Nx_Circ,Ny_Circ),
     &      Pass_d11(Nx_Circ,Ny_Circ), Pass_d12(Nx_Circ,Ny_Circ),
     &      Pass_e11(Nx_Circ,Ny_Circ), Pass_e12(Nx_Circ,Ny_Circ),
     &      Pass_f11(Nx_Circ,Ny_Circ), Pass_f12(Nx_Circ,Ny_Circ),
     &      Pass_fw(Nx_Circ,Ny_Circ),Pass_vt(Nx_Circ,Ny_Circ),
     &      Intp_Fx_Circ(Nx_Circ,Ny_Circ),Intp_Fy_Circ(Nx_Circ,Ny_Circ),
     &      Intp_ubott_Circ(Nx_Circ,Ny_Circ),
     &      Intp_Theta_Circ(Nx_Circ,Ny_Circ),
     &      Intp_Height_Circ(Nx_Circ,Ny_Circ),
     &      Intp_C_Circ(Nx_Circ,Ny_Circ),
     &      Intp_Cg_Circ(Nx_Circ,Ny_Circ),
     &      Intp_WaveNum_Circ(Nx_Circ,Ny_Circ),
     &      Intp_Diss_Circ(Nx_Circ,Ny_Circ),
     &      Intp_ibrk_Circ(Nx_Circ,Ny_Circ),
     &      Intp_Sxx_Circ(Nx_Circ,Ny_Circ),
     &      Intp_Sxy_Circ(Nx_Circ,Ny_Circ),
     &      Intp_Syy_Circ(Nx_Circ,Ny_Circ),
     &      Intp_Sxx_Surf(Nx_Circ,Ny_Circ),
     &      Intp_Sxy_Surf(Nx_Circ,Ny_Circ),
     &      Intp_Syy_Surf(Nx_Circ,Ny_Circ),
     &      Intp_Sxx_Body(Nx_Circ,Ny_Circ),
     &      Intp_Sxy_Body(Nx_Circ,Ny_Circ),
     &      Intp_Syy_Body(Nx_Circ,Ny_Circ),
     &      Intp_MassFluxU_Circ(Nx_Circ,Ny_Circ),
     &      Intp_MassFluxV_Circ(Nx_Circ,Ny_Circ),
c -- analysis
     &      ADVX(Nx_Circ,Ny_Circ),PREX(Nx_Circ,Ny_Circ),
     &      FRCX(Nx_Circ,Ny_Circ),RADX(Nx_Circ,Ny_Circ),     
     &      ADVY(Nx_Circ,Ny_Circ),PREY(Nx_Circ,Ny_Circ),
     &      FRCY(Nx_Circ,Ny_Circ),RADY(Nx_Circ,Ny_Circ))

c -- sediment module:

      ALLOCATE ( Pass_Dupdated(Nx_Sedi,Ny_Sedi),
     &      Pass_sedflux_x(Nx_Sedi,Ny_Sedi),
     &      Pass_sedflux_y(Nx_Sedi,Ny_Sedi),
     &      Pass_sedfluxcum_x(Nx_Sedi,Ny_Sedi),
     &      Pass_sedfluxcum_y(Nx_Sedi,Ny_Sedi),
c for ideal case fyshi
     &      Pass_sedqxum(Nx_Sedi,Ny_Sedi),
     &      Pass_sedqyum(Nx_Sedi,Ny_Sedi),
     &      Pass_sedqxslp(Nx_Sedi,Ny_Sedi),
     &      Pass_sedqyslp(Nx_Sedi,Ny_Sedi),
     &      Pass_sedqxsk(Nx_Sedi,Ny_Sedi),
     &      Pass_sedqysk(Nx_Sedi,Ny_Sedi),
c end
     &      Intp_sedfluxcum_x(Nx_Sedi,Ny_Sedi),
     &      Intp_sedfluxcum_y(Nx_Sedi,Ny_Sedi),
     &      Intp_U_Sedi(Nx_Sedi,Ny_Sedi),
     &      Intp_V_Sedi(Nx_Sedi,Ny_Sedi),
     &      Intp_Ub_Sedi(Nx_Sedi,Ny_Sedi),
     &      Intp_Vb_Sedi(Nx_Sedi,Ny_Sedi),
     &      Intp_ubott_Sedi(Nx_Sedi,Ny_Sedi),
     &      Intp_eta_Sedi(Nx_Sedi,Ny_Sedi),
     &      Intp_fw_Sedi(Nx_Sedi,Ny_Sedi),
     &      Intp_vt_Sedi(Nx_Sedi,Ny_Sedi),
     &      Intp_Theta_Sedi(Nx_Sedi,Ny_Sedi),
     &      Intp_Height_Sedi(Nx_Sedi,Ny_Sedi),
     &      Intp_ibrk_Sedi(Nx_Sedi,Ny_Sedi) )

c -- coordinate systems, depth and ...

       allocate ( Depth_Circ(Nx_Circ,Ny_Circ),
     &      Depth_Wave(Nx_Wave,Ny_Wave),
     &      X_Wave(Nx_Wave,Ny_Wave),Y_Wave(Nx_Wave,Ny_Wave),
     &      X_Circ(Nx_Circ,Ny_Circ),
     &      Y_Circ(Nx_Circ,Ny_Circ), 
     &      U_wind_Mast(Nx_Mast,Ny_Mast),V_wind_Mast(Nx_Mast,Ny_Mast),
     &      U_wind_Circ(Nx_Circ,Ny_Circ),V_wind_Circ(Nx_Circ,Ny_Circ),
     &      U_wind_Wave(Nx_Wave,Ny_Wave),V_wind_Wave(Nx_Wave,Ny_Wave),
     &      Depth_Sedi(Nx_Sedi,Ny_Sedi),X_Sedi(Nx_Sedi,Ny_Sedi),
     &      Y_Sedi(Nx_Sedi,Ny_Sedi),
     &      X_Mast(Nx_Mast,Ny_Mast),Y_Mast(Nx_Mast,Ny_Mast),
     &      Depth_Mast(Nx_Mast,Ny_Mast) )


        return
        END SUBROUTINE allocate_pass

       END MODULE PASS

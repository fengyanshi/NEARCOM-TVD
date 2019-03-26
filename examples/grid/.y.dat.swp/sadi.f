	Module SADI

        USE PASS
        USE Interp_coef
        implicit none
!        include 'pass.h'
!        include 'interp.h'
	integer i,j
        integer mtrim,ntrim,mtrigm,ktrim,ntide,nsta_max
        real g,rho,pi,dxi,deta
        parameter(mtrim=nx_max,ntrim=ny_max,g=9.80,ktrim=20)
        parameter(ntide=10)
        parameter(nsta_max=100)
        parameter(mtrigm=max(nx_max,ny_max))
        parameter(rho=1030.,pi=3.1415926)
        parameter(dxi=1.,deta=1.)
	integer m,n,itstep,ntime,kplot,ik,k_3d,kstation
        real uamplify,plotintra,plotstart,sta_interval
        real time_sta,time_plot
        integer i_sta(nsta_max),j_sta(nsta_max),n_station
        integer kdisplay
        real delta
        real rmanning
        real hi
        real f
        real dt
        integer ibeg(ntrim),iend(ntrim),jbeg(mtrim),jend(mtrim)
        integer ileft,iright,mtrig,im,i_tr,jm,j_tr,jbot,jtop
        integer btype_left,btype_right,btype_bot,btype_top
        real u_rj,u_H,u_v,u_beta,u_alpha,u_gamma,u_u
        real v_rj,v_H,v_u,v_beta,v_alpha,v_gamma,v_v
        real H_1,H_2,H_3,H_4,H_,HJ
        real UUH_xi,UVH_eta,Z_ctr,Z_t,Z_b,Z_eta,Eadv1,Eadv2,Eadv3,
     &       Ecor1,Ecor2,Egrd1,Egrd2,Ewind,Cd,Efrc,Ediff,Ewave,
     &       Efrc_corr,Edisp,
     &       H_L,H_R,
     &       Z_,z_rj,Res,VVH_eta,UVH_xi,Z_r,Z_l,Z_xi,H_T,H_B
        real UU(mtrim,ntrim),VV(mtrim,ntrim)
        real U(mtrim,ntrim),V(mtrim,ntrim)
        real Z(mtrim,ntrim)
        real H(mtrim,ntrim),HH(mtrim,ntrim)
        real X(mtrim,ntrim),Y(mtrim,ntrim)
        real Taux(mtrim,ntrim),Tauy(mtrim,ntrim)
        real rJ(mtrim,ntrim),Alpha(mtrim,ntrim),Gamma(mtrim,ntrim)
        real Beta(mtrim,ntrim)
        real Adu1(mtrim,ntrim),Adu2(mtrim,ntrim),Adu3(mtrim,ntrim)
        real Adv1(mtrim,ntrim),Adv2(mtrim,ntrim),Adv3(mtrim,ntrim)
        real Xxi(mtrim,ntrim),Yxi(mtrim,ntrim),Xeta(mtrim,ntrim)
        real Yeta(mtrim,ntrim)
        integer Ntype(mtrim,ntrim),Ntype1(mtrim,ntrim)
        integer Mask(mtrim,ntrim)
        real Tao_ax,Tao_ay
        real Zopl(ntrim),Zopr(ntrim),Zopt(mtrim),Zopb(mtrim)
        real Uopl(ntrim),Uopr(ntrim),Vopt(mtrim),Vopb(mtrim)
        real uconv(mtrim,ntrim),vconv(mtrim,ntrim)
        real zconv(mtrim,ntrim)

! boundary condition
        integer east_ele2_flx3_grd4_rad5_per6,jstart_E,jend_E,
     &          west_ele2_flx3_grd4_rad5_per6,jstart_W,jend_W,
     &          south_ele2_flx3_grd4_rad5_per6,istart_S,iend_S,
     &          north_ele2_flx3_grd4_rad5_per6,istart_N,iend_N,
     &          east_data1_anly2,west_data1_anly2,
     &          south_data1_anly2,north_data1_anly2

! shorecirc module
	integer disp3d,curv_grid,orth_grid

        real c_subgrid,CR,depth_min,f_cwc,
     *                       anu_tb_const,anu_t0_const

	real trd1(mtrim,ntrim),tr1(mtrim,ntrim)

	real dx,dy,time_circ

	real sq_g_0(mtrim,ntrim),g_11(mtrim,ntrim),g_12(mtrim,ntrim)
     *   ,g_22(mtrim,ntrim)
     *   ,sq_g_0_u(mtrim,ntrim),g_11_u(mtrim,ntrim),g_12_u(mtrim,ntrim)
     *   ,g_22_u(mtrim,ntrim),sq_g_0_v(mtrim,ntrim),g_11_v(mtrim,ntrim)
     *   ,g_12_v(mtrim,ntrim),g_22_v(mtrim,ntrim),sq_g_0_z(mtrim,ntrim)
     *   ,g_0_z(mtrim,ntrim),g_12_z(mtrim,ntrim),g_22_z(mtrim,ntrim)
     *   ,g_11_z(mtrim,ntrim),g_0(mtrim,ntrim)
     
        real Cf111(mtrim,ntrim),Cf112(mtrim,ntrim),Cf122(mtrim,ntrim)
     *              ,Cf211(mtrim,ntrim),Cf212(mtrim,ntrim)
     *              ,Cf222(mtrim,ntrim)

        real X_xi_1(mtrim,ntrim), X_xi_2(mtrim,ntrim)
     *               ,Y_xi_1(mtrim,ntrim), Y_xi_2(mtrim,ntrim)

        real Dup11(mtrim,ntrim),Dup12(mtrim,ntrim),Dup22(mtrim,ntrim)
     *               ,Ddn11(mtrim,ntrim),Ddn12(mtrim,ntrim)
     *               ,Ddn22(mtrim,ntrim),anu_s(mtrim,ntrim)
     *               ,anu(mtrim,ntrim),anu_t(mtrim,ntrim)
     *               ,anu_t0(mtrim,ntrim),anu_tb(mtrim,ntrim)
     *               ,anu_z(mtrim,ntrim)
	
	real Q1(mtrim,ntrim),Q2(mtrim,ntrim),zeta(mtrim,ntrim)
     *               ,Q1_v(mtrim,ntrim),Q1_z(mtrim,ntrim)
     *               ,Q1_c(mtrim,ntrim),Q2_u(mtrim,ntrim)
     *               ,Q2_z(mtrim,ntrim),Q2_c(mtrim,ntrim)
     *               ,zeta_u(mtrim,ntrim),zeta_v(mtrim,ntrim)
     *               ,zeta_c(mtrim,ntrim),depth(mtrim,ntrim)
     *               ,depth_c(mtrim,ntrim),depth_u(mtrim,ntrim)
     *               ,depth_v(mtrim,ntrim),zeta_0(mtrim,ntrim)
     *               ,Q1_0(mtrim,ntrim),Q2_0(mtrim,ntrim)
     *               ,Qx(mtrim,ntrim),Qy(mtrim,ntrim)
     *               ,Q1_p(mtrim,ntrim),Q2_p(mtrim,ntrim)
     *               ,zeta_p(mtrim,ntrim)

	real diffx(mtrim,ntrim),diffy(mtrim,ntrim)
	real dispx(mtrim,ntrim),dispy(mtrim,ntrim)

        real wavefx(mtrim,ntrim),wavefy(mtrim,ntrim)

        real Tmp1(mtrim,ntrim),Tmp2(mtrim,ntrim),Tmp3(mtrim,ntrim)	

	integer bc_specify,specify_west_east,specify_south_north
	integer bc_periodic,periodic_west_east,periodic_south_north

! --- for bottom friction

	real Tau_1(mtrim,ntrim), Tau_2(mtrim,ntrim)
     *                ,Tau_x(mtrim,ntrim), Tau_y(mtrim,ntrim)
     *                ,beta1(mtrim,ntrim),beta2(mtrim,ntrim)
     *                ,fric_x(mtrim,ntrim), fric_y(mtrim,ntrim)
     *                ,f_cw(mtrim,ntrim),U_wind(mtrim,ntrim)
     *                ,V_wind(mtrim,ntrim),Ta_x(mtrim,ntrim)
     *                ,Ta_y(mtrim,ntrim)

! --- wind
         real wind_u(mtrim,ntrim),wind_v(mtrim,ntrim)
		 real windspd(2),winddir(2),windtime(2)
		 real wind_intv

! --- tide
      real amp_west(ntrim,ntide),pha_west(ntrim,ntide)
      real amp_east(ntrim,ntide),pha_east(ntrim,ntide)
      real amp_north(ntrim,ntide),pha_north(ntrim,ntide)
      real amp_south(ntrim,ntide),pha_south(ntrim,ntide)

! --- wave and radiation stresses

	real S_xx(mtrim,ntrim),S_xy(mtrim,ntrim),S_yy(mtrim,ntrim)
     *         ,S_11(mtrim,ntrim),S_12(mtrim,ntrim),S_22(mtrim,ntrim)
     *         , Height(mtrim,ntrim),angle(mtrim,ntrim)
     *         ,C_sw(mtrim,ntrim),Qw(mtrim,ntrim),Qwx(mtrim,ntrim)
     *         ,Qwy(mtrim,ntrim),Qw1(mtrim,ntrim),Qw2(mtrim,ntrim)
     *         ,u0(mtrim,ntrim),Cg(mtrim,ntrim), Tauss_x(mtrim,ntrim)
     *         ,Tauss_y(mtrim,ntrim),Tauss_1(mtrim,ntrim)
     *         ,Tauss_2(mtrim,ntrim)
     *         ,shwave_force_x(mtrim,ntrim),shwave_force_y(mtrim,ntrim)
     *         ,shwave_force_1(mtrim,ntrim)
     *         ,shwave_force_2(mtrim,ntrim),dissipation(mtrim,ntrim)

	integer ibrk(mtrim,ntrim)


  	real period,B0

! --- 3-D dispersion terms

	real F_1(mtrim,ntrim),F_2(mtrim,ntrim),dd_11(mtrim,ntrim)
     *         ,dd_12(mtrim,ntrim)
     *         ,ee_11(mtrim,ntrim),ee_12(mtrim,ntrim)
     *         ,ff_11(mtrim,ntrim),ff_12(mtrim,ntrim)
     *         ,dd_11_cart(mtrim,ntrim),dd_12_cart(mtrim,ntrim)
     *         ,ee_11_cart(mtrim,ntrim),ee_12_cart(mtrim,ntrim)
     *         ,ff_11_cart(mtrim,ntrim),ff_12_cart(mtrim,ntrim)
     *         , A_111(mtrim,ntrim),A_112(mtrim,ntrim)
     *         ,A_121(mtrim,ntrim),A_122(mtrim,ntrim)
     *         ,A_211(mtrim,ntrim),A_212(mtrim,ntrim)
     *         ,A_221(mtrim,ntrim),A_222(mtrim,ntrim)
     *         , B_11(mtrim,ntrim),B_12(mtrim,ntrim)
     *         ,B_21(mtrim,ntrim),B_22(mtrim,ntrim)
     *         ,D_11(mtrim,ntrim),D_12(mtrim,ntrim),D_21(mtrim,ntrim)
     *         ,D_22(mtrim,ntrim),M_11(mtrim,ntrim),M_12(mtrim,ntrim)
     *         ,M_21(mtrim,ntrim),M_22(mtrim,ntrim)
     *         ,D_11_z(mtrim,ntrim),D_12_z(mtrim,ntrim)
     *         ,D_22_z(mtrim,ntrim)

! - 3D velocity
        real U_3d(mtrim,ntrim,ktrim),V_3d(mtrim,ntrim,ktrim)

! common blocks

!- for semi-implicit schemes
	common m,n,itstep,k_3d
        common uamplify
        common kdisplay,plotstart,plotintra,sta_interval
        common i_sta,j_sta,n_station,kstation  
        common delta
        common rmanning
        common hi
        common f
        common dt,ntime
        common ibeg,iend,jbeg,jend
        common UU,VV
        common U,V
        common Z
        common H,HH
        common X,Y
        common Taux,Tauy
        common rJ,Alpha,Gamma
        common Beta
        common Adu1,Adu2,Adu3
        common Adv1,Adv2,Adv3
        common Xxi,Yxi,Xeta
        common Yeta
        common Ntype,Ntype1,Mask
        common Tao_ax,Tao_ay,Cd
        common Zopl,Zopr,Zopt,Zopb
        common Uopl,Uopr,Vopt,Vopb
        common uconv,vconv,zconv

! shorecirc module
	common disp3d,curv_grid,orth_grid
        common c_subgrid,CR,depth_min,f_cwc,
     *                       anu_tb_const,anu_t0_const
	common trd1,tr1
	common dx,dy,time_circ
	common sq_g_0,g_11,g_12,g_22
     *   ,sq_g_0_u,g_11_u,g_12_u,g_22_u,sq_g_0_v,g_11_v
     *   ,g_12_v,g_22_v,sq_g_0_z,g_0_z,g_12_z,g_22_z,g_11_z
     *   ,g_0
     
        common Cf111,Cf112,Cf122
     *              ,Cf211,Cf212,Cf222

        common X_xi_1, X_xi_2
     *               ,Y_xi_1, Y_xi_2

        common Dup11,Dup12,Dup22
     *               ,Ddn11,Ddn12,Ddn22,anu_s,anu,anu_t,anu_t0
     *               ,anu_tb,anu_z
	
	common Q1,Q2,zeta,Q1_v,Q1_z
     *               ,Q1_c,Q2_u,Q2_z,Q2_c,zeta_u,zeta_v,zeta_c,depth
     *               ,depth_c,depth_u,depth_v,zeta_0,Q1_0,Q2_0
     *               ,Qx,Qy,Q1_p,Q2_p,zeta_p

        common Tmp1,Tmp2,Tmp3	

	common diffx,diffy,dispx,dispy,wavefx,wavefy

	common bc_specify,specify_west_east,specify_south_north
	common bc_periodic,periodic_west_east,periodic_south_north

! --- for bottom friction

	common Tau_1, Tau_2,Tau_x,Tau_y
     *                ,beta1,beta2
     *                ,fric_x, fric_y,f_cw,U_wind,V_wind,Ta_x,Ta_y


! --- wind
         common wind_u,wind_v
		 common windtime,windspd,winddir,wind_intv
c         common Circ_wind

! --- tide
      common amp_west,pha_west
      common amp_east,pha_east
      common amp_north,pha_north
      common amp_south,pha_south

! --- wave and radiation stresses

	common S_xx,S_xy,S_yy
     *         ,S_11,S_12,S_22, Height,angle,C_sw,Qw,Qwx,Qwy,Qw1,Qw2
     *         ,u0,Cg, Tauss_x,Tauss_y,Tauss_1,Tauss_2
     *         ,shwave_force_x,shwave_force_y,shwave_force_1
     *         ,shwave_force_2,dissipation

	common ibrk


  	common period,B0


! --- 3-D dispersion terms

	common F_1,F_2,dd_11,dd_12
     *         ,ee_11,ee_12,ff_11,ff_12, A_111,A_112,A_121,A_122
     *         ,ee_11_cart,ee_12_cart,ff_11_cart,ff_12_cart
     *         ,dd_11_cart,dd_12_cart
     *         ,A_211,A_212,A_221,A_222, B_11,B_12,B_21,B_22
     *         ,D_11,D_12,D_21,D_22,M_11,M_12,M_21,M_22
     *         ,D_11_z,D_12_z,D_22_z

! - 3D velocity
        common U_3d,V_3d

! boundary condition
        common  east_ele2_flx3_grd4_rad5_per6,jstart_E,jend_E,
     &          west_ele2_flx3_grd4_rad5_per6,jstart_W,jend_W,
     &          south_ele2_flx3_grd4_rad5_per6,istart_S,iend_S,
     &          north_ele2_flx3_grd4_rad5_per6,istart_N,iend_N,
     &          east_data1_anly2,west_data1_anly2,
     &          south_data1_anly2,north_data1_anly2

       end

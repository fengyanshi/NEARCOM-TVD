C    interpolation coef
      MODULE INTERP_COEF
      USE  SIZES
      real(SZ),ALLOCATABLE ::
     &       Sc_01(:,:),S1_01(:,:),
     &       S2_01(:,:),S3_01(:,:),
     &       Sc_02(:,:),S1_02(:,:),
     &       S2_02(:,:),S3_02(:,:),
     &       Sc_03(:,:),S1_03(:,:),
     &       S2_03(:,:),S3_03(:,:),
     &       Sc_12(:,:),S1_12(:,:),
     &       S2_12(:,:),S3_12(:,:),
     &       Sc_13(:,:),S1_13(:,:),
     &       S2_13(:,:),S3_13(:,:),
     &       Sc_21(:,:),S1_21(:,:),
     &       S2_21(:,:),S3_21(:,:),
     &       Sc_23(:,:),S1_23(:,:),
     &       S2_23(:,:),S3_23(:,:),
     &       Sc_31(:,:),S1_31(:,:),
     &       S2_31(:,:),S3_31(:,:),
     &       Sc_32(:,:),S1_32(:,:),
     &       S2_32(:,:),S3_32(:,:),
c -- for adcirc
     &       Sc_20(:,:),S1_20(:,:),
     &       S2_20(:,:),S3_20(:,:)


      integer,ALLOCATABLE :: 
     &        nx1_01(:,:),ny1_01(:,:)
     &       ,nx2_01(:,:),ny2_01(:,:)
     &       ,nx3_01(:,:),ny3_01(:,:)
     &       ,nx1_02(:,:),ny1_02(:,:)
     &       ,nx2_02(:,:),ny2_02(:,:)
     &       ,nx3_02(:,:),ny3_02(:,:)
     &       ,nx1_03(:,:),ny1_03(:,:)
     &       ,nx2_03(:,:),ny2_03(:,:)
     &       ,nx3_03(:,:),ny3_03(:,:)
     &       ,nx1_12(:,:),ny1_12(:,:)
     &       ,nx2_12(:,:),ny2_12(:,:)
     &       ,nx3_12(:,:),ny3_12(:,:)
     &       ,nx1_13(:,:),ny1_13(:,:)
     &       ,nx2_13(:,:),ny2_13(:,:)
     &       ,nx3_13(:,:),ny3_13(:,:)
     &       ,nx1_21(:,:),ny1_21(:,:)
     &       ,nx2_21(:,:),ny2_21(:,:)
     &       ,nx3_21(:,:),ny3_21(:,:)
     &       ,nx1_23(:,:),ny1_23(:,:)
     &       ,nx2_23(:,:),ny2_23(:,:)
     &       ,nx3_23(:,:),ny3_23(:,:)
     &       ,nx1_31(:,:),ny1_31(:,:)
     &       ,nx2_31(:,:),ny2_31(:,:)
     &       ,nx3_31(:,:),ny3_31(:,:)
     &       ,nx1_32(:,:),ny1_32(:,:)
     &       ,nx2_32(:,:),ny2_32(:,:)
     &       ,nx3_32(:,:),ny3_32(:,:)
c -- for adcirc
     &       ,nx1_20(:,:),ny1_20(:,:)
     &       ,nx2_20(:,:),ny2_20(:,:)
     &       ,nx3_20(:,:),ny3_20(:,:)


        CONTAINS

        SUBROUTINE ALLOCATE_INTERP_COEF
        USE PASS
      ALLOCATE (   Sc_01(Nx_Wave,Ny_Wave),S1_01(Nx_Wave,Ny_Wave),
     &       S2_01(Nx_Wave,Ny_Wave),S3_01(Nx_Wave,Ny_Wave),
     &       Sc_02(Nx_Circ,Ny_Circ),S1_02(Nx_Circ,Ny_Circ),
     &       S2_02(Nx_Circ,Ny_Circ),S3_02(Nx_Circ,Ny_Circ),
     &       Sc_03(Nx_Sedi,Ny_Sedi),S1_03(Nx_Sedi,Ny_Sedi),
     &       S2_03(Nx_Sedi,Ny_Sedi),S3_03(Nx_Sedi,Ny_Sedi),
     &       Sc_12(Nx_Circ,Ny_Circ),S1_12(Nx_Circ,Ny_Circ),
     &       S2_12(Nx_Circ,Ny_Circ),S3_12(Nx_Circ,Ny_Circ),
     &       Sc_13(Nx_Sedi,Ny_Sedi),S1_13(Nx_Sedi,Ny_Sedi),
     &       S2_13(Nx_Sedi,Ny_Sedi),S3_13(Nx_Sedi,Ny_Sedi),
     &       Sc_21(Nx_Wave,Ny_Wave),S1_21(Nx_Wave,Ny_Wave),
     &       S2_21(Nx_Wave,Ny_Wave),S3_21(Nx_Wave,Ny_Wave),
     &       Sc_23(Nx_Sedi,Ny_Sedi),S1_23(Nx_Sedi,Ny_Sedi),
     &       S2_23(Nx_Sedi,Ny_Sedi),S3_23(Nx_Sedi,Ny_Sedi),
     &       Sc_31(Nx_Wave,Ny_Wave),S1_31(Nx_Wave,Ny_Wave),
     &       S2_31(Nx_Wave,Ny_Wave),S3_31(Nx_Wave,Ny_Wave),
     &       Sc_32(Nx_Circ,Ny_Circ),S1_32(Nx_Circ,Ny_Circ),
     &       S2_32(Nx_Circ,Ny_Circ),S3_32(Nx_Circ,Ny_Circ),
CC- for adcirc
     &       Sc_20(Nx_Mast,Ny_Mast),S1_20(Nx_Mast,Ny_Mast),
     &       S2_20(Nx_Mast,Ny_Mast),S3_20(Nx_Mast,Ny_Mast)
     &       )


      ALLOCATE ( nx1_01(Nx_Wave,Ny_Wave),ny1_01(Nx_Wave,Ny_Wave)
     &       ,nx2_01(Nx_Wave,Ny_Wave),ny2_01(Nx_Wave,Ny_Wave)
     &       ,nx3_01(Nx_Wave,Ny_Wave),ny3_01(Nx_Wave,Ny_Wave)
     &       ,nx1_02(Nx_Circ,Ny_Circ),ny1_02(Nx_Circ,Ny_Circ)
     &       ,nx2_02(Nx_Circ,Ny_Circ),ny2_02(Nx_Circ,Ny_Circ)
     &       ,nx3_02(Nx_Circ,Ny_Circ),ny3_02(Nx_Circ,Ny_Circ)
     &       ,nx1_03(Nx_Sedi,Ny_Sedi),ny1_03(Nx_Sedi,Ny_Sedi)
     &       ,nx2_03(Nx_Sedi,Ny_Sedi),ny2_03(Nx_Sedi,Ny_Sedi)
     &       ,nx3_03(Nx_Sedi,Ny_Sedi),ny3_03(Nx_Sedi,Ny_Sedi)
     &       ,nx1_12(Nx_Circ,Ny_Circ),ny1_12(Nx_Circ,Ny_Circ)
     &       ,nx2_12(Nx_Circ,Ny_Circ),ny2_12(Nx_Circ,Ny_Circ)
     &       ,nx3_12(Nx_Circ,Ny_Circ),ny3_12(Nx_Circ,Ny_Circ)
     &       ,nx1_13(Nx_Sedi,Ny_Sedi),ny1_13(Nx_Sedi,Ny_Sedi)
     &       ,nx2_13(Nx_Sedi,Ny_Sedi),ny2_13(Nx_Sedi,Ny_Sedi)
     &       ,nx3_13(Nx_Sedi,Ny_Sedi),ny3_13(Nx_Sedi,Ny_Sedi)
     &       ,nx1_21(Nx_Wave,Ny_Wave),ny1_21(Nx_Wave,Ny_Wave)
     &       ,nx2_21(Nx_Wave,Ny_Wave),ny2_21(Nx_Wave,Ny_Wave)
     &       ,nx3_21(Nx_Wave,Ny_Wave),ny3_21(Nx_Wave,Ny_Wave)
     &       ,nx1_23(Nx_Sedi,Ny_Sedi),ny1_23(Nx_Sedi,Ny_Sedi)
     &       ,nx2_23(Nx_Sedi,Ny_Sedi),ny2_23(Nx_Sedi,Ny_Sedi)
     &       ,nx3_23(Nx_Sedi,Ny_Sedi),ny3_23(Nx_Sedi,Ny_Sedi)
     &       ,nx1_31(Nx_Wave,Ny_Wave),ny1_31(Nx_Wave,Ny_Wave)
     &       ,nx2_31(Nx_Wave,Ny_Wave),ny2_31(Nx_Wave,Ny_Wave)
     &       ,nx3_31(Nx_Wave,Ny_Wave),ny3_31(Nx_Wave,Ny_Wave)
     &       ,nx1_32(Nx_Circ,Ny_Circ),ny1_32(Nx_Circ,Ny_Circ)
     &       ,nx2_32(Nx_Circ,Ny_Circ),ny2_32(Nx_Circ,Ny_Circ)
     &       ,nx3_32(Nx_Circ,Ny_Circ),ny3_32(Nx_Circ,Ny_Circ) 
     &       ,nx1_20(Nx_Mast,Ny_Mast),ny1_20(Nx_Mast,Ny_Mast)
     &       ,nx2_20(Nx_Mast,Ny_Mast),ny2_20(Nx_Mast,Ny_Mast)
     &       ,nx3_20(Nx_Mast,Ny_Mast),ny3_20(Nx_Mast,Ny_Mast)
     &       )




        RETURN
        END SUBROUTINE ALLOCATE_INTERP_COEF

        END MODULE INTERP_COEF

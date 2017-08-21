recent modification

12/08/2013

in sc_main.fsc

1) subroutine update mask
   added MinDepth for flooding
2) call update_mask inside GET_Eta_U_V_HU_HV

10/23/2013
added COUPLING_NO_UV

PREVIOUS MODIFICATIONS
1) add wave breaking fraction

> mod_pass.fsc
         WaveBrFraSW, WaveBrFraGL,WaveBrFraSC,&

> main_pass.fsc

        IF(.NOT.ALLOCATED(WaveBrFraSW)) ALLOCATE(WaveBrFraSW(MXC,MYC))
        IF(.NOT.ALLOCATED(WaveBrFraSC)) ALLOCATE(WaveBrFraSC(Mloc,Nloc))

        WaveBrFraSW = ZERO
        WaveBrFraSC = ZERO

        IF(.NOT.ALLOCATED(WaveBrFraGL))ALLOCATE(WaveBrFraGL(MXCGL,MYCGL))

        WaveBrFraGL = ZERO

          WaveBrFraSW(IX,IY) = COMPDA(INDX,JQB)

! wave breaking fraction
        CALL GATHER_SWAN2GLOBAL(WaveBrFraSW,WaveBrFraGL)
# if defined (ONED_IN_X)
        DO IY=1,MYCGL
        DO IX=1,MXCGL
        WaveBrFraGL(IX,IY)=WaveBrFraGL(IX,INT(MYCGL/2))
        ENDDO
        ENDDO       
# elif defined (ONED_IN_Y)
        DO IY=1,MYCGL
        DO IX=1,MXCGL
        WaveBrFraGL(IX,IY)=WaveBrFraGL(INT(MXCGL/2),IY)
        ENDDO
        ENDDO   
# endif
        CALL DISTRIBUTE2SHORECIRC(WaveBrFraGL,WaveBrFraSC)


! wave breaking fraction
# if defined (ONED_IN_X)
        DO IY=1,MYC
        DO IX=1,MXC
        WaveBrFraSW(IX,IY)=WaveBrFraSW(IX,INT(MYC/2))
        ENDDO
        ENDDO       
# elif defined (ONED_IN_Y)
        DO IY=1,MYC
        DO IX=1,MXC
        WaveBrFraSW(IX,IY)=WaveBrFraSW(INT(MXC/2),IY)
        ENDDO
        ENDDO   
# endif
        CALL DISTRIBUTE2SHORECIRC(WaveBrFraSW,WaveBrFraSC)


>sc_io.fsc

     IF(OUT_Wdis)THEN
        TMP_NAME = TRIM(FDIR)//'Wdis_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,WaveDissSC)
        TMP_NAME = TRIM(FDIR)//'Wbrk_'//TRIM(FILE_NAME)
        call PutFile(TMP_NAME,WaveBrFraSC)
     ENDIF

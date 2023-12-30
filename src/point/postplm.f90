
SUBROUTINE POSTPLM( EMLAYS, S, ZBOT, ZTOP, PRESF, ZZF, TA, ZH, LBOT, LTOP, LFRAC )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !    Subroutine POSTPLM computes plume fractions given a top and bottom
    !    height of the plume.  It assumes a uniform distribution in pressure
    !    (mass concentration -- minor hydrostatic assumption) from bottom to top.
    !
    !  PRECONDITIONS REQUIRED:
    !    Top and bottom of plume as input, vertical grid structure defined,
    !    vertical pressure distribution and temperature provided.
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !       I/O API
    !
    !  REVISION  HISTORY:
    !    Copied from postplm.f v 1.3 in DAQM-V2 Emissions Preprocessor by
    !           M. Houyoux 3/99
    !    Replaced POLY calls with direct hydrostatic calculations by
    !           G. Pouliot 2/05
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***********************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !                System
    ! File: @(#)$Id$
    !
    ! COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
    ! All Rights Reserved
    !
    ! Carolina Environmental Program
    ! University of North Carolina at Chapel Hill
    ! 137 E. Franklin St., CB# 6116
    ! Chapel Hill, NC 27599-6116
    !
    ! smoke@unc.edu
    !
    ! Pathname: $Source$
    ! Last updated: $Date$
    !
    !***********************************************************************
    USE M3UTILIO

    IMPLICIT NONE

    !.......   INCLUDES:
    INCLUDE 'EMCNST3.h90'
    INCLUDE 'CONST3.EXT'    !  physical constants from I/O API

    !.......   SUBROUTINE ARGUMENTS
    INTEGER, INTENT (IN) :: EMLAYS               ! no. emissions layers
    INTEGER, INTENT (IN) :: S                    ! source ID
    REAL   , INTENT (IN) :: ZBOT                 ! plume bottom elevation (m)
    REAL   , INTENT (IN) :: ZTOP                 ! plume top elevation (m)
    REAL   , INTENT (IN) :: PRESF( 0:EMLAYS )    ! pressure at full-levels (mb)
    REAL   , INTENT (IN) :: ZZF  ( 0:EMLAYS )    ! elevation at full-levels (m)
    REAL   , INTENT (IN) :: TA   ( 1:EMLAYS )    ! temperature at half-levels (K)
    REAL   , INTENT (IN) :: ZH   ( 1:EMLAYS )    ! layer center  height (m)
    INTEGER, INTENT(OUT) :: LBOT                 ! plume bottom layer
    INTEGER, INTENT(OUT) :: LTOP                 ! plume top layer
    REAL   , INTENT(OUT) :: LFRAC( EMLAYS )      ! layer fractions for source

    !.......   Local variables

    INTEGER       L

    DOUBLE PRECISION    DDP
    DOUBLE PRECISION    PDIFF

    REAL          PBOT, PTOP
    REAL          TEMP

    CHARACTER(300) MESG

    !***********************************************************************
    !   begin body of subroutine POSTPLM

    !.......   Compute LBOT, LTOP so that
    !.......   ZZF( LBOT-1 ) <= ZBOT < ZZF( LBOT ) and
    !.......   ZZF( LTOP-1 ) <= ZTOP < ZZF( LTOP )

    DO L = 1, EMLAYS - 1
        IF ( ZBOT <= ZZF( L ) ) THEN
            LBOT = L
            GO TO  122       ! end loop and skip reset of LBOT
        ELSE
            LFRAC( L ) = 0.0          ! fractions below plume
        END IF
    END DO
    LBOT = EMLAYS               !  fallback

122 CONTINUE                    !  loop exit:  bottom found at LBOT

    IF ( ZTOP <= ZZF( LBOT ) ) THEN      !  plume in this layer

        LFRAC( LBOT ) = 1.0
        LTOP = LBOT

        DO L = LBOT + 1, EMLAYS      ! fractions above plume
            LFRAC( L ) = 0.0
        END DO

    !.......  Note- this check not in original algorithm, but w/o it,
    !                         can end up with fractions > 1.0
    ELSE IF( LBOT == EMLAYS ) THEN        ! plume above top layer

        LFRAC( LBOT ) = 1.0

        DO L = 1, EMLAYS-1           ! fractions below plume
            LFRAC( L ) = 0.0
        END DO

    ELSE                                   ! plume crosses layers

        DO L = LBOT + 1, EMLAYS
            IF ( ZTOP <= ZZF( L ) ) THEN
                LTOP = L
                GO TO 126      ! end loop and skip reset of LTOP
            END IF
        END DO
        LTOP = EMLAYS

126     CONTINUE

        !.......   Compute corresponding PBOT,PTOP so that
        !.......   PRESF( LBOT-1 ) <= PBOT < PRESF( LBOT ) and
        !.......   PRESF( LTOP-1 ) <= PTOP < PRESF( LTOP )

        !....... If above full layer and below half layer...
        IF( ZBOT < ZH( LBOT ) .AND. ZBOT > ZZF( LBOT-1 ) ) THEN

        !.......  Special case near ground
            IF( ZBOT < ZH( 1 ) ) THEN
                TEMP = TA( 1 )
            ELSE
                TEMP = ( TA( LBOT ) + TA( LBOT-1 ) ) / 2.
            END IF

        !.......  Otherwise, above full layer and above half layer
        ELSE
            TEMP = ( TA( LBOT ) + TA( LBOT+1 ) ) / 2.
        END IF

        !.......  Calculate bottom using hydrostatic assumption
        PBOT = PRESF( LBOT ) * EXP( GRAV / (RDGAS*TEMP) * (ZZF( LBOT ) - ZBOT ))

        !.......  If above full layer and below half layer...
        IF( ZTOP < ZH( LTOP ) .AND. ZTOP > ZZF( LTOP-1 ) ) THEN

            !.......  Special case near ground
            IF( ZTOP < ZH( 1 ) ) THEN
                TEMP = TA( 1 )
            ELSE
                TEMP = ( TA( LTOP ) + TA( LTOP-1 ) ) / 2.
            END IF

        !.......  Otherwise, above full layer and above half layer
        ELSE
            TEMP = ( TA( LTOP ) + TA( LTOP+1 ) ) / 2.
        END IF

        !.......  Calculate top using hydrostatic assumption
        PTOP = PRESF( LTOP-1 ) * EXP( -GRAV / (RDGAS*TEMP) * (ZTOP - ZZF( LTOP-1 )) )

        PDIFF = DBLE( PBOT ) - DBLE( PTOP )

        IF( PDIFF > 0. ) THEN

            DDP = 1.0D0 / ( PDIFF )
            LFRAC( LBOT ) = DDP * (DBLE( PBOT ) - DBLE( PRESF( LBOT ) ))
            LFRAC( LTOP ) = DDP * (DBLE( PRESF( LTOP-1 ) ) - DBLE( PTOP ))

        ELSE
            WRITE( MESG,94010 )    &
             'Infinitely small plume created for source ', S,    &
             CRLF() // BLANK5 //    &
             'because of inverted vertical pressure gradient' //    &
             CRLF() // BLANK5 //    &
             'All emissions put in first layer.'
            CALL M3WARN( 'POSTPLM', 0,0, MESG )

            LBOT = 1
            LTOP = 1
            LFRAC( LBOT ) = 1.0

        ENDIF

        DO L = LBOT+1, LTOP-1     !  layers in plume
            LFRAC( L ) = DDP * (DBLE( PRESF( L-1 ) ) - DBLE( PRESF( L ) ))
        END DO

        DO L = LTOP+1, EMLAYS     !  fractions above plume
            LFRAC( L ) = 0.0
        END DO

    END IF              !  if ztop in same layer as zbot, or not

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I7, :, 1X ) )

END SUBROUTINE POSTPLM


SUBROUTINE PREPLM( EMLAYS, HMIX, HTS, PSFC, TS, DDZF, QV, TA,&
&                   UW, VW, ZH, ZF, ZSTK, PRES, LSTK, LPBL, TSTK,&
&                   WSTK, DTHDZ, WSPD, ZZF )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!    Computes the values needed for the PLMRIS subroutine from the
!    meteorology data.
!
!  PRECONDITIONS REQUIRED:
!    Interpolated (to the location of a source) meteorology data as input,
!    vertical grid structure.
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!       I/O API
!
!  REVISION  HISTORY:
!       Copied from preplm.f v 1.2 in DAQM-V2 Emissions Preprocessor by
!           M. Houyoux 3/99
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
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

!...........   INCLUDES:
    INCLUDE 'CONST3.EXT'    ! physical and mathematical constants

!...........   SUBROUTINE ARGUMENTS (NOTE: All met parms are per-source)
    INTEGER, INTENT (IN) :: EMLAYS          ! no. emissions layers
    REAL   , INTENT (IN) :: HMIX            ! mixing height
    REAL   , INTENT (IN) :: HTS             ! stack height
    REAL   , INTENT (IN) :: PSFC            ! surface pressure
    REAL   , INTENT (IN) :: TS              ! surface temperature
    REAL   , INTENT (IN) :: DDZF ( EMLAYS ) ! 1/( zf(l) - zf(l-1) )
    REAL   , INTENT (IN) :: QV   ( EMLAYS ) ! mixing ratio
    REAL   , INTENT (IN) :: TA   ( EMLAYS ) ! absolute temperature
    REAL   , INTENT (IN) :: UW   ( EMLAYS ) ! x-direction winds
    REAL   , INTENT (IN) :: VW   ( EMLAYS ) ! y-direction winds
    REAL   , INTENT (IN) :: ZF   ( EMLAYS ) ! layer surface height (m)
    REAL   , INTENT (IN) :: ZH   ( EMLAYS ) ! layer center  height (m)
    REAL   , INTENT (IN) :: ZSTK ( EMLAYS ) ! zf( l,s ) - stkht(s) (m)
    REAL   , INTENT (IN) :: PRES ( 0:EMLAYS )! pressure at full-layer heights
    INTEGER, INTENT(OUT) :: LSTK            ! first L: ZF(L) > STKHT
    INTEGER, INTENT(OUT) :: LPBL            ! first L: ZF(L) > mixing layer
    REAL   , INTENT(OUT) :: TSTK            ! tmptr at top of stack (K)
    REAL   , INTENT(OUT) :: WSTK            ! wind speed @ top of stack(m/s)
    REAL   , INTENT(OUT) :: DTHDZ( EMLAYS ) ! potential temp. gradient
    REAL   , INTENT(OUT) :: WSPD ( EMLAYS ) ! wind speed (m/s)
    REAL   , INTENT(OUT) :: ZZF  ( 0:EMLAYS )! elevation at full-levels

!...........   Local variables

    INTEGER       L, M

    REAL          ES
    REAL          QSFC
    REAL          TVSFC
    REAL          THETG
    REAL          THV1
    REAL          THVK
    REAL          TV( EMLAYS )     ! virtual temperature
    REAL          TF( EMLAYS )     ! full-layer height temperatures
    REAL          P, Q
    REAL          DZZ
    REAL          DELZ

!***********************************************************************
!   begin body of subroutine PREPLM

!.........  Compute wind speed and virtual temperature
    DO L = 1, EMLAYS

        P = UW( L )
        Q = VW( L )
        WSPD( L ) = SQRT( P * P  +  Q * Q )
        TV  ( L ) = TA( L ) *&
        &            ( 1. + 0.622 * ( QV( L ) / ( 1. + QV( L ) ) ) )

    END DO

    ES    = 6.1078 * EXP( 5384.21 / CTOK - 5384.21 / TS )
    QSFC  = 0.622  * ES / ( PSFC - ES )
    TVSFC = TS   * ( 1.0 + 0.6077 * QSFC )
    THETG = TVSFC  * ( 1000.0 / PSFC )**0.286

    IF ( HMIX .LE. ZF( 1 ) ) THEN
        LPBL = 1
    END IF
    IF ( HTS .LE. ZF( 1 ) ) THEN
        LSTK = 1
    END IF

    ZZF( 0 ) = 0.0
    ZZF( 1 ) = ZF( 1 )

!.........  Compute temperatures at full-layer face heights
    DO L = 1, EMLAYS - 1
        DELZ = ZH( L+1 ) - ZH( L )
        TF( L ) = TV( L ) + ( TV( L+1 ) - TV( L ) ) *&
        &          ( ZF( L ) - ZH( L ) ) / DELZ
    END DO

    DELZ = ZH( EMLAYS ) - ZH( EMLAYS-1 )
    TF( EMLAYS ) = TV( EMLAYS ) - (TV( EMLAYS-1 ) - TV( EMLAYS )) *&
    &               ( ZF( EMLAYS ) - ZH( EMLAYS ) ) / DELZ

    THV1 = TF( 1 ) * ( 1000. / PRES( 1 ) )**0.286
    DTHDZ( 1 ) = ( THV1 - THETG ) / ZF( 1 )

    DO L = 2, EMLAYS

        IF ( HMIX > ZF( L-1 ) )  LPBL = L
        IF ( HTS > ZF( L-1 ) )  LSTK = L

        THVK = TF( L ) * ( 1000. / PRES( L ) )**0.286
        DTHDZ( L ) = DDZF( L ) * ( THVK - THV1 )
        THV1 = THVK

        ZZF( L ) = ZF( L )

    END DO

!.........  Set the 1st level vertical THETV gradient to the 2nd layer value
!           This overrides the layer 1 gradient determined above
    DTHDZ( 1 ) = DTHDZ( 2 )

    M    = MAX( 1, LSTK - 2 )
    TSTK =      POLY( HTS, ZH( M ), TA( M ), 3 )
    WSTK = MAX( POLY( HTS, ZH( M ), WSPD( M ), 3 ), 0.1 )

    RETURN

END SUBROUTINE PREPLM

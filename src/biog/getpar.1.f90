SUBROUTINE GETPAR( RSOLAR, PRES, ZEN, PARDB, PARDIF )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!
!        Based on code from Bart Brashers (10/2000), which was based on
!        code from Weiss and Norman (1985).
!
!
!  PRECONDITIONS REQUIRED:
!     Solar radiation (W/m2) and pressure (mb)
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!    3/01 : Prototype by JMV
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

!........ Inputs

    REAL , INTENT  (IN) :: RSOLAR   ! modeled or observed total radiation (W/m2)
    REAL , INTENT  (IN) :: PRES     ! atmospheric pressure (mb)
    REAL , INTENT  (IN) :: ZEN      ! solar zenith angle (radians)

!........ Outputs

    REAL, INTENT (OUT) :: PARDB     ! direct beam PAR( umol/m2-s)
    REAL, INTENT (OUT) :: PARDIF    ! diffuse PAR ( umol/m2-s)

!...........   PARAMETERS and their descriptions:

    REAL, PARAMETER :: WATT2UMOL = 4.6  ! convert W/m^2 to umol/m^2-s (4.6)

!
    REAL RATIO              ! transmission fraction for total radiation
    REAL OT                   ! optical thickness
    REAL RDVIS                ! possible direct visible beam (W/m^2)
    REAL RFVIS              ! possible visible diffuse (W/m^2)
    REAL WA                 ! water absorption in near-IR (W/m^2)
    REAL RDIR               ! direct beam in near-IR (W/m^2)
    REAL RFIR               ! diffuse near-IR (W/m^2)
    REAL RVT                ! total possible visible radiation (W/m^2)
    REAL RIRT               ! total possible near-IR radiation (W/m^2)
    REAL FVIS               ! fraction of visible to total
    REAL FVB                ! fraction of visible that is direct beam
    REAL FVD                ! fraction of visible that is diffuse

    CHARACTER(16) :: PROGNAME = 'GETPAR'   !  program name

    CHARACTER(256)  MESG

!***************************************
!   begin body of subroutine

!............ Assume that PAR = 0 if zenith angle is greater than 87 degrees
!............ or if solar radiation is zero

    IF (ZEN .GE. 1.51844 .OR. RSOLAR .LE. 0.) THEN
        PARDB  = 0.
        PARDIF = 0.
        RETURN
    ENDIF

!............ Compute clear sky (aka potential) radiation terms

    OT    = PRES / 1013.25 / COS(ZEN)              !Atmospheric Optical thickness
    RDVIS = 600. * EXP(-.185 * OT) * COS(ZEN)      !Direct visible beam, eqn (1)
    RFVIS = 0.42 * (600 - RDVIS) * COS(ZEN)        !Visible Diffuse, eqn (3)
    WA    = 1320 * .077 * (2. * OT)**0.3           !water absorption in near-IR, eqn (6)
    RDIR  = (720. * EXP(-0.06 * OT)-WA) * COS(ZEN) !Direct beam near-IR, eqn (4)
    RFIR  = 0.65 * (720. - WA - RDIR) * COS(ZEN)   !Diffuse near-IR, eqn (5)

    RVT   = RDVIS + RFVIS                    !Total visible radiation, eqn (9)
    RIRT  = RDIR + RFIR                      !Total near-IR radiation, eqn (10)
    FVIS  = RVT/(RIRT + RVT)                 !Fraction of visible to total radiation, eqn 7
    RATIO = RSOLAR /(RIRT + RVT)             !Ratio of "actual" to clear sky solar radiation

!............ Compute fraction of visible that is direct beam

    IF (RATIO .GE. 0.89) THEN
        FVB = RDVIS/RVT * 0.941124
    ELSE IF (RATIO .LE. 0.21) THEN
        FVB = RDVIS/RVT * 9.55E-3
    ELSE
        FVB = RDVIS/RVT * (1.-((0.9 - RATIO)/0.7)**0.666667)
    ENDIF
    FVD = 1. - FVB

!............ Compute PAR (direct beam and diffuse) in umol/m2-sec

    PARDB  = RSOLAR * FVIS * FVB * WATT2UMOL
    PARDIF = RSOLAR * FVIS * FVD * WATT2UMOL


    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Informational (LOG) message formats... 92xxx


!...........   Internal buffering formats............ 94xxx

END SUBROUTINE GETPAR


SUBROUTINE GETPARB( RSOLAR, PRES, COSZ, PARDB, PARDIF )

!-----------------------------------------------------------------------
! Description:
!   Compute direct and diffuse photosynthetically active radiation (PAR).
!   Based on code from Bart Brashers (10/2000), which was based on
!   code from Weiss and Norman (1985).

! Preconditions:
!   Solar radiation (W/m2) and pressure (mb)

! Subroutines and Functions Called:

! Revision History:
!   3/01 Prototype by JMV
!  10/17 J.Young: rename, optimize
!-----------------------------------------------------------------------
! Modified from:

! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling System
! File: @(#)$Id$
! COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
! All Rights Reserved
! Carolina Environmental Program
! University of North Carolina at Chapel Hill
! 137 E. Franklin St., CB# 6116
! Chapel Hill, NC 27599-6116
! smoke@unc.edu
! Pathname: $Source$
! Last updated: $Date$
!-----------------------------------------------------------------------

    IMPLICIT NONE

! Arguments:
    REAL, INTENT ( IN ) :: RSOLAR   ! modeled or observed total radiation [W/m2]
    REAL, INTENT ( IN ) :: PRES     ! atmospheric pressure [mb]
    REAL, INTENT ( IN ) :: COSZ     ! cosine of solar zenith angle
    REAL, INTENT( OUT ) :: PARDB    ! direct beam PAR [umol/m2-s]
    REAL, INTENT( OUT ) :: PARDIF   ! diffuse PAR [umol/m2-s]

! Parameters:
    REAL,             PARAMETER :: WATT2UMOL = 4.6  ! convert W/m^2 to umol/m^2-s
    REAL,             PARAMETER :: TWOTHD    = 2.0 / 3.0
    CHARACTER(  16 ), PARAMETER :: PNAME     = 'GETPAR'   ! procedure name

! Local Variables:
    REAL RATIO              ! transmission fraction for total radiation
    REAL OT                 ! optical thickness
!     REAL RDVIS              ! possible direct visible beam [(W/m^2]
!     REAL RFVIS              ! possible visible diffuse [(W/m^2]
    REAL WA                 ! water absorption in near-IR [(W/m^2]
!     REAL RDIR               ! direct beam in near-IR [(W/m^2]
!     REAL RFIR               ! diffuse near-IR [(W/m^2]
!     REAL RVT                ! total possible visible radiation [(W/m^2]
!     REAL RIRT               ! total possible near-IR radiation [(W/m^2]
!     REAL FVIS               ! fraction of visible to total
!     REAL FVB                ! fraction of visible that is direct beam
!     REAL FVD                ! fraction of visible that is diffuse
    REAL FAC                ! direct beam factor / combination factor

    REAL A, B, C, DEN, RW   ! replacement composite variables

!-----------------------------------------------------------------------

!   Original implementation ............................

! (ZEN was originally an argument, but replaced by COS(ZEN) = COSZ)

! Assume that PAR = 0 if zenith angle is greater than 87 degrees (1.51844 radians)
! or if solar radiation is zero
!     IF ( ZEN .GE. 1.51844 .OR. RSOLAR .LE. 0.0 ) THEN
!        PARDB  = 0.0
!        PARDIF = 0.0
!        RETURN
!     END IF

! Compute clear sky (aka potential) radiation terms
!     OT    = PRES / 1013.25 / COSZ                ! Atmospheric Optical thickness
!     RDVIS = 600.0 * EXP( -0.185 * OT ) * COSZ    ! Direct visible beam, eqn (1)
!     RFVIS = 0.42 * ( 600.0 - RDVIS ) * COSZ      ! Visible Diffuse, eqn (3)
!     WA    = 1320.0 * 0.077 * ( 2.0 * OT ) ** 0.3 ! water absorption in near-IR, eqn (6)
!     RDIR  = ( 720.0 * EXP( -0.06 * OT ) - WA ) * COSZ ! Direct beam near-IR, eqn (4)
!     RFIR  = 0.65 * ( 720.0 - WA - RDIR ) * COSZ  ! Diffuse near-IR, eqn (5)

!     RVT   = RDVIS + RFVIS            ! Total visible radiation, eqn (9)
!     RIRT  = RDIR + RFIR              ! Total near-IR radiation, eqn (10)
!     FVIS  = RVT / ( RIRT + RVT )     ! Fraction of visible to total radiation, eqn 7
!     RATIO = RSOLAR / ( RIRT + RVT )  ! Ratio of "actual" to clear sky solar radiation

! Compute fraction of visible that is direct beam
!     IF ( RATIO .GE. 0.89 ) THEN
!        FVB = RDVIS / RVT * 0.941124
!     ELSE IF ( RATIO .LE. 0.21 ) THEN
!        FVB = RDVIS / RVT * 9.55E-3
!     ELSE
!        FVB = RDVIS / RVT * ( 1.0 - ( ( 0.9 - RATIO ) / 0.7 ) ** 0.666667 )
!     END IF
!     FVD = 1.0 - FVB

! Compute PAR (direct beam and diffuse) in umol/m2-sec
!     PARDB  = RSOLAR * FVIS * FVB * WATT2UMOL
!     PARDIF = RSOLAR * FVIS * FVD * WATT2UMOL

!   New implementation .................................

    IF ( COSZ .LE. 0.052336 .OR. RSOLAR .LE. 0.0 ) THEN
        PARDB  = 0.0
        PARDIF = 0.0
        RETURN
    END IF

    OT = PRES / ( 101325.0 * COSZ )  ! PRES in Pa
    A = 600.0 * EXP( -0.185 * OT )
    WA = 125.1335182 * ( OT ) ** 0.3
    B = 720.0 * EXP( -0.06 * OT ) - WA
    C = ( 1.0 - 0.42 * COSZ ) * A
    DEN = 720.0 + C - 0.65 * WA + ( 1.0 - 0.65 * COSZ ) * B
    RATIO = RSOLAR / ( DEN * COSZ )
    IF ( RATIO .GE. 0.89 ) THEN
        FAC = 0.941124
    ELSE IF ( RATIO .LE. 0.21 ) THEN
        FAC = 9.55E-3
    ELSE
        FAC = 1.0 - ( 1.42857143 * ( 0.9 - RATIO ) ) ** TWOTHD
    END IF
    RW = RATIO * WATT2UMOL * COSZ
    PARDB  = RW * FAC * A
    PARDIF = RW * ( 252.0 + C ) - PARDB

    RETURN

END SUBROUTINE GETPARB


REAL FUNCTION PLUMRIS( HS, TS, VS, DS )

!***********************************************************************
!  subroutine body starts at line  82
!
!  DESCRIPTION:  Computes effective plume height using Briggs algorithm.  See:
!   1) Briggs, Gary A., 1971: Some Recent Analyses of Plume Rise Observation
!      pp 1029 - 1032 in PROCEEDINGS OF THE SECOND INTERNATIONAL CLEAN AIR
!      CONGRESS, edited by H. M. Englun and W. T. Beery. Academic Press,
!      New York.
!   2) Briggs, Gary A., 1972: Discussion on Chimney Plumes in Neutral
!      and Stable Surroundings. ATMOS. ENVIRON. 6, 507 - 510. (Jul 72).
!
!  REVISION  HISTORY:
!       Copied 8/99 from plumris.F v4.2 in SMOKE prototype
!       Adapted 10/95 by Carlie J. Coats, Jr., from UAM/EPS BEH072().
!       Discontinuity in EPS DISTF calculation at F=55 resolved so
!       that DISTF --> 595.0  as  F --> 55.0 from either side.
!       Uses standard conditions temperature T=293 deg K,
!                                pressure    P=960 mb,
!                                wind speed  U=  2 m/s,
!                       Pasquill stability KST=  2
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!***********************************************************************
!
! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
! *                System
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
!************************************************************************
    USE M3UTILIO

    IMPLICIT NONE

!...........   INCLUDES:

    INCLUDE 'CONST3.EXT'      ! physical and mathematical constants


!...........   ARGUMENTS and their descriptions:

    REAL      HS    !  physical stack height (m)
    REAL      TS    !  stack gas temperature (deg k)
    REAL      VS    !  stack gas exit velocity (m/sec)
    REAL      DS    !  inside stack diameter (m)

!...........   PARAMETERS and their descriptions:

    REAL        T           !  default ambient air temperature (deg k)
    REAL        P           !  default ambient air pressure (mb)
    REAL        U           !  default wind speed (m/sec)
    REAL        D3          !  one-third

    PARAMETER ( T    = 293.0 ,&
    &            P    = 960.0 ,&
    &            U    =   2.0 ,&
    &            D3   =   1.0 / 3.0 )

!...........   LOCAL VARIABLES:

    REAL      F     ! buoyancy flux (m**4/sec**3)
!       REAL      DELHF ! final plume rise (m)
!       REAL      DISTF ! distance of final plume rise from source (m)


!***********************************************************************
!   begin body of function PLUMRIS

    IF ( TS .LE. T ) THEN
        PLUMRIS = MAX( HS, 3.0 )
        RETURN
    END IF

    IF ( DS .LE. 0.0 ) DS = 0.2
    IF ( HS .LE. 0.0 ) HS = 3.0
    IF ( VS .LE. 0.0 ) VS = 0.5

    F = 0.25 * GRAV * VS * DS * DS * ( TS - T ) / TS

!        IF( F .LT. 55.0 ) THEN
!            DISTF = 3.5 * 13.8906395 * F**0.625
!        ELSE
!            DISTF = 3.5 * 34.22187854 * F**0.4
!        END IF
!        DELHF = ( 1.6 / U ) * ( F * DISTF * DISTF )**D3
!        PLUMRIS = HS + DELHF

    IF( F .LT. 55.0 ) THEN
        PLUMRIS = HS + 21.31311057 * F**0.75 / U
    ELSE
        PLUMRIS = HS + 38.87776061 * F**0.6 / U
    END IF

    RETURN
END FUNCTION PLUMRIS


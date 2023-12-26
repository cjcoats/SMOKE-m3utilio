
MODULE MODBEIS3

!***********************************************************************
!  Module body starts at line 42
!
!  DESCRIPTION:
!     This module contains the public variables and allocatable arrays
!     used in calculating biogenic emissions using BEIS v3.09 or v3.12.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     03/01: prototype by Jeff Vukovich
!
!***************************************************************************
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

    IMPLICIT NONE

!...........   Emission factor, vegetation types tables:

    INTEGER, PUBLIC ::    NVEG                     !  Number of veg types
    INTEGER, ALLOCATABLE, PUBLIC :: LAI( : )       !  Leaf area index
    REAL,    ALLOCATABLE, PUBLIC :: EMFAC( :, : )  !  Emission factors
    REAL,    ALLOCATABLE, PUBLIC :: LFBIO( : )     !  Dry leaf biomass
    REAL,    ALLOCATABLE, PUBLIC :: WFAC( : )      !  Winter biomass factor
    REAL,    ALLOCATABLE, PUBLIC :: SLW( : )       !  Specific leaf weight

    REAL,    ALLOCATABLE, PUBLIC :: AVGISOP( :, :, : )   ! avg isoprene
    REAL,    ALLOCATABLE, PUBLIC :: AVGOVOC( :, :, : )   ! avg other VOCs
    REAL,    ALLOCATABLE, PUBLIC :: AVGMONO( :, :, : )   ! avg monoterpenes
    REAL,    ALLOCATABLE, PUBLIC :: AVGNO( :, :, : )     ! avg nitric oxide
    REAL,    ALLOCATABLE, PUBLIC :: AVGLAI( :, :, :, : ) ! avg leaf index

    REAL,    ALLOCATABLE, PUBLIC :: AVGEMIS( :, :, :, : )  ! avg emissions (3.12)
    REAL,    ALLOCATABLE, PUBLIC ::  NOEMIS( :, :, : )     ! NO emissions (3.12)

    CHARACTER(16), ALLOCATABLE, PUBLIC :: VEGID( : )     !  Veg types

END MODULE MODBEIS3

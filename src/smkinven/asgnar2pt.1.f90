
SUBROUTINE ASGNAR2PT

!***********************************************************************
!  subroutine body starts at line 109
!
!  DESCRIPTION:
!      For each source, find the most specific area-to-point adjustment
!      that applies to that source. Do this using the grouped tables of
!      gridding cross references from RDAR2PT (and stored in MODLISTS).
!      Only one group (full FIPS & full SCC) has been defined at this time
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 11/02 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!************************************************************************
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
!***************************************************************************
    USE M3UTILIO

!...........   MODULES for public variables
!...........   This module contains the source ararys
    USE MODSOURC, ONLY: CSOURC, CSCC

!...........   This module contains the cross-reference tables
    USE MODXREF, ONLY: AR2PTTBL, AR2PTIDX, AR2PTCNT, TXCNT,&
    &                   CHRT09, CHRT08, ARPT09, ARPT08

!.........  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, NSRC, LSCCEND

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   EXTERNAL FUNCTIONS and their descriptions:

    LOGICAL, EXTERNAL :: SETSCCTYPE

!.........  Other local variables
    INTEGER          S           !  counters and indices

    INTEGER          IOS         ! i/o status
    INTEGER          F4, F5      ! tmp find indices

    LOGICAL       :: EFLAG    = .FALSE.
    LOGICAL          SCCFLAG           ! true: SCC type is different from previous

    CHARACTER(256)        MESG     ! message buffer
    CHARACTER(SRCLEN3)    CSRC     ! tmp source chars string
    CHARACTER(FIPLEN3)    CFIP     ! tmp (character) FIPS code
    CHARACTER(FPSLEN3)    CFIPSCC  ! tmp FIPS code // SCC
    CHARACTER(FPSLEN3)    CFIPSL   ! tmp FIPS code // left SCC
    CHARACTER(SCCLEN3)    TSCC     ! tmp 10-digit SCC
    CHARACTER(SCCLEN3)    TSCCL    ! tmp left digits of TSCC

    CHARACTER(16) :: PROGNAME = 'ASGNAR2PT' ! program name

!***********************************************************************
!   begin body of subroutine ASGNAR2PT

!.........  Allocate memory for source-based arrays
    ALLOCATE( AR2PTTBL( NSRC ),&
    &          AR2PTIDX( NSRC ),&
    &          AR2PTCNT( NSRC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'AR2PTTBL...AR2PTCNT', PROGNAME )

!.........  Initialize arrays
    AR2PTTBL = 0  ! array
    AR2PTIDX = 0  ! array
    AR2PTCNT = 0  ! array

!.........  Loop through the sorted sources
    DO S = 1, NSRC

! NOTE: Perhaps change this to be a generic routine for both cross-reference need?  The
! N: only problem with this is that CHRT* will be shared by both steps :(

!.............  Create selection
        SELECT CASE ( CATEGORY )

          CASE ( 'AREA', 'MOBILE' )
            CSRC    = CSOURC( S )
            TSCC    = CSCC( S )

!.................  Set type of SCC
            SCCFLAG = SETSCCTYPE( TSCC )
            TSCCL   = TSCC( 1:LSCCEND )

            CFIP    = CSRC( 1:FIPLEN3 )
            CFIPSCC = CFIP // TSCC
            CFIPSL  = CFIP // TSCCL
        END SELECT

!.................  Try for FIPS code & SCC match
        F5 = FINDC( CFIPSCC, TXCNT( 9 ), CHRT09 )
        F4 = FINDC( CFIPSL , TXCNT( 8 ), CHRT08 )

        IF( F5 .GT. 0 ) THEN
            AR2PTTBL( S ) = ARPT09( F5,1 )
            AR2PTIDX( S ) = ARPT09( F5,2 )
            AR2PTCNT( S ) = ARPT09( F5,3 )

        ELSE IF( F4 .GT. 0 ) THEN
            AR2PTTBL( S ) = ARPT08( F4,1 )
            AR2PTIDX( S ) = ARPT08( F4,2 )
            AR2PTCNT( S ) = ARPT08( F4,3 )

        END IF

    END DO        !  end loop on sources

    RETURN

END SUBROUTINE ASGNAR2PT

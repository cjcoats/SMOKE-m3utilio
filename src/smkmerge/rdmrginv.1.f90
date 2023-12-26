
SUBROUTINE RDMRGINV

!***********************************************************************
!  subroutine RDMRGINV body starts at line
!
!  DESCRIPTION:
!      The purpose of this subroutine is to read in inventory variables needed
!      for the merge program
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 8/99 by M. Houyoux
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
!****************************************************************************
    USE M3UTILIO

!.........  MODULES for public variables
!.........  This module contains the major data structure and control flags
    USE MODMERGE, ONLY: AFLAG,  MFLAG,  PFLAG,      &! source flags by category
    &                    AENAME, MENAME, PENAME,     &! inventory file names
    &                    ASDEV,  MSDEV,  PSDEV,      &! inventory file names
    &                    AIFIP,  MIFIP,  PIFIP,      &! country/state/county codes
    &                    NASRC,  NMSRC,  NPSRC       ! no. of sources

!.........  This module contains the lists of unique source characteristics
    USE MODLISTS, ONLY: NINVIFIP, INVCFIP

!...........   This module contains the source arrays
    USE MODSOURC, ONLY: CIFIP

    IMPLICIT NONE

!...........   INCLUDES:

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
    INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

!.........  Local allocatable parameters
    CHARACTER(FIPLEN3), ALLOCATABLE :: AFIPS ( : )  ! area source co/st/cy codes
    CHARACTER(FIPLEN3), ALLOCATABLE :: MFIPS ( : )  ! mobile source co/st/cy codes
    CHARACTER(FIPLEN3), ALLOCATABLE :: PFIPS ( : )  ! point source co/st/cy codes
    CHARACTER(FIPLEN3), ALLOCATABLE :: TMPFIP( : )  ! all unsort source co/st/cy codes
    INTEGER, ALLOCATABLE :: SIDX  ( : )  ! sorting index

!.........  Local parameters
    INTEGER       J, K, L, N, S      ! pointers and counters

    INTEGER       IOS          ! i/o status
    INTEGER       MXNFIP       ! max no. FIPS codes

    INTEGER    :: NAFIP = 0
    INTEGER    :: NMFIP = 0
    INTEGER    :: NPFIP = 0

    CHARACTER(FIPLEN3)    PFIP         ! tmp fip from previous iteration

    CHARACTER(33), PARAMETER :: PART1 =&
    &                           'Error reading variable CIFIP from '
    CHARACTER(15), PARAMETER :: PART3 = ' INVENTORY file'

!...........   Other local variables

    CHARACTER(300)   MESG     ! message buffer

    CHARACTER(16) :: PROGNAME = 'RDMRGINV' ! program name

!***********************************************************************
!   begin body of subroutine RDMRGINV

!.........  Read in area source country, state, county code
    IF( AFLAG ) THEN

        ALLOCATE( AIFIP( NASRC ), STAT=IOS )   ! country/state/county codes
        CALL CHECKMEM( IOS, 'AIFIP', PROGNAME )

        CALL RDINVCHR( 'AREA', AENAME, ASDEV, NASRC, 1, 'CIFIP' )
        AIFIP = CIFIP

!.............  Count area source country, state, and county codes
        CALL COUNT_FIPS( NASRC, AIFIP, NAFIP )

!.............  Allocate memory for codes
        ALLOCATE( AFIPS( NAFIP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AFIPS', PROGNAME )

!.............  Create temporary area-source FIPs code list
        CALL CREATE_FIPS( NASRC, NAFIP, AIFIP, AFIPS )

    END IF


!.........  Read in mobile source country, state, county code
    IF( MFLAG ) THEN
        ALLOCATE( MIFIP( NMSRC ), STAT=IOS )   ! country/state/county codes
        CALL CHECKMEM( IOS, 'MIFIP', PROGNAME )

        CALL RDINVCHR( 'MOBILE', MENAME, MSDEV, NMSRC, 1, 'CIFIP' )
        MIFIP = CIFIP

!.............  Count mobile source country, state, and county codes
        CALL COUNT_FIPS( NMSRC, MIFIP, NMFIP )

!.............  Allocate memory for codes
        ALLOCATE( MFIPS( NMFIP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MFIPS', PROGNAME )

!.............  Create temporary mobile-source FIPs code list
        CALL CREATE_FIPS( NMSRC, NMFIP, MIFIP, MFIPS )
    END IF

!.........  Create temporary mobile-source FIPs code list

!.........  Read in point source country, state, county code
    IF( PFLAG ) THEN
        ALLOCATE( PIFIP( NPSRC ), STAT=IOS )   ! country/state/county codes
        CALL CHECKMEM( IOS, 'PIFIP', PROGNAME )

        CALL RDINVCHR( 'POINT', PENAME, PSDEV, NPSRC, 1, 'CIFIP' )
        PIFIP = CIFIP

!.............  Count point source country, state, and county codes
        CALL COUNT_FIPS( NPSRC, PIFIP, NPFIP )

!.............  Allocate memory for codes
        ALLOCATE( PFIPS( NPFIP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PFIPS', PROGNAME )

!.............  Create temporary point-source FIPs code list
        CALL CREATE_FIPS( NPSRC, NPFIP, PIFIP, PFIPS )

    END IF

!.........  Combine temporary lists of FIPs codes into single sorted list for
!           reading state and counties file
    MXNFIP = NAFIP + NMFIP + NPFIP
    ALLOCATE(  TMPFIP( MXNFIP ),&
    &          INVCFIP( MXNFIP ),&
    &             SIDX( MXNFIP ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SIDX', PROGNAME )

!.........  Add area sources to complete FIPS list
    N = 0
    DO K = 1, NAFIP
        N = N + 1
        SIDX  ( N ) = N
        TMPFIP( N ) = AFIPS( K )
    END DO

!.........  Add mobile sources to complete FIPS list
    DO K = 1, NMFIP
        N = N + 1
        SIDX  ( N ) = N
        TMPFIP( N ) = MFIPS( K )
    END DO

!.........  Add point sources to complete FIPS list
    DO K = 1, NPFIP
        N = N + 1
        SIDX  ( N ) = N
        TMPFIP( N ) = PFIPS( K )
    END DO

    MXNFIP = N

!.........  Sort all FIPS codes
    CALL SORTIC( MXNFIP, SIDX, TMPFIP )

!.........  Store unique list across all source categories
    N = 0
    PFIP = ' '
    DO K = 1, MXNFIP

        J = SIDX( K )
        IF( TMPFIP( J ) .NE. PFIP ) THEN
            N = N + 1
            INVCFIP( N ) = TMPFIP( J )
        END IF

        PFIP = TMPFIP( J )

    END DO

    NINVIFIP = N

!.........  Deallocate local memory
    IF( ALLOCATED( AFIPS ) ) DEALLOCATE( AFIPS )
    IF( ALLOCATED( MFIPS ) ) DEALLOCATE( MFIPS )
    IF( ALLOCATED( PFIPS ) ) DEALLOCATE( PFIPS )
    IF( ALLOCATED( TMPFIP ) ) DEALLOCATE( TMPFIP )
    IF( ALLOCATED( SIDX ) ) DEALLOCATE( SIDX )

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats.............94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

!*****************  INTERNAL SUBPROGRAMS  ******************************

CONTAINS

!.............  This subprogram counts the FIPS codes in the sorted
!               input array of codes
    SUBROUTINE COUNT_FIPS( NSRC, CIFIP, NFIP )

!.............  Subprogram arguments
        INTEGER,            INTENT (IN) :: NSRC
        CHARACTER(FIPLEN3), INTENT (IN) :: CIFIP( NSRC )
        INTEGER,            INTENT(OUT) :: NFIP

!.............  Local variables
        INTEGER  K

        CHARACTER(FIPLEN3) PFIP   ! from previous iteration

!----------------------------------------------------------------------

        K = 0
        PFIP = ' '
        DO S = 1, NSRC

            IF( CIFIP( S ) .NE. PFIP ) K = K + 1
            PFIP = CIFIP( S )

        END DO
        NFIP = K

    END SUBROUTINE COUNT_FIPS

!....................................................................
!....................................................................

!.............  This subprogram creates a list of unique FIPS codes
!               from the sorted input list
    SUBROUTINE CREATE_FIPS( NSRC, NFIP, CIFIP, CFIPS )

!.............  Subprogram arguments
        INTEGER,            INTENT (IN) :: NSRC
        INTEGER,            INTENT (IN) :: NFIP
        CHARACTER(FIPLEN3), INTENT (IN) :: CIFIP( NSRC )
        CHARACTER(FIPLEN3), INTENT(OUT) :: CFIPS( NFIP )

!.............  Local variables
        INTEGER  K

        CHARACTER(FIPLEN3)  PFIP   ! from previous iteration

!----------------------------------------------------------------------

        K = 0
        PFIP = ' '
        DO S = 1, NSRC

            IF( CIFIP( S ) .NE. PFIP ) THEN
                K = K + 1
                CFIPS( K ) = CIFIP( S )
            END IF

            PFIP = CIFIP( S )

        END DO

    END SUBROUTINE CREATE_FIPS

END SUBROUTINE RDMRGINV

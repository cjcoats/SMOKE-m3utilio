
SUBROUTINE FILLTTBL( NIPOL, NXREF, ICSIZE, XTYPE, XTCNT )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine populates the temporal profile codes part of the grouped
!      temporal cross-reference tables.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!     Created 3/99 by M. Houyoux
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
!***************************************************************************

!...........   This module is for cross reference tables
    USE MODXREF, ONLY: MPRT01, WPRT01, DPRT01,&
    &                   MPRT02, WPRT02, DPRT02,&
    &                   MPRT03, WPRT03, DPRT03,&
    &                   MPRT04, WPRT04, DPRT04,&
    &                   MPRT05, WPRT05, DPRT05,&
    &                   MPRT06, WPRT06, DPRT06,&
    &                   MPRT07, WPRT07, DPRT07,&
    &                   MPRT08, WPRT08, DPRT08,&
    &                   MPRT09, WPRT09, DPRT09,&
    &                   MPRT10, WPRT10, DPRT10,&
    &                   MPRT11, WPRT11, DPRT11,&
    &                   MPRT12, WPRT12, DPRT12,&
    &                   MPRT13, WPRT13, DPRT13,&
    &                   MPRT14, WPRT14, DPRT14,&
    &                   MPRT15, WPRT15, DPRT15,&
    &                   MPRT16, WPRT16, DPRT16,&
    &                   MPRNA, WPRNA, DPRNA,&
    &                   ADDPS, INDXTA, ISPTA

    IMPLICIT NONE

!...........   SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN) :: NIPOL           ! no. pollutants
    INTEGER     , INTENT (IN) :: NXREF           ! no. ungrpd x-ref entries
    INTEGER     , INTENT (IN) :: ICSIZE( * )     ! size of x-ref groups
    INTEGER     , INTENT (IN) :: XTYPE ( NXREF ) ! group no. of x-ref entry
    INTEGER     , INTENT (IN) :: XTCNT ( NXREF ) ! pos. in x-ref group

!...........   Other local variables
    INTEGER       I, J, K, T       ! counter and indices
    INTEGER       IDIUPS           ! tmp w/ pollutant-specific indicator
    INTEGER       IMON, IWEK, IDIU ! temporary temporal profile codes
    INTEGER       ISP              ! temporary pollutant position in EINAM
    INTEGER       TDIM             ! temporary table dimension

    CHARACTER(16) :: PROGNAME = 'FILLTTBL' ! program name

!***********************************************************************
!   begin body of subroutine FILLTTBL

!.........  Store the temporal profile codes for each x-ref entry, depending
!           on the group (XTYPE) and the position in that group (XTCNT)

    DO I = 1, NXREF

        J      = INDXTA( I )
        ISP    = ISPTA ( J )
        IMON   = MPRNA ( J )
        IWEK   = WPRNA ( J )
        IDIU   = DPRNA ( J )
        IDIUPS = IDIU + ADDPS  ! for pollutant specific

        T      = XTYPE ( I )
        K      = XTCNT ( I )

!.............  Skip x-ref because it is invalid or duplicate
        IF( T .EQ. 0 ) CYCLE

        TDIM = ICSIZE( T )
!
!.............  Populate tables depending on type. Note that the pollutant-
!               specific entries are assumed to always come after the
!               non-specific ones (based on the previous sorting).
!.............  The temporal pol-specific entries are stored by adding 90000
!               to the monthly profile number (which has a maximum of 3
!               digits) so that the pol-specific can be identified later
        SELECT CASE ( T )

          CASE( 1 )
            K = 1
            CALL SET_TCODES( 1, NIPOL, MPRT01, WPRT01, DPRT01 )

          CASE( 2 )
            CALL SET_TCODES( TDIM, NIPOL, MPRT02, WPRT02, DPRT02 )

          CASE( 3 )
            CALL SET_TCODES( TDIM, NIPOL, MPRT03, WPRT03, DPRT03 )

          CASE( 4 )
            CALL SET_TCODES( TDIM, NIPOL, MPRT04, WPRT04, DPRT04 )

          CASE( 5 )
            CALL SET_TCODES( TDIM, NIPOL, MPRT05, WPRT05, DPRT05 )

          CASE( 6 )
            CALL SET_TCODES( TDIM, NIPOL, MPRT06, WPRT06, DPRT06 )

          CASE( 7 )
            CALL SET_TCODES( TDIM, NIPOL, MPRT07, WPRT07, DPRT07 )

          CASE( 8 )
            CALL SET_TCODES( TDIM, NIPOL, MPRT08, WPRT08, DPRT08 )

          CASE( 9 )
            CALL SET_TCODES( TDIM, NIPOL, MPRT09, WPRT09, DPRT09 )

          CASE( 10 )
            CALL SET_TCODES( TDIM, NIPOL, MPRT10, WPRT10, DPRT10 )

          CASE( 11 )
            CALL SET_TCODES( TDIM, NIPOL, MPRT11, WPRT11, DPRT11 )

          CASE( 12 )
            CALL SET_TCODES( TDIM, NIPOL, MPRT12, WPRT12, DPRT12 )

          CASE( 13 )
            CALL SET_TCODES( TDIM, NIPOL, MPRT13, WPRT13, DPRT13 )

          CASE( 14 )
            CALL SET_TCODES( TDIM, NIPOL, MPRT14, WPRT14, DPRT14 )

          CASE( 15 )
            CALL SET_TCODES( TDIM, NIPOL, MPRT15, WPRT15, DPRT15 )

          CASE( 16 )
            CALL SET_TCODES( TDIM, NIPOL, MPRT16, WPRT16, DPRT16 )

          CASE DEFAULT

        END SELECT

    ENDDO                            ! End Loop on sorted x-ref entries

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

!******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

!.............  This internal subprogram stores the appropriate profile codes
!               to the arrays in the argument list.  It uses a different
!               approach for pol-specific or non-pol-specific.  It uses
!               variables from the main program.
    SUBROUTINE SET_TCODES( N, M, MPRT, WPRT, DPRT )

!.............  Subprogram arguments
        INTEGER   N            ! local size for dimensioning
        INTEGER   M            ! local size for dimensioning
        INTEGER   MPRT( N,M )  ! monthly codes
        INTEGER   WPRT( N,M )  ! weekly  codes
        INTEGER   DPRT( N,M )  ! diurnal codes

!......................................................................

        IF( ISP .EQ. 0 ) THEN
            MPRT( K,: ) = IMON
            WPRT( K,: ) = IWEK
            DPRT( K,: ) = IDIU

        ELSE
            MPRT( K,ISP ) = IMON
            WPRT( K,ISP ) = IWEK
            DPRT( K,ISP ) = IDIUPS

        ENDIF

    END SUBROUTINE SET_TCODES

END SUBROUTINE FILLTTBL

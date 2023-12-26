
SUBROUTINE FILLCTBL( NIPPA, NXREF, ICSIZE, XTYPE, XTCNT )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine populates the part of the grouped control cross-
!      reference tables that is the index to the control data tables
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!     Created 3/99 by M. Houyoux
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

!.........  MODULES for public variables
!.........  This module is for cross reference tables
    USE MODXREF, ONLY: ADDPS, INDXTA, ISPTA, MPRNA,&
    &        ICTL01, ICTL02, ICTL03, ICTL04, ICTL05,&
    &        ICTL06, ICTL07, ICTL08, ICTL09, ICTL10,&
    &        ICTL11, ICTL12, ICTL13, ICTL14, ICTL15, ICTL16,&
    &        ICTL02A, ICTL02B, ICTL02C,&
    &        ICTL05A, ICTL05B, ICTL05C,&
    &        ICTL08A, ICTL08B, ICTL08C,&
    &        ICTL26, ICTL27, ICTL28, ICTL29, ICTL30, ICTL31,&
    &        ICTL32, ICTL33, ICTL34, ICTL35, ICTL36, ICTL37, ICTL38

    IMPLICIT NONE

!...........   SUBROUTINE ARGUMENTS
    INTEGER, INTENT (IN) :: NIPPA           ! no. pollutants + activities
    INTEGER, INTENT (IN) :: NXREF           ! no. ungrpd x-ref entries
    INTEGER, INTENT (IN) :: ICSIZE( * )     ! size of x-ref groups
    INTEGER, INTENT (IN) :: XTYPE ( NXREF ) ! group no. of x-ref entry
    INTEGER, INTENT (IN) :: XTCNT ( NXREF ) ! pos. in x-ref group

!...........   Other local variables
    INTEGER       I, J, K, T     ! counter and indices
    INTEGER       IDX            ! tmp index to control data table
    INTEGER       IDXPS          ! IDX with ADDPS added for pol-specific
    INTEGER       ISP            ! tmp pollutant/activity position index
    INTEGER       TDIM           ! temporary table dimension

!***********************************************************************
!   begin body of subroutine FILLCTBL

!.........  Store the temporal profile codes for each x-ref entry, depending
!           on the group (XTYPE) and the position in that group (XTCNT)

    DO I = 1, NXREF

        J      = INDXTA( I )
        ISP    = ISPTA ( J )

        T      = XTYPE ( I )
        K      = XTCNT ( I )

!.............  Skip invalid or duplicate records
        IF( T .LE. 0 ) CYCLE

        TDIM   = ICSIZE( T )

        IDX    = MPRNA( J )
        IDXPS  = IDX + ADDPS

!.............  Populate tables depending on type. Note that the pollutant/
!               activity-specific entries are assumed to always come after the
!               non-specific ones (based on the previous sorting).
!.............  The pol-specific entries are stored by adding 90000
!               to the monthly profile number (which has a maximum of 3
!               digits) so that the pol-specific can be identified later
        SELECT CASE ( T )

          CASE( 1 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL01 )

          CASE( 2 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL02 )

          CASE( 3 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL03 )

          CASE( 4 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL04 )

          CASE( 5 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL05 )

          CASE( 6 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL06 )

          CASE( 7 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL07 )

          CASE( 8 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL08 )

          CASE( 9 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL09 )

          CASE( 10 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL10 )

          CASE( 11 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL11 )

          CASE( 12 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL12 )

          CASE( 13 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL13 )

          CASE( 14 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL14 )

          CASE( 15 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL15 )

          CASE( 16 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL16 )

          CASE( 17 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL02A )

          CASE( 18 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL02B )

          CASE( 19 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL02C )

          CASE( 20 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL05A )

          CASE( 21 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL05B )

          CASE( 22 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL05C )

          CASE( 23 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL08A )

          CASE( 24 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL08B )

          CASE( 25 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL08C )

          CASE( 26 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL26 )

          CASE( 27 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL27 )

          CASE( 28 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL28 )

          CASE( 29 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL29 )

          CASE( 30 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL30 )

          CASE( 31 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL31 )

          CASE( 32 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL32 )

          CASE( 33 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL33 )

          CASE( 34 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL34 )

          CASE( 35 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL35 )

          CASE( 36 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL36 )

          CASE( 37 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL37 )

          CASE( 38 )
            CALL SET_CNTRL_INDEX( TDIM, NIPPA, ICTL38 )

          CASE DEFAULT

        END SELECT

    END DO                            ! End Loop on sorted x-ref entries

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

!******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

!.............  This internal subprogram stores the appropriate index to the
!               control factors tables.  Note that the subprogram inherits
!               variables from the main program.
    SUBROUTINE SET_CNTRL_INDEX( N, M, ITABLE )

!.............  Subprogram arguments
        INTEGER   N             ! local size for dimensioning
        INTEGER   M             ! local size for dimensioning
        INTEGER   ITABLE( -1:N,M ) ! index to control data table

!......................................................................

        IF( ISP .EQ. 0 ) THEN
            ITABLE( K,: ) = IDX

        ELSE
            ITABLE( K,ISP ) = IDXPS

        ENDIF

    END SUBROUTINE SET_CNTRL_INDEX

END SUBROUTINE FILLCTBL

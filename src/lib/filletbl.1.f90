
SUBROUTINE FILLETBL( NACTV, NXREF, ICSIZE, XTYPE, XTCNT )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine populates the part of the grouped emission factor cross-
!      reference tables that is the index to the PSI table
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 6/99 by M. Houyoux
!
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!**************************************************************************
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

!.........  MODULES for public variables
!.........  This module is for cross reference tables
    USE MODXREF, ONLY: IEFS01, IEFS02, IEFS03, IEFS04, IEFS05,&
    &                   IEFS06, IEFS07, IEFS08, IEFS09, IEFS10,&
    &                   IEFS11, IEFS12, IEFS13, IEFS14, IEFS15,&
    &                   IEFS16, ADDPS, INDXTA, ISPTA

    IMPLICIT NONE

!...........   SUBROUTINE ARGUMENTS
    INTEGER, INTENT (IN) :: NACTV           ! no. activities
    INTEGER, INTENT (IN) :: NXREF           ! no. ungrpd x-ref entries
    INTEGER, INTENT (IN) :: ICSIZE( * )     ! size of x-ref groups
    INTEGER, INTENT (IN) :: XTYPE ( NXREF ) ! group no. of x-ref entry
    INTEGER, INTENT (IN) :: XTCNT ( NXREF ) ! pos. in x-ref group

!...........   Other local variables
    INTEGER       I, J, K, T     ! counter and indices
    INTEGER       IDX            ! tmp index to emission factor xref table
    INTEGER       JPS            ! IDX with ADDPS added for pol-specific
    INTEGER       ISP            ! tmp pollutant position index
    INTEGER       TDIM           ! temporary table dimension

    CHARACTER(300)     MESG      ! message buffer

    CHARACTER(16) :: PROGNAME = 'FILLETBL' ! program name

!***********************************************************************
!   begin body of subroutine FILLETBL

!.........  Check that the number of cross-reference entries does not
!           exceed the number of digits in the ADDPS factor
    I = ( ADDPS / 9 ) - 1
    IF( NXREF .GT. I ) THEN
        WRITE( MESG,94010 ) 'INTERNAL ERROR: ADDPS factor ' //&
        &       'must be increased higher than', ADDPS, 'for ' //&
        &       PROGNAME
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Store the emission factors unsorted cross-reference table index
!           for each x-ref entry, depending on the group (XTYPE) and the
!           position in that group (XTCNT)

    DO I = 1, NXREF

        J      = INDXTA( I )
        ISP    = ISPTA ( J )

        T      = XTYPE ( I )
        K      = XTCNT ( I )

!.............  Skip x-ref because it is invalid or duplicate
        IF( T .EQ. 0 ) CYCLE

        TDIM   = ICSIZE( T )

        JPS  = J + ADDPS

!.............  Populate tables depending on type. Note that the pollutant-
!               or activity-specific entries are assumed to always come after
!               the non-specific ones (based on the previous sorting).
!.............  The pol-specific entries are stored by adding 90000
!               to the index number (which has an expected upper limit of 4
!               digits) so that the pol-specific can be identified later
        SELECT CASE ( T )

          CASE( 1 )
            CALL SET_EFSPSI_INDEX( TDIM, NACTV, IEFS01 )

          CASE( 2 )
            CALL SET_EFSPSI_INDEX( TDIM, NACTV, IEFS02 )

          CASE( 3 )
            CALL SET_EFSPSI_INDEX( TDIM, NACTV, IEFS03 )

          CASE( 4 )
            CALL SET_EFSPSI_INDEX( TDIM, NACTV, IEFS04 )

          CASE( 5 )
            CALL SET_EFSPSI_INDEX( TDIM, NACTV, IEFS05 )

          CASE( 6 )
            CALL SET_EFSPSI_INDEX( TDIM, NACTV, IEFS06 )

          CASE( 7 )
            CALL SET_EFSPSI_INDEX( TDIM, NACTV, IEFS07 )

          CASE( 8 )
            CALL SET_EFSPSI_INDEX( TDIM, NACTV, IEFS08 )

          CASE( 9 )
            CALL SET_EFSPSI_INDEX( TDIM, NACTV, IEFS09 )

          CASE( 10 )
            CALL SET_EFSPSI_INDEX( TDIM, NACTV, IEFS10 )

          CASE( 11 )
            CALL SET_EFSPSI_INDEX( TDIM, NACTV, IEFS11 )

          CASE( 12 )
            CALL SET_EFSPSI_INDEX( TDIM, NACTV, IEFS12 )

          CASE( 13 )
            CALL SET_EFSPSI_INDEX( TDIM, NACTV, IEFS13 )

          CASE( 14 )
            CALL SET_EFSPSI_INDEX( TDIM, NACTV, IEFS14 )

          CASE( 15 )
            CALL SET_EFSPSI_INDEX( TDIM, NACTV, IEFS15 )

          CASE( 16 )
            CALL SET_EFSPSI_INDEX( TDIM, NACTV, IEFS16 )

          CASE DEFAULT

        END SELECT

    ENDDO                            ! End Loop on sorted x-ref entries

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

!******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

!.............  This internal subprogram stores the appropriate index to the
!               control factors tables.  Note that the subprogram inherits
!               variables from the main program.
    SUBROUTINE SET_EFSPSI_INDEX( N, M, ITABLE )

!.............  Subprogram arguments
        INTEGER   N             ! local size for dimensioning
        INTEGER   M             ! local size for dimensioning
        INTEGER   ITABLE( N,M ) ! index to control data table

!......................................................................

        IF( ISP .EQ. 0 ) THEN
            ITABLE( K,: ) = J

        ELSE
            ITABLE( K,ISP ) = JPS

        ENDIF

    END SUBROUTINE SET_EFSPSI_INDEX

END SUBROUTINE FILLETBL

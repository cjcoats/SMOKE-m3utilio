
SUBROUTINE FILLATBL( NXREF, ICSIZE, XTYPE, XTCNT )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine populates the area-to-point table that contains
    !      the table numbers, row numbers, and counts per FIPS code.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 11/02 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !****************************************************************************/
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

    !.......   This module is for cross reference tables
    USE MODXREF, ONLY: ARPT08, ARPT09, IARPTA, INDXTA

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN) :: NXREF               ! no. ungrpd x-ref entries
    INTEGER     , INTENT (IN) :: ICSIZE( * )         ! size of x-ref groups
    INTEGER     , INTENT (IN) :: XTYPE ( NXREF )     ! group no. of x-ref entry
    INTEGER     , INTENT (IN) :: XTCNT ( NXREF )     ! pos. in x-ref group

    !.......   Other local variables
    INTEGER       I, J, K, T           ! counter and indices
    INTEGER       NTBL                 ! tmp table number
    INTEGER       IROW                 ! tmp row number
    INTEGER       ICNT                 ! tmp counter of FIPS code

    LOGICAL    :: EFLAG = .FALSE.      ! true: error has occurred

    CHARACTER(257)         MESG        ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'FILLATBL'     ! program name

    !***********************************************************************
    !   begin body of subroutine FILLATBL

    !.....  Store the temporal profile codes for each x-ref entry, depending
    !           on the group (XTYPE) and the position in that group (XTCNT)

    DO I = 1, NXREF

        J    = INDXTA ( I )
        NTBL = IARPTA( J,1 )
        IROW = IARPTA( J,2 )
        ICNT = IARPTA( J,3 )

        T      = XTYPE ( I )
        K      = XTCNT ( I )

        !.....  Populate tables depending on type. Note that tables
        !                   are not pollutant-specific
        SELECT CASE ( T )

          CASE( 0 )          ! Skip this x-ref because it is invalid or duplicate

            !            CASE( 1 )
            !                ARPT01( 1 ) = NTBL
            !                ARPT01( 2 ) = IROW
            !                ARPT01( 3 ) = ICNT

            !            CASE( 2 )
            !                ARPT02( K,1 ) = NTBL
            !                ARPT02( K,2 ) = IROW
            !                ARPT02( K,3 ) = ICNT

            !            CASE( 3 )
            !                ARPT03( K,1 ) = NTBL
            !                ARPT03( K,2 ) = IROW
            !                ARPT03( K,3 ) = ICNT

            !            CASE( 4 )
            !                ARPT04( K,1 ) = NTBL
            !                ARPT04( K,2 ) = IROW
            !                ARPT04( K,3 ) = ICNT

            !            CASE( 5 )
            !                ARPT05( K,1 ) = NTBL
            !                ARPT05( K,2 ) = IROW
            !                ARPT05( K,3 ) = ICNT

            !            CASE( 6 )
            !                ARPT06( K,1 ) = NTBL
            !                ARPT06( K,2 ) = IROW
            !                ARPT06( K,3 ) = ICNT

            !            CASE( 7 )
            !                ARPT07( K,1 ) = NTBL
            !                ARPT07( K,2 ) = IROW
            !                ARPT07( K,3 ) = ICNT

          CASE( 8 )
            ARPT08( K,1 ) = NTBL
            ARPT08( K,2 ) = IROW
            ARPT08( K,3 ) = ICNT

          CASE( 9 )
            ARPT09( K,1 ) = NTBL
            ARPT09( K,2 ) = IROW
            ARPT09( K,3 ) = ICNT

          CASE DEFAULT

            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'INTERNAL ERROR: Group', T,     &
                   'not valid in subroutine ' // TRIM(PROGNAME)
            CALL M3MSG2( MESG )

        END SELECT

    END DO                                ! End loop on sorted x-ref entries

    IF( EFLAG ) THEN
        MESG = 'Problem processing area-to-point records.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats.... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE FILLATBL

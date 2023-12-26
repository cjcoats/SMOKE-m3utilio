
SUBROUTINE FILLMTBL( NXREF, ICSIZE, XTYPE, XTCNT )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine populates the mobile data part of the grouped VMT Mix
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 2/2000 by M. Houyoux
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
    USE MODXREF, ONLY: IMVS01, IMVS02, IMVS03, IMVS04, IMVS05,&
                       IMVS06, IMVS07, IMVS08, IMVS09, IMVS10,&
                       IMVS11, IMVS12, INDXTA

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
    INTEGER       IMVS                 ! tmp for gridding surrogate code

    LOGICAL    :: EFLAG = .FALSE.      ! true: error has occurred

    CHARACTER(300)         MESG        ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'FILLMTBL'     ! program name

    !***********************************************************************
    !   begin body of subroutine FILLMTBL

    !.....  Store the position for each unsorted x-ref entry, depending
    !           on the group (XTYPE) and the position in that group (XTCNT)

    DO I = 1, NXREF

        J    = INDXTA ( I )

        T      = XTYPE ( I )
        K      = XTCNT ( I )

        !.....  Populate tables depending on type. Note that tables
        !                   are not pollutant-specific
        SELECT CASE ( T )

          CASE( 0 )          ! Skip this x-ref because it is invalid or duplicate

          CASE( 1 )
            IMVS01 = J

          CASE( 2 )
            IMVS02( K ) = J

          CASE( 3 )
            IMVS03( K ) = J

          CASE( 4 )
            IMVS04( K ) = J

          CASE( 5 )
            IMVS05( K ) = J

          CASE( 6 )
            IMVS06( K ) = J

          CASE( 7 )
            IMVS07( K ) = J

          CASE( 8 )
            IMVS08( K ) = J

          CASE( 9 )
            IMVS09( K ) = J

          CASE( 10 )
            IMVS10( K ) = J

          CASE( 11 )
            IMVS11( K ) = J

          CASE( 12 )
            IMVS12( K ) = J

          CASE DEFAULT

            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'INTERNAL ERROR: Group', T,     &
                   'not valid in subroutine '// TRIM(PROGNAME)
            CALL M3MESG( MESG )

        END SELECT

    END DO                                ! End Loop on sorted x-ref entries

    IF( EFLAG ) THEN
        MESG = 'Problem processing cross-reference records.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats.... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE FILLMTBL

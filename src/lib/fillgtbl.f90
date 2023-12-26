
SUBROUTINE FILLGTBL( NXREF, ICSIZE, XTYPE, XTCNT )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine populates the gridding surrogate codes part of the
    !      grouped gridding cross-reference tables.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 4/99 by M. Houyoux
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
    USE MODXREF, ONLY: ISRG01, ISRG02, ISRG03, ISRG04, ISRG05,&
                       ISRG06, ISRG07, ISRG08, ISRG09,&
                       INDXTA, ISRGCDA

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
    INTEGER       ISRG                 ! tmp for gridding surrogate code

    LOGICAL    :: EFLAG = .FALSE.      ! true: error has occurred

    CHARACTER(300)         MESG        ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'FILLGTBL'     ! program name

    !***********************************************************************
    !   begin body of subroutine FILLGTBL

    !.....  Store the temporal profile codes for each x-ref entry, depending
    !           on the group (XTYPE) and the position in that group (XTCNT)

    DO I = 1, NXREF

        J    = INDXTA ( I )
        ISRG = ISRGCDA( J )

        T      = XTYPE ( I )
        K      = XTCNT ( I )

        !.....  Populate tables depending on type. Note that tables
        !                   are not pollutant-specific
        SELECT CASE ( T )

          CASE( 0 )      ! Skip this x-ref because it is invalid or duplicate

          CASE( 1 )
            ISRG01 = ISRG

          CASE( 2 )
            ISRG02( K ) = ISRG

          CASE( 3 )
            ISRG03( K ) = ISRG

          CASE( 4 )
            ISRG04( K ) = ISRG

          CASE( 5 )
            ISRG05( K ) = ISRG

          CASE( 6 )
            ISRG06( K ) = ISRG

          CASE( 7 )
            ISRG07( K ) = ISRG

          CASE( 8 )
            ISRG08( K ) = ISRG

          CASE( 9 )
            ISRG09( K ) = ISRG

          CASE DEFAULT

            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'INTERNAL ERROR: Group', T,     &
                   'not valid in subroutine ' // TRIM(PROGNAME)
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

END SUBROUTINE FILLGTBL

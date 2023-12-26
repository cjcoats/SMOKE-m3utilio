
SUBROUTINE FILLPTBL( NXREF, ICSIZE, XTYPE, XTCNT )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine populates the speed profile codes part of the
    !      grouped speed cross-reference tables.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 1/03 by C. Seppanen (copied from FILLGTBL)
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
    USE MODXREF, ONLY: ISPD01, ISPD02, ISPD03, ISPD04, ISPD05,  &
                       ISPD06, ISPD07, ISPD08, ISPD09,          &
                       INDXTA, ISPDCDA

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
    INTEGER       ISPD                 ! tmp for gridding surrogate code

    LOGICAL    :: EFLAG = .FALSE.      ! true: error has occurred

    CHARACTER(300)         MESG        ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'FILLPTBL'     ! program name

    !***********************************************************************
    !   begin body of subroutine FILLPTBL

    !.....  Store the speed profile code for each x-ref entry, depending
    !           on the group (XTYPE) and the position in that group (XTCNT)

    DO I = 1, NXREF

        J    = INDXTA ( I )
        ISPD = ISPDCDA( J )

        T      = XTYPE ( I )
        K      = XTCNT ( I )

        !.....  Populate tables depending on type. Note that tables
        !                   are not pollutant-specific
        SELECT CASE ( T )

          CASE( 0 )      ! Skip this x-ref because it is invalid or duplicate

          CASE( 1 )
            ISPD01 = ISPD

          CASE( 2 )
            ISPD02( K ) = ISPD

          CASE( 3 )
            ISPD03( K ) = ISPD

          CASE( 4 )
            ISPD04( K ) = ISPD

          CASE( 5 )
            ISPD05( K ) = ISPD

          CASE( 6 )
            ISPD06( K ) = ISPD

          CASE( 7 )
            ISPD07( K ) = ISPD

          CASE( 8 )
            ISPD08( K ) = ISPD

          CASE( 9 )
            ISPD09( K ) = ISPD

          CASE DEFAULT

            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'INTERNAL ERROR: Group', T, &
                   'not valid in subroutine', PROGNAME
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

END SUBROUTINE FILLPTBL

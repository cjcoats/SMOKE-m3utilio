
SUBROUTINE TMNAMUNT

    !***********************************************************************
    !  subroutine body starts at line 81
    !
    !  DESCRIPTION:
    !       This program creates the temporal emissions output file variable names
    !       and associated activities.  It also sets the units and conversion
    !       factors for creating the output emission values.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 10/99 by M. Houyoux
    !
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !*************************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !        System
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

    !.......  MODULES for public variables
    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: NIPPA, NIPOL, EANAM, EACNV, EAUNIT, EINAM

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters
    INCLUDE 'M6CNST3.h90'       !  Mobile6 constants

    !.......   EXTERNAL FUNCTIONS and their descriptions:
    CHARACTER(NAMLEN3), EXTERNAL :: MULTUNIT
    REAL              , EXTERNAL :: UNITFAC

    !.......   Other local variables
    INTEGER         I, J, K, L, L2, M         !  counters and indices
    INTEGER         IOS                   !  i/o status
    REAL            FAC1, FAC2            ! tmp conversion factors

    LOGICAL         EFLAG                 ! true: error found

    CHARACTER(NAMLEN3)  CURRUNIT       !  current unit
    CHARACTER(NAMLEN3)  CURRVNAME      !  current variable name
    CHARACTER(NAMLEN3)  CBUF           !  tmp variable name
    CHARACTER(300)      MESG           !  message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'TMNAMUNT'     ! program name

    !***********************************************************************
    !   begin body of subroutine TMNAMUNT
    EFLAG = .FALSE.
    !.......  Allocate memory for units conversions for inventory pollutants and
    !        activities (stored in MODINFO)
    ALLOCATE( EACNV( NIPPA ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EACNV', PROGNAME )

    EACNV  = 1.       ! array

    !.......  Now loop through pollutants and create units and conversion factors
    DO I = 1, NIPOL

        M = INDEX1( EINAM( I ), NIPPA, EANAM )

        CBUF = EAUNIT ( M )
        FAC1 = UNITFAC( CBUF, 'tons', .TRUE. )
        FAC2 = UNITFAC( EAUNIT( M ), '1/yr', .FALSE. )

        IF ( FAC1 .LT. 0. ) FAC1 = 1.
        IF ( FAC2 .LT. 0. ) FAC2 = 1.

        !.......  keep the orig miles/hr for Movesmrg to process MOVES lookup tables.
        IF( INDEX( CBUF,'miles' ) > 0 ) THEN
            EAUNIT( M ) = 'miles/hr'
        ELSE
            EAUNIT( M ) = 'tons/hr'
        END IF

        EACNV ( M ) = FAC1 / FAC2

    END DO

    !.......  Abort if error was found
    IF ( EFLAG ) THEN
        MESG = 'Problem with emission types or emission factors'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    RETURN

END SUBROUTINE TMNAMUNT

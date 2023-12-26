
SUBROUTINE ALOCCHRT( ICSIZE )

    !***********************************************************************
    !  subroutine body starts at line 60
    !
    !  DESCRIPTION:
    !      This subroutine allocates memory for the portion of the grouped
    !      cross-reference tables that contain the source characteristics.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 3/99 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
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

    !.......   This module is for cross reference tables
    USE MODXREF, ONLY: CHRT02, CHRT03, CHRT04, CHRT05, CHRT06,  &
        CHRT07, CHRT08, CHRT09, CHRT10, CHRT11,                 &
        CHRT12, CHRT13, CHRT14, CHRT15, CHRT16,                 &
        CHRT02A, CHRT02B, CHRT02C,                              &
        CHRT05A, CHRT05B, CHRT05C,                              &
        CHRT08A, CHRT08B, CHRT08C,                              &
        CHRT26, CHRT27, CHRT28, CHRT29, CHRT30, CHRT31,         &
        CHRT32, CHRT33, CHRT34, CHRT35, CHRT36, CHRT37, CHRT38

    IMPLICIT NONE

    !.......   SUBROUTINE ARGUMENTS
    INTEGER, INTENT(IN) :: ICSIZE( * )      ! size of x-ref groups

    !.......   Other local variables

    INTEGER       IOS                  ! i/o status

    CHARACTER(16), PARAMETER :: PROGNAME = 'ALOCCHRT'     ! program name

    !***********************************************************************
    !   begin body of subroutine ALOCCHRT

    !.....  First deallocate if these have previously been allocated
    IF ( ALLOCATED( CHRT02 ) ) THEN

        DEALLOCATE( CHRT02, CHRT03, CHRT04, CHRT05, CHRT06 )
        DEALLOCATE( CHRT07, CHRT08, CHRT09, CHRT10, CHRT11 )
        DEALLOCATE( CHRT12, CHRT13, CHRT14, CHRT15, CHRT16 )
        DEALLOCATE( CHRT02A, CHRT02B, CHRT02C )
        DEALLOCATE( CHRT05A, CHRT05B, CHRT05C )
        DEALLOCATE( CHRT08A, CHRT08B, CHRT08C )
        DEALLOCATE( CHRT26, CHRT27, CHRT28, CHRT29, CHRT30, CHRT31 )
        DEALLOCATE( CHRT32, CHRT33, CHRT34, CHRT35, CHRT36, CHRT37 )
        DEALLOCATE( CHRT38 )

    END IF

    ALLOCATE( CHRT02( MAX( 1, ICSIZE(  2 ) ) ),        &
              CHRT03( MAX( 1, ICSIZE(  3 ) ) ),        &
              CHRT04( MAX( 1, ICSIZE(  4 ) ) ),        &
              CHRT05( MAX( 1, ICSIZE(  5 ) ) ),        &
              CHRT06( MAX( 1, ICSIZE(  6 ) ) ),        &
              CHRT07( MAX( 1, ICSIZE(  7 ) ) ),        &
              CHRT08( MAX( 1, ICSIZE(  8 ) ) ),        &
              CHRT09( MAX( 1, ICSIZE(  9 ) ) ),        &
              CHRT10( MAX( 1, ICSIZE( 10 ) ) ),        &
              CHRT11( MAX( 1, ICSIZE( 11 ) ) ),        &
              CHRT12( MAX( 1, ICSIZE( 12 ) ) ),        &
              CHRT13( MAX( 1, ICSIZE( 13 ) ) ),        &
              CHRT14( MAX( 1, ICSIZE( 14 ) ) ),        &
              CHRT15( MAX( 1, ICSIZE( 15 ) ) ),        &
              CHRT16( MAX( 1, ICSIZE( 16 ) ) ),        &
             CHRT02A( MAX( 1, ICSIZE( 17 ) ) ),        &
             CHRT02B( MAX( 1, ICSIZE( 18 ) ) ),        &
             CHRT02C( MAX( 1, ICSIZE( 19 ) ) ),        &
             CHRT05A( MAX( 1, ICSIZE( 20 ) ) ),        &
             CHRT05B( MAX( 1, ICSIZE( 21 ) ) ),        &
             CHRT05C( MAX( 1, ICSIZE( 22 ) ) ),        &
             CHRT08A( MAX( 1, ICSIZE( 23 ) ) ),        &
             CHRT08B( MAX( 1, ICSIZE( 24 ) ) ),        &
             CHRT08C( MAX( 1, ICSIZE( 25 ) ) ),        &
              CHRT26( MAX( 1, ICSIZE( 26 ) ) ),        &
              CHRT27( MAX( 1, ICSIZE( 27 ) ) ),        &
              CHRT28( MAX( 1, ICSIZE( 28 ) ) ),        &
              CHRT29( MAX( 1, ICSIZE( 29 ) ) ),        &
              CHRT30( MAX( 1, ICSIZE( 30 ) ) ),        &
              CHRT31( MAX( 1, ICSIZE( 31 ) ) ),        &
              CHRT32( MAX( 1, ICSIZE( 32 ) ) ),        &
              CHRT33( MAX( 1, ICSIZE( 33 ) ) ),        &
              CHRT34( MAX( 1, ICSIZE( 34 ) ) ),        &
              CHRT35( MAX( 1, ICSIZE( 35 ) ) ),        &
              CHRT36( MAX( 1, ICSIZE( 36 ) ) ),        &
              CHRT37( MAX( 1, ICSIZE( 37 ) ) ),        &
              CHRT38( MAX( 1, ICSIZE( 38 ) ) ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CHRT02...CHRT38', PROGNAME )

    RETURN

END SUBROUTINE ALOCCHRT


SUBROUTINE ALOCETBL( NACTV, ICSIZE )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine allocates memory for the portion of the emission factor
    !      cross-reference tables that contain the indices to the unsorted PSIs, and
    !      it initializes these to missing.  The subroutine arguments are the number
    !      of inventory pollutants and an array that contains the dimensions for
    !      each of the different groups of the cross-reference.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 6/99 by M. Houyoux
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
    USE MODXREF, ONLY: IEFS01, IEFS02, IEFS03, IEFS04, IEFS05,  &
                       IEFS06, IEFS07, IEFS08, IEFS09, IEFS10,  &
                       IEFS11, IEFS12, IEFS13, IEFS14, IEFS15,  &
                       IEFS16

    IMPLICIT NONE

    !.......   SUBROUTINE ARGUMENTS
    INTEGER, INTENT(IN) :: NACTV            ! number of pollutants
    INTEGER, INTENT(IN) :: ICSIZE( * )      ! size of x-ref groups

    !.......   Other local variables
    INTEGER       IOS       ! i/o status

    CHARACTER(16), PARAMETER :: PROGNAME = 'ALOCETBL'     ! program name

    !***********************************************************************
    !   begin body of subroutine ALOCETBL

    !.....  First deallocate if these have previously been allocated
    IF ( ALLOCATED( IEFS02 ) ) THEN

        DEALLOCATE( IEFS02, IEFS03, IEFS04, IEFS05, IEFS06 )
        DEALLOCATE( IEFS07, IEFS08, IEFS09, IEFS10, IEFS11 )
        DEALLOCATE( IEFS12, IEFS13, IEFS14, IEFS15, IEFS16 )

    END IF

    ALLOCATE( IEFS01(              NACTV ),     &
              IEFS02( ICSIZE(  2 ),NACTV ),     &
              IEFS03( ICSIZE(  3 ),NACTV ),     &
              IEFS04( ICSIZE(  4 ),NACTV ),     &
              IEFS05( ICSIZE(  5 ),NACTV ),     &
              IEFS06( ICSIZE(  6 ),NACTV ),     &
              IEFS07( ICSIZE(  7 ),NACTV ),     &
              IEFS08( ICSIZE(  8 ),NACTV ),     &
              IEFS09( ICSIZE(  9 ),NACTV ),     &
              IEFS10( ICSIZE( 10 ),NACTV ),     &
              IEFS11( ICSIZE( 11 ),NACTV ),     &
              IEFS12( ICSIZE( 12 ),NACTV ),     &
              IEFS13( ICSIZE( 13 ),NACTV ),     &
              IEFS14( ICSIZE( 14 ),NACTV ),     &
              IEFS15( ICSIZE( 15 ),NACTV ),     &
              IEFS16( ICSIZE( 16 ),NACTV ), STAT=IOS )
    CALL CHECKMEM( IOS, 'IEFS01IEFS16', PROGNAME )

    IEFS01 = IMISS3
    IEFS02 = IMISS3
    IEFS03 = IMISS3
    IEFS04 = IMISS3
    IEFS05 = IMISS3
    IEFS06 = IMISS3
    IEFS07 = IMISS3
    IEFS08 = IMISS3
    IEFS09 = IMISS3
    IEFS10 = IMISS3
    IEFS11 = IMISS3
    IEFS12 = IMISS3
    IEFS13 = IMISS3
    IEFS14 = IMISS3
    IEFS15 = IMISS3
    IEFS16 = IMISS3

    RETURN

END SUBROUTINE ALOCETBL

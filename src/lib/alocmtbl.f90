
SUBROUTINE ALOCMTBL( ICSIZE )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine allocates memory for the portion of the VMT mix
    !      table that contain the various mobile-source characteristics. The
    !      subroutine argument is an array that contains the dimensions for each
    !      of the different groups of the cross-reference.
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
    USE MODXREF, ONLY: IMVS01, IMVS02, IMVS03, IMVS04, IMVS05,  &
                       IMVS06, IMVS07, IMVS08, IMVS09, IMVS10,  &
                       IMVS11, IMVS12

    IMPLICIT NONE

    !.......   SUBROUTINE ARGUMENTS
    INTEGER, INTENT(IN) :: ICSIZE( * )      ! size of x-ref groups

    !.......   Other local variables
    INTEGER       IOS                  ! i/o status

    CHARACTER(16), PARAMETER :: PROGNAME = 'ALOCMTBL'     ! program name

    !***********************************************************************
    !   begin body of subroutine ALOCMTBL

    !.....  First deallocate if these have previously been allocated
    IF ( ALLOCATED( IMVS02 ) ) THEN

        DEALLOCATE( IMVS02, IMVS03, IMVS04, IMVS05 )
        DEALLOCATE( IMVS06, IMVS07, IMVS08, IMVS09 )

    END IF

    ALLOCATE( IMVS02( ICSIZE( 2 ) ),        &
              IMVS03( ICSIZE( 3 ) ),        &
              IMVS04( ICSIZE( 4 ) ),        &
              IMVS05( ICSIZE( 5 ) ),        &
              IMVS06( ICSIZE( 6 ) ),        &
              IMVS07( ICSIZE( 7 ) ),        &
              IMVS08( ICSIZE( 8 ) ),        &
              IMVS09( ICSIZE( 9 ) ),        &
              IMVS10( ICSIZE(10 ) ),        &
              IMVS11( ICSIZE(11 ) ),        &
              IMVS12( ICSIZE(12 ) ), STAT=IOS )
    CALL CHECKMEM( IOS, 'IMVS02...IMVS12', PROGNAME )

    IMVS02 = IMISS3        ! array
    IMVS03 = IMISS3        ! array
    IMVS04 = IMISS3        ! array
    IMVS05 = IMISS3        ! array
    IMVS06 = IMISS3        ! array
    IMVS07 = IMISS3        ! array
    IMVS08 = IMISS3        ! array
    IMVS09 = IMISS3        ! array
    IMVS09 = IMISS3        ! array
    IMVS09 = IMISS3        ! array
    IMVS09 = IMISS3        ! array

    RETURN

END SUBROUTINE ALOCMTBL

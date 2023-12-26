
SUBROUTINE ALOCPTBL( ICSIZE )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine allocates memory for the portion of the speed
    !      cross-reference tables that contain the speed profile numbers, and
    !      it initializes these to missing.  The subroutine argument is
    !      an array that contains the dimensions for each of the different
    !      groups of the cross-reference.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 1/03 by C. Seppanen
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
                       ISPD06, ISPD07, ISPD08, ISPD09

    IMPLICIT NONE

    !.......   SUBROUTINE ARGUMENTS
    INTEGER, INTENT(IN) :: ICSIZE( * )      ! size of x-ref groups

    !.......   Other local variables
    INTEGER       IOS                  ! i/o status

    CHARACTER(16), PARAMETER :: PROGNAME = 'ALOCPTBL'     ! program name

    !***********************************************************************
    !   begin body of subroutine ALOCPTBL

    !.....  First deallocate if these have previously been allocated
    IF( ALLOCATED( ISPD02 ) ) THEN

        DEALLOCATE( ISPD02, ISPD03, ISPD04, ISPD05 )
        DEALLOCATE( ISPD06, ISPD07, ISPD08, ISPD09 )

    END IF
        ! SCC=left, FIP=0
    ALLOCATE( ISPD02( ICSIZE( 2 ) ),    &
              ISPD03( ICSIZE( 3 ) ),    &
              ISPD04( ICSIZE( 4 ) ),    &
              ISPD05( ICSIZE( 5 ) ),    &
              ISPD06( ICSIZE( 6 ) ),    &
              ISPD07( ICSIZE( 7 ) ),    &
              ISPD08( ICSIZE( 8 ) ),    &
              ISPD09( ICSIZE( 9 ) ), STAT=IOS )
    CALL CHECKMEM( IOS, 'ISPD09', PROGNAME )

    ISPD01 = IMISS3         ! scalar
    ISPD02 = IMISS3         ! array
    ISPD03 = IMISS3         ! array
    ISPD04 = IMISS3         ! array
    ISPD05 = IMISS3         ! array
    ISPD06 = IMISS3         ! array
    ISPD07 = IMISS3         ! array
    ISPD08 = IMISS3         ! array
    ISPD09 = IMISS3         ! array

    RETURN

END SUBROUTINE ALOCPTBL

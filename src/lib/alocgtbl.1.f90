
SUBROUTINE ALOCGTBL( ICSIZE )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine allocates memory for the portion of the gridding
!      cross-reference tables that contain the gridding surrogate code numbers,
!      and it initializes these to missing.  The subroutine argument is
!      an array that contains the dimensions for each of the different groups
!      of the cross-reference.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!     Created 3/99 by M. Houyoux
!
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

!...........   This module is for cross reference tables
    USE MODXREF, ONLY: ISRG01, ISRG02, ISRG03, ISRG04, ISRG05,&
    &                   ISRG06, ISRG07, ISRG08, ISRG09

    IMPLICIT NONE

!...........   SUBROUTINE ARGUMENTS
    INTEGER, INTENT(IN) :: ICSIZE( * )  ! size of x-ref groups

!...........   Other local variables
    INTEGER       IOS              ! i/o status

    CHARACTER(16), PARAMETER :: PROGNAME = 'ALOCGTBL' ! program name

!***********************************************************************
!   begin body of subroutine ALOCGTBL

!.........  First deallocate if these have previously been allocated
    IF ( ALLOCATED( ISRG02 ) ) THEN

        DEALLOCATE( ISRG02, ISRG03, ISRG04, ISRG05,&
        &            ISRG06, ISRG07, ISRG08, ISRG09 )

    END IF

    ALLOCATE( ISRG02( ICSIZE( 2 ) ),&
    &          ISRG03( ICSIZE( 3 ) ),&
    &          ISRG04( ICSIZE( 4 ) ),&
    &          ISRG05( ICSIZE( 5 ) ),&
    &          ISRG06( ICSIZE( 6 ) ),&
    &          ISRG07( ICSIZE( 7 ) ),&
    &          ISRG08( ICSIZE( 8 ) ),&
    &          ISRG09( ICSIZE( 9 ) ), STAT=IOS )
    CALL CHECKMEM( IOS, 'ISRG02...ISRG09', PROGNAME )

    ISRG01 = IMISS3     ! is scalar
    ISRG02 = IMISS3     ! array
    ISRG03 = IMISS3     ! array
    ISRG04 = IMISS3     ! array
    ISRG05 = IMISS3     ! array
    ISRG06 = IMISS3     ! array
    ISRG07 = IMISS3     ! array
    ISRG08 = IMISS3     ! array
    ISRG09 = IMISS3     ! array

    RETURN

END SUBROUTINE ALOCGTBL

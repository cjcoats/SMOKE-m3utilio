
SUBROUTINE ALOCATBL( ICSIZE )

!***********************************************************************
!  subroutine body starts at line 46
!
!  DESCRIPTION:
!      This subroutine allocates memory for the portion of the area-
!      to-point file that contains the table number, row number, and count,
!      and it initializes these to missing.  The subroutine argument is
!      an array that contains the dimensions for each of the different groups
!      of the cross-reference. Only group 9 (full FIPS, full SCC) is
!      currently valid for the area-to-point assignments.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 11/02 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
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
    USE MODXREF, ONLY: ARPT08, ARPT09

    IMPLICIT NONE

!...........   SUBROUTINE ARGUMENTS
    INTEGER, INTENT(IN) :: ICSIZE( * )  ! size of x-ref groups

!...........   Other local variables
    INTEGER       J, K     ! counter and indices
    INTEGER       IOS              ! i/o status

    CHARACTER(16), PARAMETER :: PROGNAME = 'ALOCATBL' ! program name

!***********************************************************************
!   begin body of subroutine ALOCATBL

!.........  First deallocate if these have previously been allocated
    IF ( ALLOCATED( ARPT08 ) )  DEALLOCATE( ARPT08, ARPT09 )

    J = ICSIZE( 8 )
    K = ICSIZE( 9 )                                ! SCC=left, FIP=all
    ALLOCATE( ARPT08( J,3 ),&
    &          ARPT09( K,3 ), STAT=IOS )
    CALL CHECKMEM( IOS, 'ARPT08,ARPT09', PROGNAME )
    ARPT09 = IMISS3    ! array

    RETURN

END SUBROUTINE ALOCATBL


SUBROUTINE PADNZERO( N, STRING )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine right-justifies a string and replaces the leading spaces
!      with zeros for a string of size N.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!     Subroutines: Models-3 subroutines
!     Functions: Models-3 functions
!
!  REVISION  HISTORY:
!     Created 3/02 by G. Cano
!
!****************************************************************************/
!
! Project Title: EDSS Tools Library
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
! Pathname: @(#)$Id$
! Last updated: $Date$
!
!***************************************************************************

    IMPLICIT NONE

!...........   SUBROUTINE ARGUMENTS
    INTEGER,      INTENT(IN)     :: N      ! length of the string to pad
    CHARACTER(N), INTENT(IN OUT) :: STRING ! character string to adjust

    INTEGER         I, L

    CHARACTER(16) :: PROGNAME = 'PADNZERO' ! program name

!***********************************************************************
!   begin body of subroutine PADZERO

    STRING = ADJUSTR( STRING )

    DO I = 1, N

        IF( STRING( I:I ) .EQ. ' ' ) THEN

            STRING( I:I ) = '0'

        ELSE

            EXIT

        ENDIF

    ENDDO

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE PADNZERO

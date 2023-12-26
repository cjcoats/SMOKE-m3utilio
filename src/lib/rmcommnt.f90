
SUBROUTINE RMCOMMNT( CMT_DELIM, LINE )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!    The RMCOMMNT routine will delete the comment part of a line, if it is
!    there.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!
!***********************************************************************
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
!***********************************************************************

    IMPLICIT NONE

!...........   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN)     :: CMT_DELIM   ! comment delimeter
    CHARACTER(*), INTENT (IN OUT) :: LINE        ! line of data

!...........   Local variables
    INTEGER      J

!***********************************************************************
!   begin body of subroutine RMCOMMNT

    J = INDEX( LINE, CMT_DELIM )

!.........  If comment delimeter is found, then remove remainder of line
!           MRH update 1/2/09: remove the comment field too and allow for
!           length of comment delimiter greater than 1 character (e.g.,
!           Smkreport's ## in-line comment delimiter.
    IF( J .GT. 0 ) THEN

        J = J - 1
        LINE = LINE( 1:J )

    END IF

    RETURN

END SUBROUTINE RMCOMMNT



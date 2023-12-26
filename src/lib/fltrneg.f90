
SUBROUTINE FLTRNEG( STRING )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine filters the valid "missing" terms that can be in a
    !      cross-reference file changes the argument value to blank when it is
    !      one of the missing values.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !
    !  REVISION  HISTORY:
    !       Started 6/99 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
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
    ! Pathname: $Source$
    ! Last updated: $Date$
    !
    !***************************************************************************

    IMPLICIT NONE

    !.....  Subroutine arguments
    CHARACTER(*), INTENT (IN OUT) :: STRING       ! string for filter '-9'

    !***********************************************************************
    !   Begin body of subroutine FLTRNEG

    IF( INDEX( STRING, '-9' ) .GT.  0  .OR.     &
        STRING                .EQ. '0'      ) STRING = ' '

    RETURN

END SUBROUTINE FLTRNEG

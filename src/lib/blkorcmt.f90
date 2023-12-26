
LOGICAL FUNCTION BLKORCMT( LINE )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !    The BLKORCMT routine will return true of the buffer given it is a blank
    !    or a comment.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version ??/???? by ???
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
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

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: LINE            ! line of data

    !.......   Local variables
    INTEGER      J

    !***********************************************************************
    !   begin body of subroutine BLKORCMT


    IF( LINE .EQ. ' ' .OR. LINE( 1:1 ) .EQ. CINVHDR ) THEN
        BLKORCMT = .TRUE.

    ELSE
        BLKORCMT = .FALSE.

    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats.... 94xxx

94010 FORMAT( 10( A, :, I10, :, 1X ) )

END FUNCTION BLKORCMT



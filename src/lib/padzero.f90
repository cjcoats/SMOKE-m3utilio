
SUBROUTINE PADZERO( STRING )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine right-justifies a string and replaces the leading spaces
    !      with zeros.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !     Subroutines: Models-3 subroutines
    !     Functions: Models-3 functions
    !
    !  REVISION  HISTORY:
    !     Created 1/99 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, ".f90" source format, and
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
    ! Pathname: @(#)$Id$
    ! Last updated: $Date$
    !
    !***************************************************************************

    IMPLICIT NONE

    !......   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT(IN OUT) :: STRING     ! character string to adjust

    INTEGER         I, L

    CHARACTER(16) :: PROGNAME = 'PADZERO'     ! program name

    !***********************************************************************
    !   begin body of subroutine PADZERO

    STRING = ADJUSTR( STRING )

    L = LEN( STRING )

    DO I = 1, L

        IF( STRING( I:I ) .EQ. ' ' ) THEN

            STRING( I:I ) = '0'

        ELSE

            EXIT

        ENDIF

    ENDDO

    RETURN

END  SUBROUTINE PADZERO

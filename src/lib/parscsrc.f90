
SUBROUTINE PARSCSRC( INSTRING, MAXN, STARTS, ENDS, OUTCOL, NARR, STRARRAY )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine splits a CSOURC entry into a character
    !      array and returns only non-blank strings and the array length
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version ??/???? by ???
    !       Version 11/2023 by CJC:  USE M3UTILIO, ".f90" source format, and
    !       related changes
    !**************************************************************************
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

    IMPLICIT NONE

    !......   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !......   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: INSTRING           ! input string
    INTEGER     , INTENT (IN) :: MAXN               ! max number of fields
    INTEGER     , INTENT (IN) :: STARTS( MAXN )     ! starting field pos'ns
    INTEGER     , INTENT (IN) :: ENDS  ( MAXN )     ! ending field pos'ns
    LOGICAL     , INTENT (IN) :: OUTCOL( MAXN )     ! true if column is valid
    INTEGER     , INTENT(OUT) :: NARR               ! no. of non-blank chars
    CHARACTER(*), INTENT(OUT) :: STRARRAY( * )      ! output array of strings

    !......   Other local variables
    INTEGER         I, J, L1, L2

    CHARACTER(300)  BUFFER

    !***********************************************************************
    !   begin body of subroutine PARSCSRC

    J = 0
    DO I = 1, MAXN
        IF( OUTCOL( I ) ) THEN
            J = J + 1

            !......  Retrieve and left-justify contents of this part of string
            L1 = STARTS( I )
            L2 = ENDS  ( I )
            BUFFER = ADJUSTL( INSTRING( L1:L2 ) )

            !......  Convert missing entries to blanks
            IF( BUFFER .EQ. EMCMISS3 ) BUFFER = ' '

            STRARRAY( J ) = BUFFER

        ENDIF
    ENDDO

    NARR = J

    RETURN

END SUBROUTINE PARSCSRC


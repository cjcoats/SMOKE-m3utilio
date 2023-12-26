
CHARACTER(16) FUNCTION VERCHAR( BUFFER )

    !***********************************************************************
    !  function body starts at line
    !
    !  DESCRIPTION:
    !      This function returns the version number from a string filled in
    !      by SCCS operations
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version 11/2023 by CJC
    !
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

    !...........   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT( IN ) :: BUFFER         ! string FIPS code

    !...........   Other local variables
    INTEGER         J

    !***********************************************************************
    !   begin body of function VERCHAR

    VERCHAR = ADJUSTL( BUFFER )

    IF( VERCHAR( 1:1 ) .EQ. '%' ) THEN
        VERCHAR = 'Dev'
    ELSE
        J  = INDEX  ( VERCHAR, ' ' )
        VERCHAR = ADJUSTL( VERCHAR( J+1: ) )
    ENDIF

    RETURN

END FUNCTION VERCHAR


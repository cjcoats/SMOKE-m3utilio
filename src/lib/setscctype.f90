
LOGICAL FUNCTION SETSCCTYPE( TSCC )

    !***********************************************************************
    !  function body starts at line 68
    !
    !  DESCRIPTION:
    !       Checks SCC code and resets parameters based on type
    !
    !  PRECONDITIONS REQUIRED:
    !       CATEGORY type must be set in MODINFO
    !       SCC must be 10-digits long and right-justified
    !       8-digit SCCs must start with '00'
    !
    !  SUBROUTINES AND FUNCTIONS CALLED: none
    !
    !  REVISION  HISTORY:
    !     7/03: Created by C. Seppanen
    !       12/2023 by CJC:  ".f90" code format changes
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

    !.........  MODULES for public variables
    !.........  This module contains the information about the source category
    USE MODINFO, ONLY: LSCCEND, RSCCBEG, SCCLEV1, SCCLEV2, SCCLEV3, CATEGORY

    IMPLICIT NONE

    !...........   INCLUDES:
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !........  Function arguments
    CHARACTER(*), INTENT (IN) :: TSCC       ! SCC code

    !***********************************************************************
    !   begin body of function SETSCCTYPE

    SETSCCTYPE = .FALSE.

    !.........  Don't change any parameters if category is mobile
    IF( CATEGORY == 'MOBILE' ) RETURN

    !.........  Check if first two digits of SCC are zero
    IF( TSCC( SCCEXPLEN3+1:SCCEXPLEN3+2 ) == '00' ) THEN

        !.............  Only set new values if needed and set flag
        IF( LSCCEND /= SCCEXPLEN3 + 5 ) THEN
            SETSCCTYPE = .TRUE.          ! flag indicates that values have been changed
            LSCCEND = SCCEXPLEN3 + 5
            RSCCBEG = SCCEXPLEN3 + 6
            SCCLEV1 = SCCEXPLEN3 + 3
            SCCLEV2 = SCCEXPLEN3 + 5
            SCCLEV3 = SCCEXPLEN3 + 8
        END IF
    ELSE
        IF( LSCCEND /= SCCEXPLEN3 + 7 ) THEN
            SETSCCTYPE = .TRUE.
            LSCCEND = SCCEXPLEN3 + 7
            RSCCBEG = SCCEXPLEN3 + 8
            SCCLEV1 = SCCEXPLEN3 + 2
            SCCLEV2 = SCCEXPLEN3 + 4
            SCCLEV3 = SCCEXPLEN3 + 7
        END IF
    END IF

    RETURN

END FUNCTION SETSCCTYPE

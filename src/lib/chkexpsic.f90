
LOGICAL FUNCTION CHKEXPSIC( CSIC )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      Returns true if SIC is expanded.
    !
    !  REVISION  HISTORY:
    !       Created 10/13 by C. Seppanen
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !****************************************************************************/
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !                System
    ! File: @(#)$Id$
    !
    ! COPYRIGHT (C) 2005, Environmental Modeling for Policy Development
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

    !.......   INCLUDES:
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......  EXTERNAL FUNCTIONS and their descriptions:
    LOGICAL, EXTERNAL :: CHKINT

    !....  Function arguments
    CHARACTER(SICLEN3), INTENT(IN) :: CSIC       ! SIC code

    !***********************************************************************
    !   begin body of function CHKEXPSIC

    CHKEXPSIC = .FALSE.

    !.....  Check if any of the expanded characters are not zero
    IF( CSIC( 1:SICEXPLEN3 ) /= REPEAT( '0', SICEXPLEN3 ) ) THEN
        CHKEXPSIC = .TRUE.
        RETURN
    END IF

    !.....  Check if last 4 characters are not digits
    IF( .NOT. CHKINT( CSIC( SICLEN3-3:SICLEN3 ) ) ) THEN
        CHKEXPSIC = .TRUE.
        RETURN
    END IF

    RETURN

END FUNCTION CHKEXPSIC


LOGICAL FUNCTION CHKEXPSCC( TSCC )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      Returns true if SCC is expanded.
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

    !....  Function arguments
    CHARACTER(SCCLEN3), INTENT(IN) :: TSCC       ! SCC code

    !***********************************************************************
    !   begin body of function CHKEXPSCC

    CHKEXPSCC = .FALSE.

    IF( TSCC( 1:SCCEXPLEN3 ) /= '0000000000' ) THEN
        CHKEXPSCC = .TRUE.
    END IF

    RETURN

END FUNCTION CHKEXPSCC

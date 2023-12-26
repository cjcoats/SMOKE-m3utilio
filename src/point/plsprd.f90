
SUBROUTINE PLSPRD( DTHDZ, ZF, KZ, CEFSTK, HTMIX, PLTOP, PLBOT )

    !***********************************************************************
    !  subroutine body starts at line 70
    !
    !  DESCRIPTION:
    !       Calculates the initial vertical spread of a plume; modified
    !       from Gillani's model.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !
    !  REVISION  HISTORY:
    !       Created from code by Jim Godowitch at EPA, 9/03
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

    !.......   ARGUMENTS and their descriptions:

    REAL,    INTENT  (IN) :: DTHDZ( KZ )        ! potential temperature lapse rate (K/m)
    REAL,    INTENT  (IN) :: ZF( 0:KZ )         ! full-layer heights (m)
    INTEGER, INTENT  (IN) :: KZ                 ! number of emissions layers
    REAL,    INTENT  (IN) :: CEFSTK             ! effective stack height (m)
    REAL,    INTENT  (IN) :: HTMIX              ! mixing height (m)
    REAL,    INTENT (OUT) :: PLTOP              ! plume top (m)
    REAL,    INTENT (OUT) :: PLBOT              ! plume bottom (m)

    !.......   PARAMETERS and their descriptions:
    REAL, PARAMETER :: SZ0FAC = 3.545        ! factor used to derive plume depth
    REAL, PARAMETER :: SPRFAC = 15.          ! empirical coefficient for vertical spread
    REAL, PARAMETER :: GAMA   = -0.0098      ! adiabatic lapse rate (K/m)

    !.......   Local variables
    INTEGER    K

    REAL       SIGZ0
    REAL       DTDZ
    REAL       DPTH

    !***********************************************************************
    !   begin body of subroutine  PLSPRD

    !......  Get ambient temperature above plume rise height (effective stack height)
    K = 0
    DO
        K = K + 1
        IF( K == KZ .OR. CEFSTK <= ZF( K ) ) EXIT
    END DO
    DTDZ  = DTHDZ( K ) + GAMA

    !......  Compute initial vertical spread
    SIGZ0 = MAX( 10.0, SPRFAC * EXP( -117. * DTDZ ) )
    DPTH  = SZ0FAC * SIGZ0

    !......  Compute plume top and bottom heights; plume is either completely
    !          within or outside mixing layer
    PLTOP = CEFSTK + DPTH/2.
    PLBOT = CEFSTK - DPTH/2.

    !......  Make sure plume bottom is at least zero
    PLBOT = MAX( 0.0, PLBOT )

    !......  Make sure that plume top and bottom heights are less than
    !          the top layer's top and bottom heights
    PLTOP = MIN( ZF( KZ ), PLTOP )
    PLBOT = MIN( ZF( KZ ) - 1., PLBOT )

    RETURN

END SUBROUTINE PLSPRD


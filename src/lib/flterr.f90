LOGICAL FUNCTION FLTERR( P, Q )

    !***********************************************************************
    !  function body starts at line 52
    !
    !  DESCRIPTION:
    !      These functions determine if two REAL or REAL*8 numbers are
    !      "definitely unequal", assuming that they may come from different
    !      computer syste,s (including IBM mainframe and WMO GRIB)
    !      It compares the square of the
    !      normalized difference against the square of the tolerance and
    !      returns TRUE if and only if the numbers are significantly
    !      different.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created by C. Seppanen 6/04 from FLTERR.h90
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",
    !       put in more-carefully-tuned tolerance, and related changes
    !**************************************************************************
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
    ! Pathname: $Source: /afs/isis/depts/cep/emc/apps/archive/edss_tools/edss_tools/src/lib/ioapi_grd_size.f
    ! Last updated: $Date$
    !
    !***************************************************************************

    IMPLICIT NONE

    !.....  Function arguments
    REAL, INTENT(IN) :: P
    REAL, INTENT(IN) :: Q

    !....  function body  .......

    FLTERR = ((P - Q)**2 .GT. 1.0E-10*( P*P + Q*Q + 1.0E-5 ))

    RETURN

END FUNCTION FLTERR

    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

LOGICAL FUNCTION DBLERR( PD, QD )

    IMPLICIT NONE

    !.....  Function arguments
    REAL(8), INTENT(IN) :: PD
    REAL(8), INTENT(IN) :: QD

    !....  function body  .......

    DBLERR = ((PD - QD)**2 .GT. 1.0D-10*( PD*PD + QD*QD + 1.0D-5 ))

    RETURN

END FUNCTION DBLERR

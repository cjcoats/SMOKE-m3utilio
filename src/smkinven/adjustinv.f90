
SUBROUTINE ADJUSTINV( NRAWBP, UDEV, YDEV, CDEV, LDEV )

    !**************************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine reads the area-to-point file and processes the
    !      sources as needed. It also reads and assigns the non-HAP inclusions/exclusions.

    !  PRECONDITIONS REQUIRED:
    !      CSOURC, POLVAL, and NPCNT arrays allocated and populated in sorted order
    !      NSRC is set.
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 11/02 by C. Seppanen
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
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

    !.......   MODULES for public variables

    IMPLICIT NONE

    !.......   SUBROUTINE ARGUMENTS
    INTEGER , INTENT (IN) :: NRAWBP      ! no. raw records by pollutant
    INTEGER , INTENT (IN) :: UDEV        ! unit no. for non-HAP inclusions/exclusions
    INTEGER , INTENT (IN) :: YDEV        ! unit no. for ar-to-point
    INTEGER , INTENT (IN) :: CDEV        ! SCC descriptions unit no.
    INTEGER , INTENT (IN) :: LDEV        ! log file unit no.

    !***********************************************************************
    !   begin body of subroutine ADJUSTINV

    !......  If area-to-point factors file is present...
    IF( YDEV .GT. 0 ) THEN

        !.......  Read and preprocess area-to-point factors file
        !.......  Result of this call is that the NAR2PT and AR2PTABL
        !         arrays from MODAR2PT and the CHRT09 and ARPT09 arrays
        !         from MODLISTS will be populated.
        CALL RDAR2PT( YDEV, CDEV, LDEV )

        !.......  Assign area-to-point cross-reference entries to sources
        !.......  Result of this call is that the AR2PTTBL, AR2PTIDX, and
        !         AR2PTCNT arrays from MODLISTS will be populated
        CALL ASGNAR2PT

        !.......  Process area-to-point sources
        CALL PROCAR2PT( NRAWBP )

    END IF

    !......  If non-HAP inclusions/exclusions file is present...
    IF( UDEV .GT. 0 ) THEN

        !.......  Read and preprocess NONHAPVOC inclusions/exclusions x-ref
        !.......  Only the CHRT* arrays of the MODXREF will be populated,
        !         because we only need to identify the sources, not assign
        !         anything to them.
        CALL RDXCLUDE( UDEV )

        !.......   Assign array for non-HAP inclusions/exclusions
        CALL ASGNNHAPX

    END IF

    RETURN

END SUBROUTINE ADJUSTINV

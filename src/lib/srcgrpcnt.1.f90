
SUBROUTINE SRCGRPCNT( NSRC, NMAT1, NX, IX )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!     This subroutine matches source groups to grid cells to count the
!     total number of data points. This subroutine only exists because
!     of the way the gridding matrix is stored in memory.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!     Created 7/2013 by C. Seppanen
!
!***************************************************************************
!
! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
!                System
! File: @(#)$Id$
!
! COPYRIGHT (C) 2013, Environmental Modeling for Policy Development
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
!**************************************************************************

!.........  MODULES for public variables
    USE MODMERGE, ONLY: ISRCGRP, NSGOUTPUT, GRPCNT, PFLAG

!.........  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: NGRID

!.........  This module contains arrays for plume-in-grid and major sources
    USE MODELEV, ONLY: LMAJOR, LPING

    IMPLICIT NONE

!...........   SUBROUTINE ARGUMENTS
    INTEGER, INTENT(IN) :: NSRC         ! number of sources
    INTEGER, INTENT(IN) :: NMAT1        ! dim 1 for gridding matrix
    INTEGER, INTENT(IN) :: NX( NGRID )  ! no. of sources per cell
    INTEGER, INTENT(IN) :: IX( NMAT1 )  ! list of sources per cell

!...........   Other local variables
    INTEGER C, J, K, SRC, GIDX, CNT   ! counters and indices

!**********************************************************************
!   begin body of subroutine SRCGRPCNT

    K = 0
    DO C = 1, NGRID
        DO J = 1, NX( C )
            K = K + 1
            SRC = IX( K )

!.................  Skip elevated sources
            IF( PFLAG ) THEN
                IF( LMAJOR( SRC ) .OR. LPING( SRC ) ) CYCLE
            END IF

            GIDX = ISRCGRP( SRC )
            CNT = GRPCNT( C, GIDX )
            IF( CNT == 0 ) THEN
                NSGOUTPUT = NSGOUTPUT + 1
            END IF
            GRPCNT( C, GIDX ) = CNT + 1
        END DO
    END DO

    RETURN

END SUBROUTINE SRCGRPCNT

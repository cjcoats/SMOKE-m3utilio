
SUBROUTINE APPLUMAT( NSRC, NMATX, VAL, NU, IU, CU, VALBYSRC )

    !***********************************************************************
    !  subroutine APPLUMAT body starts at line 72
    !
    !  DESCRIPTION:
    !      Applies the "ungridding" matrix to gridded data to compute a per-source
    !      value of the data.  If the ungridding matrix has no factors for a source,
    !      a missing value is returned for that source (i.e., BADVAL3)
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION HISTORY:
    !       Version ??/???? by ???
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***************************************************************************
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
    !****************************************************************************
    USE M3UTILIO

    !.....  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: NCOLS, YOFF, XDIFF, XOFF

    IMPLICIT NONE

    !.......   SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN) :: NSRC              ! no. sources
    INTEGER     , INTENT (IN) :: NMATX             ! ugridding factors dimension
    REAL        , INTENT (IN) :: VAL( * )          ! gridded data
    INTEGER     , INTENT (IN) :: NU( NSRC  )       ! no. cells per source
    INTEGER     , INTENT (IN) :: IU( NMATX )       ! grid cell IDs
    REAL        , INTENT (IN) :: CU( NMATX )       ! ungridding coefficients
    REAL        , INTENT(OUT) :: VALBYSRC( NSRC )     ! ungridded data

    !.......   Other local variables
    INTEGER     C, J, K, S         ! counters and indices
    INTEGER     COL                ! subgrid column number
    INTEGER     ROW                ! subgrid row number

    REAL        RDUM               ! tmp value for summing over cells

    CHARACTER(16), PARAMETER :: PROGNAME = 'APPLUMAT'     ! program name

    !***********************************************************************
    !   begin body of subroutine APPLUMAT

    !.....  Apply ungridding matrix from a (possible) subgrid to data on base
    !           grid.  If no subgrid, then XOFF and YOFF will be 1 and no problem.
    K = 0
    DO S = 1, NSRC

        IF( NU( S ) .GT. 0 ) THEN
            RDUM = 0.0
        ELSE
            RDUM = BADVAL3
        END IF

        DO J = 1, NU( S )
            K = K + 1

            !.....  Get column and row from subgrid
            C = IU( K )
            ROW = 1 + ( C - 1 ) / NCOLS                      ! note: integer math
            COL = C - ( ROW-1 ) * NCOLS

            !.....  Compute cell number of base grid
            C = ( ROW + YOFF - 1 ) * ( NCOLS + XDIFF ) + COL + XOFF

            RDUM = RDUM + VAL( C ) * CU( K )
        END DO

        VALBYSRC( S ) = RDUM

    END DO

    RETURN

END SUBROUTINE APPLUMAT

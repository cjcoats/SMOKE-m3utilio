
SUBROUTINE GRDFIPS( NSRC, CNTY, VAL, VALBYSRC, TFLAG )

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
    !       Created ??/???? by ???
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

    !.......   Modules for public variables
    !.......   This module contains the gridding surrogates tables
    USE MODSURG, ONLY: NSRGFIPS, SRGFIPS, NCELLS, FIPCELL

    !.....  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: NCOLS

    !.....  This module is the derived meteorology data for emission factors
    USE MODMET, ONLY: MINTSRC, MAXTSRC, MAXTDAY, MINTDAY

    IMPLICIT NONE

    !.......   SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN) :: NSRC                 ! no. sources
    CHARACTER(*), INTENT (IN) :: CNTY( NSRC  )        ! no. counties
    REAL        , INTENT (IN) :: VAL( * )             ! gridded data
    REAL        , INTENT(OUT) :: VALBYSRC( NSRC )     ! ungridded data
    LOGICAL     , INTENT (IN) :: TFLAG                ! processing temp. variable

    !.......   Other local variables
    INTEGER     C, I, J, K, L, N, S      ! counters and indices
    INTEGER     COL                ! subgrid column number
    INTEGER     ROW                ! subgrid row number

    REAL        CNTOT              ! tmp value for summing over cells
    REAL        TEMPVAL            ! temperature value in Farenheight

    !***********************************************************************
    !   begin body of subroutine APPLUMAT

    !.....  Apply ungridding matrix from a (possible) subgrid to data on base
    !           grid.  If no subgrid, then XOFF and YOFF will be 1 and no problem.
    DO S = 1, NSRC

        L = FINDC( CNTY( S ), NSRGFIPS, SRGFIPS )

        IF( L < 1 ) CYCLE

        IF( NCELLS( L ) .GT. 0 ) THEN
            CNTOT = 0.0
        ELSE
            CNTOT = BADVAL3
        END IF

        K = 0
        DO J = 1, NCELLS( L )
            K = K + 1

            !.....  Get column and row from subgrid
            N   = FIPCELL( J,L ) - 1
            ROW = 1 + N / NCOLS                      ! note: integer math
            COL = 1 + MOD ( N, NCOLS )

            ! bbh           C = (ROW-1)*NCOLS + COL               ! org cell equations using row,col

            !.....  Convert K to F
            IF( TFLAG ) THEN
                TEMPVAL = 1.8 * VAL( C ) - 459.67

                MAXTSRC( S ) = MAX( TEMPVAL, MAXTSRC( S ) )
                MINTSRC( S ) = MIN( TEMPVAL, MINTSRC( S ) )

                MAXTDAY( S ) = MAX( TEMPVAL, MAXTDAY( S ) )
                MINTDAY( S ) = MIN( TEMPVAL, MINTDAY( S ) )
            END IF

            CNTOT = CNTOT + VAL( C )

        END DO

        VALBYSRC( S ) = CNTOT / K        ! averaged

    END DO

    RETURN

END SUBROUTINE GRDFIPS

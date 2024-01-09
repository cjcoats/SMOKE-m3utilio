
SUBROUTINE BLDSRCCELL( NSRC, NGRID, NMAT, NX, IX, CX )

    !***********************************************************************
    !  subroutine BLDSRCCELL body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine uses the gridding matrix to build a list of grid
    !      cells and associated fractions for each source.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 4/10 by C. Seppanen
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***********************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !         System
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

    !.......  MODULES for public variables
    !.......  This module contains data structures and flags specific to Movesmrg
    USE MODMVSMRG, ONLY: NSRCCELLS, SRCCELLS, SRCCELLFRACS

    IMPLICIT NONE

    !.......  SUBROUTINE ARGUMENTS
    INTEGER, INTENT(IN) :: NSRC            ! number of sources
    INTEGER, INTENT(IN) :: NGRID           ! number of grid cells
    INTEGER, INTENT(IN) :: NMAT            ! dimension for matrixes
    INTEGER, INTENT(IN) :: NX( NGRID )     ! number of sources per cell
    INTEGER, INTENT(IN) :: IX( NMAT )      ! list of sources per cell
    REAL,    INTENT(IN) :: CX( NMAT )      ! list of source fractions per cell

    !.......  LOCAL VARIABLES and their descriptions:

    !.......  Other local variables
    INTEGER   INDX, NG, NS           ! indexes and counters
    INTEGER   MXNSRCCELLS            ! max. no. cells per source
    INTEGER   CELL                   ! cell number
    INTEGER   SRC                    ! source number
    INTEGER   IOS                    ! error status

    CHARACTER(16), PARAMETER :: PROGNAME = 'BLDSRCCELL'     ! program name

    !***********************************************************************
    !   begin body of subroutine BLDSRCCELL

    !.......  Determine maximum number of cells per source
    ALLOCATE( NSRCCELLS( NSRC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'NSRCCELLS', PROGNAME )
    NSRCCELLS = 0       ! array

    MXNSRCCELLS = 0
    INDX        = 0
    DO NG = 1, NGRID

        !.......  Loop over sources for current grid cell
        DO NS = 1, NX( NG )

            INDX = INDX + 1
            SRC = IX( INDX )
            NSRCCELLS( SRC ) = NSRCCELLS( SRC ) + 1
            IF( NSRCCELLS( SRC ) > MXNSRCCELLS ) THEN
                MXNSRCCELLS = NSRCCELLS( SRC )
            END IF

        END DO

    END DO

    !.......  Store list of cells and fractions for each source
    ALLOCATE(     SRCCELLS( MXNSRCCELLS, NSRC ),    &
              SRCCELLFRACS( MXNSRCCELLS, NSRC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SRCCELLFRACS', PROGNAME )
    NSRCCELLS = 0           ! array
    SRCCELLS = 0            ! array
    SRCCELLFRACS = 0.0      ! array

    INDX = 0
    DO NG = 1, NGRID

        DO NS = 1, NX( NG )

            INDX = INDX + 1
            SRC = IX( INDX )
            NSRCCELLS( SRC ) = NSRCCELLS( SRC ) + 1

            CELL = NSRCCELLS( SRC )
            SRCCELLS    ( CELL, SRC ) = NG
            SRCCELLFRACS( CELL, SRC ) = CX( INDX )

        END DO

    END DO

END SUBROUTINE BLDSRCCELL

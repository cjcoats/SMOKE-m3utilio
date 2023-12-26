
LOGICAL FUNCTION INGRID( XX, YY, NCOLS, NROWS, COL, ROW )

!***************************************************************************
!  function body starts at line
!
!  DESCRIPTION:
!      This function determines if coordinates are inside the grid or not.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created by M. Houyoux 9/2000
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
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
!***************************************************************************

    USE M3UTILIO

    IMPLICIT NONE

!...........   INCLUDES

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   EXTERNAL FUNCTIONS and their descriptions:
    LOGICAL, EXTERNAL :: DSCM3GRD

!...........   SUBROUTINE ARGUMENTS
    REAL   , INTENT (IN) :: XX    ! X-coordinate in projection
    REAL   , INTENT (IN) :: YY    ! Y-ccordinate in projection
    INTEGER, INTENT(OUT) :: NCOLS ! no. of columns
    INTEGER, INTENT(OUT) :: NROWS ! no. of rows
    INTEGER, INTENT(OUT) :: COL   ! column number
    INTEGER, INTENT(OUT) :: ROW   ! row number

!...........   Local variables

    INTEGER          C, I  !  indices and counters.
    INTEGER, SAVE :: NC    !  saved no. columns
    INTEGER, SAVE :: NR    !  saved no. rows

    REAL             XDUM, YDUM ! tmp X and Y coordinates
    REAL, SAVE    :: XX0, YY0   ! X and Y origin in coordinates of grid
    REAL, SAVE    :: XX1, YY1   ! X and Y upper right position of grid
    REAL, SAVE    :: DDX, DDY   ! Inverse cell length in X and Y directions

    LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! true: first time routine called

    CHARACTER(16)   COORD     !  coordinate system name
    CHARACTER(16)   COORUNIT  !  coordinate system projection units
    CHARACTER(16)   GRDNM     !  grid name
    CHARACTER(80)   GDESC     !  grid description
    CHARACTER(300)  MESG      !  message buffer

    CHARACTER(16) :: PROGNAME = 'INGRID' ! program name

!***********************************************************************
!   begin body of function INGRID

!.........  Get grid name from the environment and read grid parameters
    IF( FIRSTIME ) THEN

        IF( .NOT. DSCM3GRD( GRDNM, GDESC, COORD, GDTYP3D, COORUNIT,&
        &                P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,&
        &                XORIG3D, YORIG3D, XCELL3D, YCELL3D,&
        &                NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            MESG = 'Could not get Models-3 grid description'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

!.............  Store grid parameters for later processing
        ELSE

            XX0 = SNGL( XORIG3D )
            YY0 = SNGL( YORIG3D )
            XX1 = XX0 + FLOAT( NCOLS3D ) * SNGL( XCELL3D )
            YY1 = YY0 + FLOAT( NROWS3D ) * SNGL( YCELL3D )
            DDX = 1.0 / SNGL( XCELL3D )
            DDY = 1.0 / SNGL( YCELL3D )
            NR  = NROWS3D
            NC  = NCOLS3D

        END IF

        FIRSTIME = .FALSE.

    END IF

    NROWS = NR
    NCOLS = NC

!.........  Initialize function
    INGRID = .TRUE.

!.........  Check X starting coordinate
    XDUM = DDX * ( XX - XX0 )
    IF ( XDUM .LT. 0.0  )  THEN
        INGRID = .FALSE.
        RETURN
    END IF

!.........  Check X cell
    COL = 1 + INT( XDUM )
    IF ( COL .GT. NCOLS .OR. COL .LE. 0 )  THEN
        INGRID = .FALSE.
        RETURN
    END IF

!.........  Check Y starting coordinate
    YDUM = DDY * ( YY - YY0 )
    IF ( YDUM .LT. 0.0  )  THEN
        INGRID = .FALSE.
        RETURN
    END IF

!.........  Check Y cell
    ROW = 1 + INT( YDUM )
    IF ( ROW .GT. NROWS .OR. ROW .LE. 0 )  THEN
        INGRID = .FALSE.
        RETURN
    END IF

    RETURN

END FUNCTION INGRID


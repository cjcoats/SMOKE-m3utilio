
MODULE MODGRDLIB

    !***********************************************************************
    !
    !   DESCRIPTION:
    !   Generic-F90 routines for coordinate conversion and related tasks
    !
    !       CONVRTLL:  convert from grid-projection coords to Lat-Lon
    !
    !       CONVRTXY:  convert from Lat-Lon to grid-projection coords
    !
    !       LNK2GRD:  Intersect a link (XBEG,YBEG)--(XEND,YEND) with a grid,
    !       and return the number NCEL of cells intersected by the link,
    !       the cell-indices K=col + (row-1)*ncols), and the fractions
    !       AFRAC( I ) of the link in each cell:link intersections.
    !
    !       UNGRIDBV:  compute "ungridding" matrices for I/O API routine BMATVEC(),
    !       based on the grid cell coordinates in a variable grid
    !
    !  REVISION  HISTORY:
    !       Prototype 12/2017-1/2018 by Carlie J. Coats, Jr., UNC IE
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !****************************************************************************/
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

    USE M3UTILIO, M3U_INITSPHERES => INITSPHERES
    USE MODGCTP

    IMPLICIT NONE

    PUBLIC :: CONVRTLL, CONVRTXY, DBLERR, FLTERR, INGRID, LNK2GRD, UNGRIDBV

        !    !......  Generic interfaces:

    INTERFACE CONVRTLL
        MODULE PROCEDURE  CONVRTLLD1, CONVRTLLR1, CONVRTLLD2, CONVRTLLR2
    END INTERFACE CONVRTLL

    INTERFACE CONVRTXY
        MODULE PROCEDURE  CONVRTXYD1, CONVRTXYR1, CONVRTXYD2, CONVRTXYR2
    END INTERFACE CONVRTXY


        !    !  FLTERR(), DBLERR():  "it is not reasonable that the two arguments
        !    !   are intended for the same thing, even though they may have been
        !    !  calculated on different machines, and using different precisions."
    INTERFACE FLTERR
        MODULE PROCEDURE  RRERR, RDERR, DRERR, DDERR
    END INTERFACE FLTERR

    INTERFACE DBLERR
        MODULE PROCEDURE  RRERR, RDERR, DRERR, DDERR
    END INTERFACE DBLERR

    INTERFACE INGRID
        MODULE PROCEDURE  INGRIDD, INGRIDR
    END INTERFACE INGRID

    INTERFACE LNK2GRD
        MODULE PROCEDURE  LNK2GRDD, LNK2GRDR
    END INTERFACE LNK2GRD

    INTERFACE UNGRIDBV
        MODULE PROCEDURE  UNGRIDBVD, UNGRIDBVR
    END INTERFACE UNGRIDBV


    PRIVATE             !    !  everything else

    LOGICAL :: FIRSTIME = .TRUE.


CONTAINS        !    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


    SUBROUTINE CONVRTLLD1( NSRC, CTYPE, GDNAM, P_ALP, P_BET,    &
                           P_GAM, XCENT, YCENT, XVALS, YVALS )

        !.......   INCLUDES

        INCLUDE 'EMCNST3.h90'           !  emissions constant parameters

        !.......  SUBROUTINE ARGUMENTS
        INTEGER,            INTENT (IN) :: NSRC                  !  actual source count
        INTEGER,            INTENT (IN) :: CTYPE                 !  coord sys type
        CHARACTER(NAMLEN3), INTENT (IN) :: GDNAM                 !  grid name
        REAL(8),            INTENT (IN) :: P_ALP                 !  first, second, third map
        REAL(8),            INTENT (IN) :: P_BET                 !  projection descriptive
        REAL(8),            INTENT (IN) :: P_GAM                 !  parameters
        REAL(8),            INTENT (IN) :: XCENT                 !  lon for coord-system X=0
        REAL(8),            INTENT (IN) :: YCENT                 !  lat for coord-system Y=0
        REAL(8),            INTENT(IN OUT) :: XVALS( NSRC )         !  x location (input grid coord)
        REAL(8),            INTENT(IN OUT) :: YVALS( NSRC )         !  y location (input grid coord)

        !.......   Other local variables

        REAL(8)     XLOC( NSRC ), YLOC( NSRC )          ! tmp x and y coordinates

        !.......   Subroutine body ............

        IF ( FIRSTIME ) THEN
            IF ( INITSPHERES() ) THEN
                FIRSTIME = .FALSE.
            ELSE
                CALL M3EXIT( 'MODGRDLIB/CONVRTLL', 0,0, 'Failure in INITSPHERES()', 2 )
            END IF
        END IF

        CALL XY2XY( LATGRD3, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
                    CTYPE, P_ALP, P_BET, P_GAM, XCENT, YCENT,   &
                    NSRC, XVALS, YVALS, XLOC, YLOC )

        XVALS = XLOC
        YVALS = YLOC

        RETURN

    END SUBROUTINE CONVRTLLD1


    !***********************************************************************


    SUBROUTINE CONVRTLLR1( NSRC, CTYPE, GDNAM, P_ALP, P_BET,    &
                           P_GAM, XCENT, YCENT, XVALS, YVALS )

        !.......   INCLUDES

        INCLUDE 'EMCNST3.h90'           !  emissions constant parameters

        !.......  SUBROUTINE ARGUMENTS
        INTEGER,            INTENT (IN) :: NSRC                  !  actual source count
        INTEGER,            INTENT (IN) :: CTYPE                 !  coord sys type
        CHARACTER(NAMLEN3), INTENT (IN) :: GDNAM                 !  grid name
        REAL(8),            INTENT (IN) :: P_ALP                 !  first, second, third map
        REAL(8),            INTENT (IN) :: P_BET                 !  projection descriptive
        REAL(8),            INTENT (IN) :: P_GAM                 !  parameters
        REAL(8),            INTENT (IN) :: XCENT                 !  lon for coord-system X=0
        REAL(8),            INTENT (IN) :: YCENT                 !  lat for coord-system Y=0
        REAL,               INTENT(IN OUT) :: XVALS( NSRC )         !  x location (input grid coord)
        REAL,               INTENT(IN OUT) :: YVALS( NSRC )         !  y location (input grid coord)

        !.......   Other local variables

        REAL(8)     XVAL( NSRC ), YVAL( NSRC )          ! tmp x and y coordinates
        REAL(8)     XLOC( NSRC ), YLOC( NSRC )          ! tmp x and y coordinates

        !.......   Subroutine body ............

        IF ( FIRSTIME ) THEN
            IF ( INITSPHERES() ) THEN
                FIRSTIME = .FALSE.
            ELSE
                CALL M3EXIT( 'MODGRDLIB/CONVRTLL', 0,0, 'Failure in INITSPHERES()', 2 )
            END IF
        END IF

        XVAL = XVALS
        YVAL = YVALS

        CALL XY2XY( LATGRD3, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
                    CTYPE, P_ALP, P_BET, P_GAM, XCENT, YCENT,   &
                    NSRC, XVAL, YVAL, XLOC, YLOC )

        XVALS = XLOC
        YVALS = YLOC

        RETURN

    END SUBROUTINE CONVRTLLR1


    !***********************************************************************


    SUBROUTINE CONVRTLLD2( NCOL, NROW, CTYPE, GDNAM, P_ALP, P_BET,  &
                           P_GAM, XCENT, YCENT, XVALS, YVALS )

        !.......   INCLUDES

        INCLUDE 'EMCNST3.h90'           !  emissions constant parameters

        !.......  SUBROUTINE ARGUMENTS
        INTEGER,            INTENT (IN) :: NCOL, NROW            !  actual source dims
        INTEGER,            INTENT (IN) :: CTYPE                 !  coord sys type
        CHARACTER(NAMLEN3), INTENT (IN) :: GDNAM                 !  grid name
        REAL(8),            INTENT (IN) :: P_ALP                 !  first, second, third map
        REAL(8),            INTENT (IN) :: P_BET                 !  projection descriptive
        REAL(8),            INTENT (IN) :: P_GAM                 !  parameters
        REAL(8),            INTENT (IN) :: XCENT                 !  lon for coord-system X=0
        REAL(8),            INTENT (IN) :: YCENT                 !  lat for coord-system Y=0
        REAL(8),            INTENT(IN OUT) :: XVALS( NCOL,NROW )         !  x location (input grid coord)
        REAL(8),            INTENT(IN OUT) :: YVALS( NCOL,NROW )         !  y location (input grid coord)

        !.......   Other local variables

        REAL(8)     XLOC( NCOL,NROW ), YLOC( NCOL,NROW )          ! tmp x and y coordinates

        !.......   Subroutine body ............

        IF ( FIRSTIME ) THEN
            IF ( INITSPHERES() ) THEN
                FIRSTIME = .FALSE.
            ELSE
                CALL M3EXIT( 'MODGRDLIB/CONVRTLL', 0,0, 'Failure in INITSPHERES()', 2 )
            END IF
        END IF

        CALL XY2XY( LATGRD3, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
                    CTYPE, P_ALP, P_BET, P_GAM, XCENT, YCENT,   &
                    NCOL, NROW, XVALS, YVALS, XLOC, YLOC )

        XVALS = XLOC
        YVALS = YLOC

        RETURN

    END SUBROUTINE CONVRTLLD2


    !***********************************************************************


    SUBROUTINE CONVRTLLR2( NCOL, NROW, CTYPE, GDNAM, P_ALP, P_BET,  &
                           P_GAM, XCENT, YCENT, XVALS, YVALS )

        !.......   INCLUDES

        INCLUDE 'EMCNST3.h90'           !  emissions constant parameters

        !.......  SUBROUTINE ARGUMENTS
        INTEGER,            INTENT (IN) :: NCOL, NROW            !  actual source dims
        INTEGER,            INTENT (IN) :: CTYPE                 !  coord sys type
        CHARACTER(NAMLEN3), INTENT (IN) :: GDNAM                 !  grid name
        REAL(8),            INTENT (IN) :: P_ALP                 !  first, second, third map
        REAL(8),            INTENT (IN) :: P_BET                 !  projection descriptive
        REAL(8),            INTENT (IN) :: P_GAM                 !  parameters
        REAL(8),            INTENT (IN) :: XCENT                 !  lon for coord-system X=0
        REAL(8),            INTENT (IN) :: YCENT                 !  lat for coord-system Y=0
        REAL,               INTENT(IN OUT) :: XVALS( NCOL,NROW )         !  x location (input grid coord)
        REAL,               INTENT(IN OUT) :: YVALS( NCOL,NROW )         !  y location (input grid coord)

        !.......   Other local variables

        REAL(8)     XVAL( NCOL,NROW ), YVAL( NCOL,NROW )          ! tmp x and y coordinates
        REAL(8)     XLOC( NCOL,NROW ), YLOC( NCOL,NROW )          ! tmp x and y coordinates

        !.......   Subroutine body ............

        IF ( FIRSTIME ) THEN
            IF ( INITSPHERES() ) THEN
                FIRSTIME = .FALSE.
            ELSE
                CALL M3EXIT( 'MODGRDLIB/CONVRTLL', 0,0, 'Failure in INITSPHERES()', 2 )
            END IF
        END IF

        XVAL = XVALS
        YVAL = YVALS

        CALL XY2XY( LATGRD3, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
                    CTYPE, P_ALP, P_BET, P_GAM, XCENT, YCENT,   &
                    NCOL, NROW, XVAL, YVAL, XLOC, YLOC )

        XVALS = XLOC
        YVALS = YLOC

        RETURN

    END SUBROUTINE CONVRTLLR2


    !***********************************************************************


    SUBROUTINE CONVRTXYD1( NSRC, CTYPE, GDNAM, P_ALP, P_BET,    &
                           P_GAM, XCENT, YCENT, XVALS, YVALS )

        !.......  SUBROUTINE ARGUMENTS
        INTEGER,      INTENT    (IN) :: NSRC                  !  actual source count
        INTEGER,      INTENT    (IN) :: CTYPE                 !  coord sys type
        CHARACTER(*), INTENT    (IN) :: GDNAM                 !  grid name
        REAL(8),      INTENT    (IN) :: P_ALP                 !  first, second, third map
        REAL(8),      INTENT    (IN) :: P_BET                 !  projection descriptive
        REAL(8),      INTENT    (IN) :: P_GAM                 !  parameters
        REAL(8),      INTENT    (IN) :: XCENT                 !  lon for coord-system X=0
        REAL(8),      INTENT    (IN) :: YCENT                 !  lat for coord-system Y=0
        REAL(8),      INTENT(IN OUT) :: XVALS( NSRC )         !  x location (input lon)
        REAL(8),      INTENT(IN OUT) :: YVALS( NSRC )         !  y location (input lat)

        !.......   Other local variables

        REAL(8)     XLOC( NSRC ), YLOC( NSRC )          ! tmp x and y coordinates

        !.......   Subroutine body ............

        IF ( FIRSTIME ) THEN
            IF ( INITSPHERES() ) THEN
                FIRSTIME = .FALSE.
            ELSE
                CALL M3EXIT( 'MODGRDLIB/CONVRTXY', 0,0, 'Failure in INITSPHERES()', 2 )
            END IF
        END IF

        CALL XY2XY( CTYPE, P_ALP, P_BET, P_GAM, XCENT, YCENT,   &
                    LATGRD3, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
                    NSRC, XVALS, YVALS, XLOC, YLOC )

        XVALS = XLOC
        YVALS = YLOC

        RETURN

    END SUBROUTINE CONVRTXYD1


    !***********************************************************************


    SUBROUTINE CONVRTXYR1( NSRC, CTYPE, GDNAM, P_ALP, P_BET,    &
                           P_GAM, XCENT, YCENT, XVALS, YVALS )

        !.......  SUBROUTINE ARGUMENTS
        INTEGER,      INTENT    (IN) :: NSRC                  !  actual source count
        INTEGER,      INTENT    (IN) :: CTYPE                 !  coord sys type
        CHARACTER(*), INTENT    (IN) :: GDNAM                 !  grid name
        REAL(8),      INTENT    (IN) :: P_ALP                 !  first, second, third map
        REAL(8),      INTENT    (IN) :: P_BET                 !  projection descriptive
        REAL(8),      INTENT    (IN) :: P_GAM                 !  parameters
        REAL(8),      INTENT    (IN) :: XCENT                 !  lon for coord-system X=0
        REAL(8),      INTENT    (IN) :: YCENT                 !  lat for coord-system Y=0
        REAL,         INTENT(IN OUT) :: XVALS( NSRC )         !  x location (input lon)
        REAL,         INTENT(IN OUT) :: YVALS( NSRC )         !  y location (input lat)

        !.......   Other local variables

        REAL(8)     XVAL( NSRC ), YVAL( NSRC )          ! tmp x and y coordinates
        REAL(8)     XLOC( NSRC ), YLOC( NSRC )          ! tmp x and y coordinates

        !.......   Subroutine body ............

        IF ( FIRSTIME ) THEN
            IF ( INITSPHERES() ) THEN
                FIRSTIME = .FALSE.
            ELSE
                CALL M3EXIT( 'MODGRDLIB/CONVRTXY', 0,0, 'Failure in INITSPHERES()', 2 )
            END IF
        END IF

        XVAL = XVALS
        YVAL = YVALS

        CALL XY2XY( CTYPE, P_ALP, P_BET, P_GAM, XCENT, YCENT,   &
                    LATGRD3, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
                    NSRC, XVAL, YVAL, XLOC, YLOC )

        XVALS = XLOC
        YVALS = YLOC

        RETURN

    END SUBROUTINE CONVRTXYR1


    !***********************************************************************


    SUBROUTINE CONVRTXYD2( NCOL, NROW, CTYPE, GDNAM, P_ALP, P_BET,  &
                           P_GAM, XCENT, YCENT, XVALS, YVALS )

        !.......  SUBROUTINE ARGUMENTS
        INTEGER,      INTENT    (IN) :: NCOL, NROW            !  actual source dims
        INTEGER,      INTENT    (IN) :: CTYPE                 !  coord sys type
        CHARACTER(*), INTENT    (IN) :: GDNAM                 !  grid name
        REAL(8),      INTENT    (IN) :: P_ALP                 !  first, second, third map
        REAL(8),      INTENT    (IN) :: P_BET                 !  projection descriptive
        REAL(8),      INTENT    (IN) :: P_GAM                 !  parameters
        REAL(8),      INTENT    (IN) :: XCENT                 !  lon for coord-system X=0
        REAL(8),      INTENT    (IN) :: YCENT                 !  lat for coord-system Y=0
        REAL(8),      INTENT(IN OUT) :: XVALS( NCOL,NROW )         !  x location (input lon)
        REAL(8),      INTENT(IN OUT) :: YVALS( NCOL,NROW )         !  y location (input lat)

        !.......   Other local variables

        REAL(8)     XLOC( NCOL,NROW ), YLOC( NCOL,NROW )          ! tmp x and y coordinates

        !.......   Subroutine body ............

        IF ( FIRSTIME ) THEN
            IF ( INITSPHERES() ) THEN
                FIRSTIME = .FALSE.
            ELSE
                CALL M3EXIT( 'MODGRDLIB/CONVRTXY', 0,0, 'Failure in INITSPHERES()', 2 )
            END IF
        END IF

        CALL XY2XY( CTYPE, P_ALP, P_BET, P_GAM, XCENT, YCENT,   &
                    LATGRD3, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
                    NCOL, NROW, XVALS, YVALS, XLOC, YLOC )

        XVALS = XLOC
        YVALS = YLOC

        RETURN

    END SUBROUTINE CONVRTXYD2


    !***********************************************************************


    SUBROUTINE CONVRTXYR2( NCOL, NROW, CTYPE, GDNAM, P_ALP, P_BET,  &
                           P_GAM, XCENT, YCENT, XVALS, YVALS )

        !.......  SUBROUTINE ARGUMENTS
        INTEGER,      INTENT    (IN) :: NCOL, NROW            !  actual source dims
        INTEGER,      INTENT    (IN) :: CTYPE                 !  coord sys type
        CHARACTER(*), INTENT    (IN) :: GDNAM                 !  grid name
        REAL(8),      INTENT    (IN) :: P_ALP                 !  first, second, third map
        REAL(8),      INTENT    (IN) :: P_BET                 !  projection descriptive
        REAL(8),      INTENT    (IN) :: P_GAM                 !  parameters
        REAL(8),      INTENT    (IN) :: XCENT                 !  lon for coord-system X=0
        REAL(8),      INTENT    (IN) :: YCENT                 !  lat for coord-system Y=0
        REAL,         INTENT(IN OUT) :: XVALS( NCOL,NROW )         !  x location (input lon)
        REAL,         INTENT(IN OUT) :: YVALS( NCOL,NROW )         !  y location (input lat)

        !.......   Other local variables

        REAL(8)     XLOC( NCOL,NROW ), YLOC( NCOL,NROW )          ! tmp x and y coordinates
        REAL(8)     XVAL( NCOL,NROW ), YVAL( NCOL,NROW )          ! tmp x and y coordinates

        !.......   Subroutine body ............

        IF ( FIRSTIME ) THEN
            IF ( INITSPHERES() ) THEN
                FIRSTIME = .FALSE.
            ELSE
                CALL M3EXIT( 'MODGRDLIB/CONVRTXY', 0,0, 'Failure in INITSPHERES()', 2 )
            END IF
        END IF

        XVAL = XVALS
        YVAL = YVALS

        CALL XY2XY( CTYPE, P_ALP, P_BET, P_GAM, XCENT, YCENT,   &
                    LATGRD3, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
                    NCOL, NROW, XVAL, YVAL, XLOC, YLOC )

        XVALS = XLOC
        YVALS = YLOC

        RETURN

    END SUBROUTINE CONVRTXYR2


    !***********************************************************************


    LOGICAL FUNCTION RRERR( P, Q )
        REAL, INTENT(IN) :: P
        REAL, INTENT(IN) :: Q

        RRERR = DDERR( DBLE( P ), DBLE( Q ) )
        RETURN

    END FUNCTION RRERR


    !***********************************************************************


    LOGICAL FUNCTION RDERR( P, Q )
        REAL   , INTENT(IN) :: P
        REAL(8), INTENT(IN) :: Q

        RDERR = DDERR( DBLE( P ), Q )
        RETURN

    END FUNCTION RDERR


    !***********************************************************************


    LOGICAL FUNCTION DRERR( P, Q )
        REAL(8), INTENT(IN) :: P
        REAL   , INTENT(IN) :: Q

        DRERR = DDERR( P, DBLE( Q ) )
        RETURN

    END FUNCTION DRERR


    !***********************************************************************


    LOGICAL FUNCTION DDERR( PD, QD )
        REAL(8), INTENT(IN) :: PD
        REAL(8), INTENT(IN) :: QD

        DDERR = ((PD - QD)**2 .GT. 1.0D-11*( PD*PD + QD*QD + 1.0D-5 ))
        RETURN

    END FUNCTION DDERR


    !***********************************************************************

    LOGICAL FUNCTION INGRIDD( XX, YY, NCOLS, NROWS, COL, ROW )

        !.......   SUBROUTINE ARGUMENTS
        REAL(8), INTENT (IN) :: XX            ! X-coordinate in projection
        REAL(8), INTENT (IN) :: YY            ! Y-ccordinate in projection
        INTEGER, INTENT(OUT) :: NCOLS         ! no. of columns
        INTEGER, INTENT(OUT) :: NROWS         ! no. of rows
        INTEGER, INTENT(OUT) :: COL           ! column number
        INTEGER, INTENT(OUT) :: ROW           ! row number

        !.......   Local declarations

        LOGICAL, EXTERNAL :: DSCM3GRD

        INTEGER          C, I          !  indices and counters.
        INTEGER, SAVE :: NC            !  saved no. columns
        INTEGER, SAVE :: NR            !  saved no. rows

        REAL(8)             XDUM, YDUM         ! tmp X and Y coordinates
        REAL(8), SAVE    :: XX0, YY0           ! X and Y origin in coordinates of grid
        REAL(8), SAVE    :: XX1, YY1           ! X and Y upper right position of grid
        REAL(8), SAVE    :: DDX, DDY           ! Inverse cell length in X and Y directions

        LOGICAL, SAVE :: FIRSTIME = .TRUE.           ! true: first time routine called

        CHARACTER(16)   COORD             !  coordinate system name
        CHARACTER(16)   COORUNIT          !  coordinate system projection units
        CHARACTER(16)   GRDNM             !  grid name
        CHARACTER(80)   GDESC             !  grid description
        CHARACTER(300)  MESG              !  message buffer

        CHARACTER(16), PARAMETER :: PROGNAME = 'INGRID'         ! program name

        !.......  begin body of INGRIDD  .........
        IF( FIRSTIME ) THEN

            IF( .NOT. DSCM3GRD( GRDNM, GDESC, COORD, GDTYP3D, COORUNIT,     &
                            P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,    &
                            XORIG3D, YORIG3D, XCELL3D, YCELL3D,             &
                            NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

                MESG = 'Could not get Models-3 grid description'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            !.......  Store grid parameters for later processing
            ELSE

                XX0 = XORIG3D
                YY0 = YORIG3D
                XX1 = XX0 + DBLE( NCOLS3D ) * XCELL3D
                YY1 = YY0 + DBLE( NROWS3D ) * YCELL3D
                DDX = 1.0 / XCELL3D
                DDY = 1.0 / YCELL3D
                NR  = NROWS3D
                NC  = NCOLS3D

            END IF

            FIRSTIME = .FALSE.

        END IF

        NROWS = NR
        NCOLS = NC

        !.......  Initialize function
        INGRIDD = .TRUE.

        !.......  Check X starting coordinate
        XDUM = DDX * ( XX - XX0 )
        IF ( XDUM .LT. 0.0D0  )  THEN
            INGRIDD = .FALSE.
            RETURN
        END IF

        !.......  Check X cell
        COL = 1 + INT( XDUM )
        IF ( COL .GT. NCOLS .OR. COL .LE. 0 )  THEN
            INGRIDD = .FALSE.
            RETURN
        END IF

        !.......  Check Y starting coordinate
        YDUM = DDY * ( YY - YY0 )
        IF ( YDUM .LT. 0.0D0  )  THEN
            INGRIDD = .FALSE.
            RETURN
        END IF

        !.......  Check Y cell
        ROW = 1 + INT( YDUM )
        IF ( ROW .GT. NROWS .OR. ROW .LE. 0 )  THEN
            INGRIDD = .FALSE.
            RETURN
        END IF

        RETURN

    END FUNCTION INGRIDD


    !***********************************************************************


    LOGICAL FUNCTION INGRIDR( XX, YY, NCOLS, NROWS, COL, ROW )

        !.......   SUBROUTINE ARGUMENTS
        REAL   , INTENT (IN) :: XX        ! X-coordinate in projection
        REAL   , INTENT (IN) :: YY        ! Y-ccordinate in projection
        INTEGER, INTENT(OUT) :: NCOLS     ! no. of columns
        INTEGER, INTENT(OUT) :: NROWS     ! no. of rows
        INTEGER, INTENT(OUT) :: COL       ! column number
        INTEGER, INTENT(OUT) :: ROW       ! row number

        INGRIDR = INGRIDD( DBLE( XX ), DBLE( YY ), NCOLS, NROWS, COL, ROW )

        RETURN

    END FUNCTION INGRIDR


    !***********************************************************************


    SUBROUTINE LNK2GRDD( NDIM, XBEGIN, YBEGIN, XENDIN, YENDIN,  &
                         NCEL, ACEL, AFRAC, ALEN, EFLAG )


        !.......   MODULES for public variables
        !.......  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: XCELL, YCELL, XORIG, YORIG, NCOLS, NROWS

        IMPLICIT NONE

        !.......   ARGUMENTS and their descriptions:

        INTEGER, INTENT (IN) :: NDIM                  !  max number of intersections
        REAL(8), INTENT (IN) :: XBEGIN                !  start points of link
        REAL(8), INTENT (IN) :: YBEGIN                !  start points of link
        REAL(8), INTENT (IN) :: XENDIN                !  end points of link
        REAL(8), INTENT (IN) :: YENDIN                !  end points of link
        INTEGER, INTENT(OUT) :: NCEL                  !  number of intersections
        INTEGER, INTENT(OUT) :: ACEL ( NDIM )         !  cell #:  col + (row-1)*ncols
        REAL   , INTENT(OUT) :: AFRAC( NDIM )         !  length of link in cell
        REAL   , INTENT(OUT) :: ALEN                  !  link length
        LOGICAL, INTENT(OUT) :: EFLAG                 !  error flag

        !.......  Local arrays dimensioned with subroutine arguments
        INTEGER :: XCOL(  NDIM+1 )            !  subscript for this grid intersection
        INTEGER :: YROW(  NDIM+1 )            !  subscript for this grid intersection
        REAL(8) :: XFAC(0:NDIM+1 )            !  frac of link traversed
        REAL(8) :: YFAC(0:NDIM+1 )            !  frac of link traversed

        !.......   LOCAL VARIABLES and their descriptions:

        REAL(8)  DDX, DDY                        !  1/cellsize
        REAL(8)  DXLNK                           !  One over x-dir link length
        REAL(8)  DYLNK                           !  One over y-dir link length
        REAL(8)  ENDS                            !  scratch ending coordinate
        REAL(8)  FF                              !  scratch factor vbles
        REAL(8)  FAC                             !  fraction of link traversed so far
        REAL(8)  RSAV
        REAL(8)  START                           !  start starting coordinate
        REAL(8)  XX, YY
        REAL(8)  XBEG, YBEG                      !  local beginning coords
        REAL(8)  XEND, YEND                      !  local ending coords
        REAL(8)  XLNK                            !  link length in x direction
        REAL(8)  YLNK                            !  link length in y direction

        INTEGER  CCC                             !  ending   x-cell of link
        INTEGER  COL                             !  starting x-cell of link
        INTEGER  IX, IY
        INTEGER  IOS                             !  i/o status
        INTEGER  ISAV
        INTEGER  J, IR, IC, JZ
        INTEGER  NX, NY
        INTEGER  ROW                             !  starting y-cell of link
        INTEGER  RRR                             !  ending   y-cell of link
        INTEGER  XINC, YINC                      !  merge counters increments
        INTEGER  XB, XE                          !  pntr to calc fracs @ link beg  end
        INTEGER  YB, YE                          !  pntr to calc fracs @ link beg  end

        CHARACTER(16), PARAMETER :: PNAME = 'LNK2GRDD'          ! subroutine name

        !***********************************************************************
        !   begin body of subroutine  LNK2GRDD

        !.......   Initializations
        EFLAG = .FALSE.

        XBEG = XBEGIN
        XEND = XENDIN
        YBEG = YBEGIN
        YEND = YENDIN

        XCOL = 0           ! array
        YROW = 0           ! array
        XFAC = 0.          ! array
        YFAC = 0.          ! array

        XLNK = XEND - XBEG
        YLNK = YEND - YBEG

        DDX   = 1.0 / XCELL
        DDY   = 1.0 / YCELL

        COL  = 1 + INT( DDX * ( XBEG - XORIG ) )              !  starting cell
        ROW  = 1 + INT( DDY * ( YBEG - YORIG ) )

        CCC  = 1 + INT( DDX * ( XEND - XORIG ) )              !  ending cell
        RRR  = 1 + INT( DDY * ( YEND - YORIG ) )

        ALEN = SQRT( XLNK**2 + YLNK**2 )                        !  link length

        !.......  Initialize the number of cells to zero
        NCEL = 0

        !.......   Check for 0-D, 1-D cases, using tolerance of
        !.......   1.0E-5 * cellsize (on the order of 0.5-10.0 meters).
        !.......   Within a cell case:
        IF ( ( CCC .EQ. COL  .AND. RRR .EQ. ROW ) .OR.  &
             ( ALEN .LT. 1.0E-5 * XCELL )            ) THEN

            !.......  Make sure 1-cell link is inside domain
            IF ( COL .GE. 1  .AND.  COL .LE. NCOLS .AND.    &
                 ROW .GE. 1  .AND.  ROW .LE. NROWS      ) THEN

                NCEL = 1
                ACEL ( 1 ) = COL + NCOLS * ( ROW - 1 )
                AFRAC( 1 ) = 1.0

            END IF

            RETURN

        END IF                     !  O-D or within-a-cell case

        !.......   Calculate once and only once for remaining cases

        DXLNK = 1.0 / ABS( XLNK )          !  safe "/" -- 0-D case already done.
        DYLNK = 1.0 / ABS( YLNK )          !  safe "/" -- 0-D case already done.

        !.......   1-D problems:

        !.......  Case:  (Near)constant X -- 1-D problem in Y only

        IF ( ( COL .EQ. CCC )  .OR.                 &
             ( ABS( XLNK ) .LT. 1.0E-5 * XCELL ) ) THEN

            !.......  Case:  Link is outside domain
            IF ( COL .LT. 1 .OR. COL .GT. NCOLS ) THEN
                NCEL = 0
                RETURN
            END IF

            !.......  For cases in which YBEG is not located in a grid cell lower
            !.......  than YEND, then flip the coordinates and use the same algorithm.

            IF ( ROW .GT. RRR ) THEN

                ISAV = RRR
                RRR  = ROW
                ROW  = ISAV

                RSAV = YEND
                YEND = YBEG
                YBEG = RSAV

            ENDIF

            !.......  Process intersections moving from bottom to top of domain

            !.......  Case: Link outside of domain
            IF ( RRR .LT. 1  .OR.  ROW .GT. NROWS ) THEN
                NCEL = 0
                RETURN

            !.......  Case: Link crosses domain at some point
            ELSE

                !.......  Initialize for all cells (inside domain) intersecting link
                J  = 0
                FF = DYLNK * YCELL
                DO  11  IR = MAX( 1,ROW ) , MIN( RRR,NROWS )

                    IF( IR .GE. 1 .AND. IR .LE. NROWS ) THEN
                        J  = J + 1
                        ACEL ( J ) = COL + NCOLS * ( IR - 1 )
                        AFRAC( J ) = FF
                    ENDIF

11              CONTINUE

                NCEL = J

                !.......  Reset first fraction if it starts on interior of domain
                IF( YBEG .GT. YORIG ) THEN
                    YY = YORIG + YCELL * FLOAT( ROW )
                    AFRAC( 1 ) = DYLNK * ( YY - YBEG )
                    IF( AFRAC( 1 ) < 0.0 ) AFRAC( 1 ) = 0.0
                ENDIF

                IF( NCEL .GT. NDIM ) CALL REPORT_BAD_CELL( NCEL, PNAME, EFLAG )

                !.......  Reset last  fraction if it ends   on interior of domain
                IF( YEND .LT. YORIG + YCELL * NROWS ) THEN
                    YY = YORIG + YCELL * FLOAT( RRR - 1 )
                    AFRAC( NCEL ) = DYLNK * ( YEND - YY )
                    IF( AFRAC( NCEL ) < 0.0 ) AFRAC( NCEL ) = 0.0
                ENDIF

                RETURN

            END IF              !  if link crosses over domain at all

        END IF                  !  if link crosses cells in Y direction only


        !.......  Case: (Near)-constant Y -- 1-D problem in X only:

        IF ( ( ROW .EQ. RRR ) .OR.              &
             ( ABS( YLNK ) .LT. 1.0E-5 * YCELL ) ) THEN

            !.......  Case:  Link is outside domain
            IF ( ROW .LT. 1  .OR. ROW .GT. NROWS ) THEN
                NCEL = 0
                RETURN
            END IF

            !.......  For cases in which XBEG is not located in a grid cell to the left
            !.......  of XEND, then flip the coordinates and use the same algorithm.

            IF ( COL .GT. CCC ) THEN

                ISAV = CCC
                CCC  = COL
                COL  = ISAV

                RSAV = XEND
                XEND = XBEG
                XBEG = RSAV

            ENDIF

            !.......  Case:  Link is outside domain
            IF ( CCC .LT. 1  .OR.  COL .GT. NCOLS ) THEN
                NCEL = 0
                RETURN

            !.......  Case: Link crosses domain at some point
            ELSE

                !.......  Initialize for all cells (inside domain) intersecting link
                J  = 0
                FF = DXLNK * XCELL
                JZ = NCOLS * ( ROW - 1 )
                DO  22  IC = MAX( 1,COL ) , MIN( CCC,NCOLS )
                    J  = J + 1
                    ACEL ( J ) = IC + JZ
                    AFRAC( J ) = FF
22              CONTINUE
                NCEL = J

                IF( NCEL .GT. NDIM ) CALL REPORT_BAD_CELL( NCEL, PNAME, EFLAG )

                !.......  Reset first fraction if it starts on interior of domain
                IF( XBEG .GT. XORIG ) THEN
                    XX = XORIG + XCELL * FLOAT( COL )
                    AFRAC( 1 ) = DXLNK * ( XX - XBEG )
                    IF( AFRAC( 1 ) < 0.0 ) AFRAC( 1 ) = 0.0
                ENDIF

                !.......  Reset last  fraction if it ends   on interior of domain
                IF( XEND .LT. XORIG + XCELL * NCOLS ) THEN
                    XX = XORIG + XCELL * FLOAT( CCC - 1 )
                    AFRAC( NCEL ) = DXLNK * ( XEND - XX )
                    IF( AFRAC( NCEL ) < 0.0 ) AFRAC( NCEL ) = 0.0
                ENDIF

                RETURN

            END IF              !  if link crosses over domain at all

        END IF                  !  if link crosses cells in X direction only

        !...............
        !.......  2-D CASE:
        !.......  (by construction, we know that | XLNK |, | YLNK | > 1.0E5*CELL)
        !.......  Find the intersections of link with the grid, starting
        !.......  from XBEG, YBEG:
        !...............

        !...............
        !.......  Set fractions of link between vertical grid lines by
        !.......  using the X fractions
        !...............

        !.......  Establish loop bounds and starting and ending coordinates
        !.......  Also, precalculate pointers (XB  XE) for adjusting fractions
        !.......  based on orientation of link.
        IX    = COL
        NX    = CCC
        START = XBEG
        ENDS  = XEND

        IF ( COL .LT. CCC ) THEN
            XINC = 1
            XB   = IX
            XE   = NX - 1
        ELSE
            XINC = -1
            XB   = IX - 1
            XE   = NX
        END IF

        J         = 0                 !  cell counter
        XFAC( J ) = 9.0E36            !  sentinel on front end-of-list

        !.......  Initialize for all cells (inside domain) intersecting link
        FF = DXLNK * XCELL

        DO  IC = IX, NX, XINC

            !.......  Skip cell ID is outside domain
            IF( IC < 1 ) CYCLE

            !.......  First X-cell of link is inside domain
            IF( IC .EQ. IX .AND. START .GE. XORIG .AND.     &
                IC .GE. 1  .AND. IC    .LE. NCOLS       ) THEN
                J         = J + 1
                XX        = XORIG + XCELL * FLOAT( XB )
                XCOL( J ) = IC
                XFAC( J ) = DXLNK * FLOAT( XINC ) * ( XX - START )

            !.......  Initialize first cell when link starts outside domain
            ELSEIF( IC .EQ. 1 ) THEN
                J         = J + 1
                XCOL( J ) = IC
                XFAC( J ) = FF

            !.......  Last X-cell of link is inside domain
            ELSEIF( IC   .EQ. NX .AND.                          &
                    ENDS .LE. XORIG + XCELL * NCOLS .AND.       &
                    IC   .GE. 1  .AND. IC .LE. NCOLS          ) THEN
                J         = J + 1
                XX        = XORIG + XCELL * FLOAT( XE )
                XCOL( J ) = IC
                XFAC( J ) = XFAC( J-1 ) + DXLNK * FLOAT( XINC ) * ( ENDS - XX )

            !.......  Set fractions for interior of domain and non-end of link
            ELSEIF( IC .GT. 1 .AND. IC .LE. NCOLS ) THEN
                J         = J + 1
                XCOL( J ) = IC
                XFAC( J ) = XFAC( J-1 ) + FF

            ENDIF

            IF( XFAC( J ) < 0.0 ) XFAC( J ) = 0.0

        END DO

        NX           = J                 !  total number of columns intersected
        XFAC( NX+1 ) = 9.0E36            !  sentinel on tail end-of-list

        !.......  Case: Link is outside domain
        IF ( NX .EQ. 0 ) THEN
            NCEL = 0
            RETURN
        END IF


        !...............
        !.......  Set fractions of link between horizontal grid lines by
        !.......  using the Y fractions
        !...............

        !.......  Establish loop bounds and starting and ending coordinates
        !.......  Also, precalculate pointers (YB  YE) for adjusting fractions
        !.......  based on orientation of link.  Note that YINC is used for
        !.......  stepping through loop _AND_ for changing sign in the YFAC calc
        IY    = ROW
        NY    = RRR
        START = YBEG
        ENDS  = YEND

        IF ( ROW .LT. RRR ) THEN
            YINC = 1
            YB   = IY
            YE   = NY - 1
        ELSE
            YINC = -1
            YB   = IY - 1
            YE   = NY
        END IF

        J         = 0                 !  cell counter
        YFAC( J ) = 9.0E36            !  sentinel on front end-of-list

        !.......  Calculate fractions for cells (inside domain) intersecting link
        FF = DYLNK * YCELL

        DO  IR = IY, NY, YINC

            !.......  Skip cell ID is outside domain
            IF( IR < 1 ) CYCLE

            !.......  First Y-cell of link is inside domain
            IF( IR .EQ. IY .AND. START .GE. YORIG .AND.     &
                IR .GE. 1  .AND. IR    .LE. NROWS      ) THEN
                J         = J + 1
                YY        = YORIG + YCELL * FLOAT( YB )
                YROW( J ) = IR
                YFAC( J ) = DYLNK * FLOAT( YINC ) * ( YY - START )

            !.......  Initialize first cell when link starts outside domain
            ELSEIF( IR .EQ. 1 ) THEN
                J         = J + 1
                YROW( J ) = IR
                YFAC( J ) = FF

            !.......  Last Y-cell of link is inside domain
            ELSEIF( IR   .EQ. NY .AND.                      &
                    ENDS .LE. YORIG + YCELL * NROWS .AND.   &
                    IR   .GE. 1 .AND. IR .LE. NROWS           ) THEN
                J         = J + 1
                YY        = YORIG + YCELL * FLOAT( YE )
                YROW( J ) = IR
                YFAC( J ) = YFAC( J-1 ) + DYLNK * FLOAT( YINC ) * ( ENDS - YY )

            !.......  Set fractions for interior of domain and non-end of link
            ELSEIF( IR .GT. 1 .AND. IR .LE. NROWS ) THEN
                J         = J + 1
                YROW( J ) = IR
                YFAC( J ) = YFAC( J-1 ) + FF

            ENDIF

            IF( YFAC( J ) < 0.0 ) YFAC( J ) = 0.0

        END DO

        NY           = J                 !  total number of columns intersected
        YFAC( NY+1 ) = 9.0E36            !  sentinel on tail end-of-list

        !.......  Case: Link is outside domain
        IF ( NY .EQ. 0 ) THEN
            NCEL = 0
            RETURN
        END IF

        !...............
        !......   Now merge the two intersection lists:
        !......   This algorithm starts at one end of the link,
        !......   and continues to the other.  The XFAC and YFAC variables
        !......   contain the total fraction of the link up to the point
        !......   where the link intersects column IX (for XFAC) and/or
        !......   row IY (for YFAC).   The next intersection produces the
        !......   next AFRAC for cell (IC,IR).  If the link passes through
        !......   a cell corner, then both the IX counter and IY counter
        !......   are incremented by 1.
        !...............

        !.......  Perform the merge.  Loop terminates when both lists hit
        !......  sentinels on either end of fractions arrays

        J  = 0
        FF = 0
        IX = 1
        IY = 1
        XINC = 1
        YINC = 1
133     CONTINUE

        !.......  Intersect new column and new row at same time
        !.......  or both arrays have reached end of lists
        IF( ABS( XFAC( IX ) - YFAC( IY ) ) .LT. 1.0E-05 ) THEN
            FAC = XFAC( IX ) - FF
            FF  = XFAC( IX )
            IC  = XCOL( IX )
            IR  = YROW( IY )
            IX  = IX + XINC
            IY  = IY + YINC

        !.......  Intersect new column next
        ELSE IF ( XFAC( IX ) .LT. YFAC( IY ) ) THEN
            FAC = XFAC( IX ) - FF
            FF  = XFAC( IX )
            IC  = XCOL( IX )
            IR  = YROW( IY )
            IX  = IX + XINC

        !.......  Intersect new row next
        ELSE
            FAC = YFAC( IY ) - FF
            FF  = YFAC( IY )
            IC  = XCOL( IX )
            IR  = YROW( IY )
            IY  = IY + YINC
        END IF
        IF( IR > 0 .AND. IC > 0 ) THEN
            IF ( FF .LE. 2.0 ) THEN           ! Check for sentinel
                J = J + 1
                ACEL ( J ) = IC +  NCOLS * ( IR - 1 )
                AFRAC( J ) = MIN( 1.0, MAX( 0.0, FAC ) )
                GO TO  133                      !  to head of loop
            END IF
        END IF

        !.......   Merge complete.  Return the number of cells found:

        NCEL = J
        RETURN

    END SUBROUTINE LNK2GRDD


    !***********************************************************************


    SUBROUTINE LNK2GRDR( NDIM, XBEGIN, YBEGIN, XENDIN, YENDIN,  &
                        NCEL, ACEL, AFRAC, ALEN, EFLAG )

        !.......   ARGUMENTS and their descriptions:

        INTEGER, INTENT (IN) :: NDIM                  !  max number of intersections
        REAL   , INTENT (IN) :: XBEGIN                !  start points of link
        REAL   , INTENT (IN) :: YBEGIN                !  start points of link
        REAL   , INTENT (IN) :: XENDIN                !  end points of link
        REAL   , INTENT (IN) :: YENDIN                !  end points of link
        INTEGER, INTENT(OUT) :: NCEL                  !  number of intersections
        INTEGER, INTENT(OUT) :: ACEL ( NDIM )         !  cell #:  col + (row-1)*ncols
        REAL   , INTENT(OUT) :: AFRAC( NDIM )         !  length of link in cell
        REAL   , INTENT(OUT) :: ALEN                  !  link length
        LOGICAL, INTENT(OUT) :: EFLAG                 !  error flag

        !.......   LOCAL VARIABLES and their descriptions:

        REAL(8)     X1, X2, Y1, Y2

        !.......   body ..............

        X1 = DBLE( XBEGIN )
        Y1 = DBLE( YBEGIN )
        X2 = DBLE( XENDIN )
        Y2 = DBLE( YENDIN )

        CALL LNK2GRDD( NDIM, X1, Y1, X2, Y2,&
                       NCEL, ACEL, AFRAC, ALEN, EFLAG )

        RETURN

    END SUBROUTINE LNK2GRDR


    !***********************************************************************


    SUBROUTINE REPORT_BAD_CELL( NCEL, PNAME, EFLAG )            !    !  for LNK2GRD*()

        INTEGER,      INTENT(IN   ) :: NCEL
        CHARACTER(*), INTENT(IN   ) :: PNAME
        LOGICAL,      INTENT(INOUT) :: EFLAG

        !.......  Local variables
        CHARACTER(300) MESG

        !..............

        WRITE( MESG,94010 ) 'INTERNAL ERROR: Bad number of cells ', NCEL, 'computed in ' // PNAME
        CALL M3MSG2( MESG )
        EFLAG = .TRUE.

        !.......   Internal buffering formats...... 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

    END SUBROUTINE REPORT_BAD_CELL



    !***********************************************************************


    SUBROUTINE UNGRIDBVD( NC, NR, XREFS, YREFS, NPTS, XLOC, YLOC, NU, CU )

        IMPLICIT NONE

        !.......  SUBROUTINE ARGUMENTS
        INTEGER,      INTENT (IN) :: NC                   ! number of columns in grid
        INTEGER,      INTENT (IN) :: NR                   ! number of rows in grid
        REAL,         INTENT (IN) :: XREFS( NC,NR )         ! grid cell center x coordinates
        REAL,         INTENT (IN) :: YREFS( NC,NR )         ! grid cell center y coordinates
        INTEGER,      INTENT (IN) :: NPTS                 ! number of point-source locations
        REAL(8),      INTENT (IN) :: XLOC( NPTS )         ! X point coordinates
        REAL(8),      INTENT (IN) :: YLOC( NPTS )         ! Y point coordinates
        INTEGER,      INTENT(OUT) :: NU( 4,NPTS )         ! surrounding grid cells for each source
        REAL,         INTENT(OUT) :: CU( 4,NPTS )         ! fraction of each grid cell

        !.......  Local variables
        INTEGER          I, K, S             ! counters and indices
        INTEGER          COL                 ! column for current point
        INTEGER          ROW                 ! row for current point

        REAL(8)          LEFT                ! left edge of matching grid cell
        REAL(8)          RIGHT               ! right edge of matching grid cell
        REAL(8)          BOTTOM              ! bottom edge of matching grid cell
        REAL(8)          TOP                 ! top edge of matching grid cell
        REAL(8)          X, Y                ! normalized difference between point location
                !   and grid cell center
        REAL(8)          P, Q                ! fractions used for calculating coefficients

        CHARACTER(16) :: PNAME = 'UNGRIDBVD'         ! program name

        !***********************************************************************
        !   begin body of subroutine UNGRIDBVD

        !.......  Loop through point-source locations
        DO S = 1, NPTS

            !.......  Find column

            !.......  Check if point is to the left of the grid
            IF( XLOC( S ) < XREFS( 1,1 ) ) THEN
                LEFT = 1
                RIGHT = 1

                X = 0.

            !.......  Check if point is to the right of the grid
            ELSE IF( XLOC( S ) > XREFS( NC, 1 ) ) THEN
                LEFT = NC
                RIGHT = NC

                X = 0.

            !.......  Loop through columns to find the correct one
            ELSE
                DO I = 2, NC
                    IF( XLOC( S ) < XREFS( I,1 ) ) THEN
                        LEFT = I - 1
                        RIGHT = I

                        X = ( XLOC( S )    - XREFS( I-1,1 ) ) / &
                            ( XREFS( I,1 ) - XREFS( I-1,1 ) )
                        EXIT
                    END IF
                END DO
            END IF

            !.......  Find row

            !.......  Check if point is below the grid
            IF( YLOC( S ) < YREFS( 1,1 ) ) THEN
                BOTTOM = 1
                TOP = 1

                Y = 0.

            !.......  Check if point is above the grid
            ELSE IF( YLOC( S ) > YREFS( 1,NR ) ) THEN
                BOTTOM = NR
                TOP = NR

                Y = 0.

            !.......  Loop through rows to find the correct one
            ELSE
                DO I = 2, NR
                    IF( YLOC( S ) < YREFS( 1,I ) ) THEN
                        BOTTOM = I - 1
                        TOP = I

                        Y = ( YLOC( S )    - YREFS( 1,I-1 ) ) /    &
                            ( YREFS( 1,I ) - YREFS( 1,I-1 ) )
                        EXIT
                    END IF
                END DO
            END IF

            !.......  Set grid cells surrounding point location
            NU( 1,S ) = ( BOTTOM - 1 ) * NC + LEFT
            NU( 2,S ) = ( BOTTOM - 1 ) * NC + RIGHT
            NU( 3,S ) = ( TOP - 1 ) * NC + LEFT
            NU( 4,S ) = ( TOP - 1 ) * NC + RIGHT

            !.......  Calculate fractions for each surrounding grid cell
            P = 1. - X
            Q = 1. - Y
            CU( 1,S ) = P * Q
            CU( 2,S ) = X * Q
            CU( 3,S ) = P * Y
            CU( 4,S ) = X * Y

        END DO

        RETURN

    END SUBROUTINE UNGRIDBVD



    !***********************************************************************


    SUBROUTINE UNGRIDBVR( NC, NR, XREFS, YREFS, NPTS, XLOC, YLOC, NU, CU )

        IMPLICIT NONE

        !.......  SUBROUTINE ARGUMENTS
        INTEGER,      INTENT (IN) :: NC                   ! number of columns in grid
        INTEGER,      INTENT (IN) :: NR                   ! number of rows in grid
        REAL,         INTENT (IN) :: XREFS( NC,NR )     ! grid cell center x coordinates
        REAL,         INTENT (IN) :: YREFS( NC,NR )     ! grid cell center y coordinates
        INTEGER,      INTENT (IN) :: NPTS             ! number of point-source locations
        REAL,         INTENT (IN) :: XLOC( NPTS )     ! X point coordinates
        REAL,         INTENT (IN) :: YLOC( NPTS )     ! Y point coordinates
        INTEGER,      INTENT(OUT) :: NU( 4,NPTS )     ! surrounding grid cells for each source
        REAL,         INTENT(OUT) :: CU( 4,NPTS )     ! fraction of each grid cell

        CALL UNGRIDBVD( NC, NR, XREFS, YREFS,   &
                        NPTS, DBLE( XLOC ), DBLE( YLOC ), NU, CU )

        RETURN

    END SUBROUTINE UNGRIDBVR


END MODULE MODGRDLIB

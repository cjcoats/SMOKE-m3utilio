
SUBROUTINE GENPTVCEL( NRECS, NGRID, XLOCA, YLOCA, NEXCLD, NX, INDX, GN, SN )

    !***************************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine creates the scratch gridding matrix, which contains
    !      the cell IDs for each source using a variable grid.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created by C. Seppanen 7/04 based on genptcel.f
    !
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

    !.....  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: GRDNM, GDTYP, NCOLS, NROWS,          &
                       P_ALP, P_BET, P_GAM, XCENT, YCENT,   &
                       XORIG, YORIG, XCELL, YCELL

    USE MODGRDLIB

    IMPLICIT NONE

    !.......   SUBROUTINE ARGUMENTS
    INTEGER, INTENT (IN) :: NRECS              ! no. records w/ coordinates
    INTEGER, INTENT (IN) :: NGRID              ! no. grid cells
    REAL(8), INTENT (IN) :: XLOCA( NRECS )     ! X-coordinate in projection
    REAL(8), INTENT (IN) :: YLOCA( NRECS )     ! Y-ccordinate in projection
    INTEGER, INTENT(OUT) :: NEXCLD             ! no. sources excluded
    INTEGER, INTENT(OUT) :: NX   ( NGRID )     ! no. of sources per cell
    INTEGER, INTENT(OUT) :: INDX ( NRECS )     ! sorting index
    INTEGER, INTENT(OUT) :: GN   ( NRECS )     ! cell number
    INTEGER, INTENT(OUT) :: SN   ( NRECS )     ! record number

    !.......   Local variables
    REAL*8          XVALS( NCOLS+1, NROWS+1 )     ! x values for grid cell boundaries
    REAL*8          YVALS( NCOLS+1, NROWS+1 )     ! y values for grid cell boundaries

    INTEGER         C, I, J      !  indices and counters

    INTEGER         COL       ! tmp column number
    INTEGER         ROW       ! tmp row number
    INTEGER         CSAV      ! saved value of C
    INTEGER         IOS       ! i/o status
    INTEGER         NCDOT     ! number of columns in dot file
    INTEGER         NRDOT     ! number of rows in dot file

    REAL            STEP          ! distance to edge of column or row
    REAL*8          XX, YY        ! tmp X and Y coordinates
    REAL            XMIN          ! minimum location for column
    REAL            XMAX          ! maximum location for column
    REAL            YMIN          ! minimum location for row
    REAL            YMAX          ! maximum location for row

    LOGICAL ::      EFLAG = .FALSE.     ! true: an error has occured

    CHARACTER(16)   GNAME         !  logical file name for GRID_DOT_2D file
    CHARACTER(300)  MESG          !  message buffer

    CHARACTER(16) :: PROGNAME = 'GENPTVCEL'     ! program name

    !***********************************************************************
    !   begin body of subroutine GENPTVCEL

    !.....  Open 2-D grid parameters file to get cell coordinates
    !           Use the GRID_DOT_2D file to get cell boundaries
    GNAME = PROMPTMFILE( 'Enter name for DOT-POINT SURFACE GRID file',  &
                         FSREAD3, 'GRID_DOT_2D', PROGNAME )

    IF( .NOT. DESC3( GNAME ) ) THEN
        MESG = 'Could not get description of file "' // TRIM( GNAME ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.....  Check grid against previously set grid description
    CALL CHKGRID( GNAME, 'DOT', 0, EFLAG )

    IF( EFLAG ) THEN
        MESG = TRIM( GNAME ) // ' does not match previously set ' //&
               'grid description'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    NCDOT = NCOLS + 1
    NRDOT = NROWS + 1

    !.....  Read grid cell coordinates
    IF( .NOT. READ3( GNAME, 'LON', 1, 0, 0, XVALS ) ) THEN
        MESG = 'Could not read LON from file "' // TRIM( GNAME ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF( .NOT. READ3( GNAME, 'LAT', 1, 0, 0, YVALS ) ) THEN
        MESG = 'Could not read LAT from file "' // TRIM( GNAME ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.....  Convert coordinates to map projection units
    CALL CONVRTXY( NCDOT, NRDOT, GDTYP, GRDNM, P_ALP, P_BET, P_GAM,     &
                   XCENT, YCENT, XVALS, YVALS )

    !.....  Initialize number of sources per cell
    NX = 0       ! array

    !.....  Initialize scratch gridding matrix - before sparse storage
    NEXCLD = 0
    CSAV = 0
    DO I = 1, NRECS
        GN  ( I ) = 0
        SN  ( I ) = 0
        INDX( I ) = I

        XX = XLOCA( I )
        YY = YLOCA( I )

        IF( .NOT. INGRID( XLOCA(I), YLOCA(I), NCOLS, NROWS, COL, ROW ) ) THEN
            NEXCLD = NEXCLD + 1
            CYCLE          ! To end of loop
        END IF

        !.....  Find correct column for point
        DO J = 1, NCOLS

            IF( XX >= XVALS( J,1 ) .AND. XX <= XVALS( J+1,1 ) ) EXIT

        END DO

        IF( J <= NCOLS ) THEN
            COL = J
        END IF

        !.....  Find correct row for point
        DO J = 1, NROWS

            IF( YY >= YVALS( 1,J ) .AND. YY <= YVALS( 1,J+1 ) ) EXIT

        END DO

        IF( J <= NROWS ) THEN
            ROW = J
        END IF

        !.....  Compute grid cell number based on column and row
        C = COL + NCOLS * ( ROW - 1 )

        IF( C .LE. NGRID ) THEN
            NX( C ) = NX( C ) + 1

        ELSE IF( C .GT. CSAV ) THEN
            CSAV = C

        END IF

        GN( I ) = C
        SN( I ) = I

    END DO            !  end loop on sources I, computing gridding matrix.

    IF( NCOLS * NROWS .NE. NGRID ) THEN

        MESG = 'INTERNAL ERROR: Number of cells in "' // TRIM( PROGNAME ) //    &
                   '" inconsistent  with calling program'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF ( CSAV .GT. NGRID ) THEN
        WRITE( MESG,94010 )     &
               'INTERNAL ERROR: Number of grid cells C=', CSAV, 'exceeds NGRID=', NGRID
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats.... 94xxx

94010 FORMAT( 10( A, :, I9, :, 1X ) )

END SUBROUTINE GENPTVCEL


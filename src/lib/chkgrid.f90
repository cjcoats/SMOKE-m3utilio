
SUBROUTINE CHKGRID( DATDESC, FTYPE, CHKLEVEL, EFLAG )

    !***********************************************************************
    !  subroutine body starts at line 91
    !
    !  DESCRIPTION:
    !      This subroutine updates the grid information and compares to
    !      the existing information, if it has been previously set.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version  ??/????  by  ???
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  
    !       standard I/O API routines FILCHK3, GRDCHK3, 
    !       Bug-fix:  missing length; use TRIM( FILEDESC ), and related changes
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
    ! Pathname: $Source$
    ! Last updated: $Date$
    !
    !***************************************************************************

    USE M3UTILIO

    !.....  MODULES for public variables
    !.....  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: NCOLS, NROWS, XORIG, YORIG, XOFF, YOFF,  &
                       GDTYP, XCELL, YCELL, XCENT, YCENT,       &
                       P_ALP, P_BET, P_GAM, OFFLAG, GRDNM,      &
                       XOFF_A, YOFF_A, XDIFF, YDIFF, NGRID

    IMPLICIT NONE

    !.....  EXTERNAL FUNCTIONS and their descriptions:
    INTEGER     , EXTERNAL :: GETIFDSC

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT(IN   ) :: DATDESC      ! data descriptions
    CHARACTER(*), INTENT(IN   ) :: FTYPE        ! GMAT|GRID file type of interest
    INTEGER     , INTENT(IN   ) :: CHKLEVEL     ! strigency of check
    LOGICAL     , INTENT(INOUT) :: EFLAG        ! true: comparison failed

    !.....  INCLUDES:
    INCLUDE 'IOCNST3.h90'       !  emissions constant parameters

    !.......   Local parameters
    INTEGER, PARAMETER :: CHK_ALL = 0         ! check all of the grid settings
    INTEGER, PARAMETER :: CHK_SUBGRID = 1     ! check all but allow subgrids
    INTEGER, PARAMETER :: CHK_TMPSUBG = 2     ! check all, allow temporary subgrids

    CHARACTER(16), PARAMETER :: PROGNAME = 'CHKGRID'     ! program name

    !.......   Local variables
    INTEGER       NC          ! tmp number of columns
    INTEGER       NR          ! tmp number of rows
    INTEGER       XO          ! tmp x-offset
    INTEGER       YO          ! tmp y-offset

    REAL(8)       CHK_X       ! tmp val for checking subgrid even with grid
    REAL(8)       CHK_Y       ! tmp val for checking subgrid even with grid
    REAL(8)       X0          ! tmp x origin
    REAL(8)       Y0          ! tmp y origin

    LOGICAL, SAVE :: GFLAG  = .FALSE.     ! true: grid settings have been init
    LOGICAL       :: SFLAG  = .FALSE.     ! true: local error

    CHARACTER(25)   FILDESC      ! description of input file
    CHARACTER(300)  MESG         ! message buffer

    !***********************************************************************
    !   begin body of function CHKGRID

    !.....  Initialize local error flag
    EFLAG = .FALSE.

    !.....  Set tmp rows, columns, and total cells depending on file type
    X0 = XORIG3D
    Y0 = YORIG3D

    IF( FTYPE .EQ. 'GMAT' ) THEN
        NCOLS3D = GETIFDSC( FDESC3D, '/NCOLS3D/', .TRUE. )
        NROWS3D = GETIFDSC( FDESC3D, '/NROWS3D/', .TRUE. )
        FILDESC = 'gridding matrix'
    ELSEIF( FTYPE .EQ. 'GROUPS' ) THEN
        NCOLS3D = GETIFDSC( FDESC3D, '/NCOLS3D/', .TRUE. )
        NROWS3D = GETIFDSC( FDESC3D, '/NROWS3D/', .TRUE. )
        FILDESC = 'stack groups file'
    ELSEIF( FTYPE .EQ. 'GRID' ) THEN
        FILDESC = 'gridded file'
    ELSEIF( FTYPE .EQ. 'GRIDDESC' ) THEN
        FILDESC = 'grid description file'
    ELSEIF( FTYPE .EQ. 'SURROGATES' ) THEN
        FILDESC = 'surrogates file'
    ELSEIF( FTYPE .EQ. 'LANDUSE' ) THEN
        FILDESC = 'landuse file'
    ELSEIF( FTYPE .EQ. 'DOT' ) THEN
        NCOLS3D = NCOLS3D - 1
        NROWS3D = NROWS3D - 1
        XORIG3D = XORIG3D + 0.5 * XCELL3D
        YORIG3D = YORIG3D + 0.5 * YCELL3D
        FILDESC = 'dot-gridded file'
    ELSE
        MESG = 'INTERNAL ERROR: File type "' // FTYPE //    &
                '" not known in call to ' // PROGNAME
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ENDIF

    !.....  If grid information has already been initialized, then compare
    !               existing to this file.

    IF( GFLAG ) THEN

        IF( OFFLAG ) THEN
            X0      = XORIG3D
            Y0      = YORIG3D
            XORIG3D = XORIG3D + DBLE( XOFF ) * XCELL3D
            YORIG3D = YORIG3D + DBLE( YOFF ) * YCELL3D
        END IF

        IF ( .NOT.FILCHK3( FILDESC, FTYPE3D,                    &
                           NCOLS, NROWS, NLAYS3D, NTHIK3D ) ) THEN

            EFLAG = .TRUE.
            CALL M3MSG2( 'Inconsistent dimensions for '//FILDESC )

        ELSE IF ( .NOT.GRDCHK3( FILDESC,                            &
                                P_ALP, P_BET, P_GAM, XCENT, YCENT,  &
                                XORIG, YORIG, XCELL, YCELL,         &
                                NLAYS3D, VGTYP3D, VGTOP3D, VGLVS3D ) ) THEN

            EFLAG = .TRUE.
            CALL M3MSG2( 'Inconsistent dimensions for '//FILDESC )

        ELSE

            EFLAG = .FALSE.

        END IF

    ELSE        !    ! Initialize grid information

        CHK_X  = ( XORIG3D  - XORIG ) / XCELL
        CHK_Y  = ( YORIG3D  - YORIG ) / YCELL

        IF ( DBLERR( CHK_X, DBLE( NINT( CHK_X ) ) ) .OR.    &
             DBLERR( CHK_Y, DBLE( NINT( CHK_Y ) ) ) ) THEN

            EFLAG = .TRUE.
            CALL M3MSG2( 'Inconsistent alignment for subgrid in '//FILDESC )
            RETURN

        ELSE IF ( DBLERR( CHK_X, 0.0D0 ) .OR.   &
                  DBLERR( CHK_Y, 0.0D0 ) ) THEN

            OFFLAG = .TRUE.
            XOFF = NINT( CHK_X  )
            YOFF = NINT( CHK_Y )
            XORIG3D = XORIG3D + DBLE( XOFF ) * XCELL3D
            YORIG3D = YORIG3D + DBLE( YOFF ) * YCELL3D

        END IF

        GFLAG = .TRUE.

        GRDNM = GDNAM3D
        GDTYP = GDTYP3D
        P_ALP = P_ALP3D
        P_BET = P_BET3D
        P_GAM = P_GAM3D
        XCENT = XCENT3D
        YCENT = YCENT3D
        XORIG = XORIG3D
        YORIG = YORIG3D
        XCELL = XCELL3D
        YCELL = YCELL3D
        NCOLS = NCOLS3D
        NROWS = NROWS3D
        NGRID = NCOLS * NROWS

        MESG = 'NOTE: Grid settings initialized using ' //  &
               DATDESC // ' in ' // CRLF()// BLANK10 // TRIM( FILDESC ) // '.'

        CALL M3MSG2( MESG )

        EFLAG = .FALSE.

    ENDIF

    RETURN

CONTAINS        !    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

     !    !  Double-precision "definitely unequal"


    LOGICAL FUNCTION DBLERR( P, Q )
        REAL*8, INTENT( IN ) :: P, Q
        DBLERR = ( (P - Q)**2  .GT.  1.0D-10*( P*P + Q*Q + 1.0D-5 ) )
    END FUNCTION DBLERR

END  SUBROUTINE CHKGRID


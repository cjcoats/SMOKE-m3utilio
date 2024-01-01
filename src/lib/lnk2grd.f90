
SUBROUTINE LNK2GRD( NDIM, XBEGIN, YBEGIN, XENDIN, YENDIN,   &
                    NCEL, ACEL, AFRAC, ALEN, EFLAG )

    !***********************************************************************
    !  subroutine body starts at line 111
    !
    !  FUNCTION:  Given a link with end points XBEG, YBEG, XEND, YEND:
    !       Compute the number NCEL of cells intersected by the link,
    !       the cell-numbers (in terms of storage order col + (row-1)*ncols),
    !       the total length ALEN of the link, and the fraction
    !       AFRAC( I ) of each cell:link intersections.
    !
    !  PRECONDITIONS REQUIRED:
    !       Grid description stored in MODGRID prior to call
    !       right handed coord system (XCELL and YCELL positive)
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Copied from lnk2grd.F 3.2 by mhouyoux
    !       Updated for new SMOKE in 5/99
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***********************************************************************
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

    !.......   MODULES for public variables
    !.....  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: XCELL, YCELL, XORIG, YORIG, NCOLS, NROWS

    IMPLICIT NONE

    !.......   ARGUMENTS and their descriptions:

    INTEGER, INTENT (IN) :: NDIM              !  max number of intersections
    REAL   , INTENT (IN) :: XBEGIN            !  end points of link
    REAL   , INTENT (IN) :: YBEGIN            !  end points of link
    REAL   , INTENT (IN) :: XENDIN            !  end points of link
    REAL   , INTENT (IN) :: YENDIN            !  end points of link
    INTEGER, INTENT(OUT) :: NCEL              !  number of intersections
    INTEGER, INTENT(OUT) :: ACEL ( NDIM )     !  cell #:  col + (row-1)*ncols
    REAL   , INTENT(OUT) :: AFRAC( NDIM )     !  length of link in cell
    REAL   , INTENT(OUT) :: ALEN              !  link length
    LOGICAL, INTENT(OUT) :: EFLAG             !  error flag

    !.....  Local arrays dimensioned with subroutine arguments
    INTEGER, ALLOCATABLE, SAVE :: XCOL( : )        !  subscript for this grid intersection
    INTEGER, ALLOCATABLE, SAVE :: YROW( : )        !  subscript for this grid intersection
    REAL   , ALLOCATABLE, SAVE :: XFAC( : )        !  frac of link traversed
    REAL   , ALLOCATABLE, SAVE :: YFAC( : )        !  frac of link traversed

    !.......   LOCAL VARIABLES and their descriptions:

    REAL     DDX, DDY                    !  1/cellsize
    REAL     DXLNK                       !  One over x-dir link length
    REAL     DYLNK                       !  One over y-dir link length
    REAL     ENDS                        !  scratch ending coordinate
    REAL     FF                          !  scratch factor vbles
    REAL     FAC                         !  fraction of link traversed so far
    REAL     RSAV
    REAL     START                       !  start starting coordinate
    REAL     XX, YY
    REAL     XBEG, YBEG                  !  local beginning coords
    REAL     XEND, YEND                  !  local ending coords
    REAL     XLNK                        !  link length in x direction
    REAL     YLNK                        !  link length in y direction

    INTEGER  CCC                         !  ending   x-cell of link
    INTEGER  COL                         !  starting x-cell of link
    INTEGER  IX, IY
    INTEGER  IOS                         !  i/o status
    INTEGER  ISAV
    INTEGER  J, IR, IC, JZ
    INTEGER  NX, NY
    INTEGER  ROW                         !  starting y-cell of link
    INTEGER  RRR                         !  ending   y-cell of link
    INTEGER  XINC, YINC                  !  merge counters increments
    INTEGER  XB, XE                      !  pntr to calc fracs @ link beg  end
    INTEGER  YB, YE                      !  pntr to calc fracs @ link beg  end

    LOGICAL :: FIRSTIME = .TRUE.         !  true: first time routine called

    CHARACTER(16), PARAMETER :: PROGNAME = 'LNK2GRD'      ! progname name

    !***********************************************************************
    !   begin body of subroutine  LNK2GRD

    !.......   Allocate memory first time routine is called
    IF( FIRSTIME ) THEN
        ALLOCATE( XCOL(   NDIM+1 ),     &
                  YROW(   NDIM+1 ),     &
                  XFAC( 0:NDIM+1 ),     &
                  YFAC( 0:NDIM+1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'YFAC', PROGNAME )
        FIRSTIME = .FALSE.
    END IF

    !.......   Initializations
    XBEG = XBEGIN
    XEND = XENDIN
    YBEG = YBEGIN
    YEND = YENDIN

    XCOL = 0       ! array
    YROW = 0       ! array
    XFAC = 0.      ! array
    YFAC = 0.      ! array

    XLNK = XEND - XBEG
    YLNK = YEND - YBEG

    DDX   = 1.0 / XCELL
    DDY   = 1.0 / YCELL

    COL  = 1 + INT( DDX * ( XBEG - XORIG ) )          !  starting cell
    ROW  = 1 + INT( DDY * ( YBEG - YORIG ) )

    CCC  = 1 + INT( DDX * ( XEND - XORIG ) )          !  ending cell
    RRR  = 1 + INT( DDY * ( YEND - YORIG ) )

    ALEN = SQRT( XLNK**2 + YLNK**2 )                    !  link length

    !.......  Initialize the number of cells to zero
    NCEL = 0

    !.......   Check for 0-D, 1-D cases, using tolerance of
    !.......   1.0E-5 * cellsize (on the order of 0.5-10.0 meters).
    !.......   Within a cell case:
    IF ( ( CCC .EQ. COL  .AND. RRR .EQ. ROW ) .OR. ( ALEN .LT. 1.0E-5 * XCELL ) ) THEN

        !.....  Make sure 1-cell link is inside domain
        IF ( COL .GE. 1  .AND.  COL .LE. NCOLS .AND.        &
             ROW .GE. 1  .AND.  ROW .LE. NROWS      ) THEN

            NCEL = 1
            ACEL ( 1 ) = COL + NCOLS * ( ROW - 1 )
            AFRAC( 1 ) = 1.0

        END IF

        RETURN

    END IF                 !  O-D or within-a-cell case

    !.......   Calculate once and only once for remaining cases
    IF( XLNK == 0.0 ) XLNK = 1.0
    IF( YLNK == 0.0 ) YLNK = 1.0
    DXLNK = 1.0 / ABS( XLNK )      !  safe "/" -- 0-D case already done.
    DYLNK = 1.0 / ABS( YLNK )      !  safe "/" -- 0-D case already done.

    !.......   1-D problems:

    !.....  Case:  (Near)constant X -- 1-D problem in Y only

    IF ( ( COL .EQ. CCC )  .OR. ( ABS( XLNK ) .LT. 1.0E-5 * XCELL ) ) THEN

        !.....  Case:  Link is outside domain
        IF ( COL .LT. 1 .OR. COL .GT. NCOLS ) THEN
            NCEL = 0
            RETURN
        END IF

        !.....  For cases in which YBEG is not located in a grid cell lower
        !.....  than YEND, then flip the coordinates and use the same algorithm.

        IF ( ROW .GT. RRR ) THEN

            ISAV = RRR
            RRR  = ROW
            ROW  = ISAV

            RSAV = YEND
            YEND = YBEG
            YBEG = RSAV

        ENDIF

        !.....  Process intersections moving from bottom to top of domain

        !.....  Case: Link outside of domain
        IF ( RRR .LT. 1  .OR.  ROW .GT. NROWS ) THEN
            NCEL = 0
            RETURN

            !.....  Case: Link crosses domain at some point
        ELSE

            !.....  Initialize for all cells (inside domain) intersecting link
            J  = 0
            FF = DYLNK * YCELL
            DO  11  IR = MAX( 1,ROW ) , MIN( RRR,NROWS )

                IF( IR .GE. 1 .AND. IR. LE. NROWS ) THEN
                    J  = J + 1
                    ACEL ( J ) = COL + NCOLS * ( IR - 1 )
                    AFRAC( J ) = FF
                ENDIF

11          CONTINUE

            NCEL = J

            !.....  Reset first fraction if it starts on interior of domain
            IF( YBEG .GT. YORIG ) THEN
                YY = YORIG + YCELL * FLOAT( ROW )
                AFRAC( 1 ) = DYLNK * ( YY - YBEG )
                IF( AFRAC( 1 ) < 0.0 ) AFRAC( 1 ) = 0.0
            ENDIF

            IF( NCEL .GT. NDIM ) CALL REPORT_BAD_CELL

            !.....  Reset last  fraction if it ends   on interior of domain
            IF( YEND .LT. YORIG + YCELL * NROWS ) THEN
                YY = YORIG + YCELL * FLOAT( RRR - 1 )
                AFRAC( NCEL ) = DYLNK * ( YEND - YY )
                IF( AFRAC( NCEL ) < 0.0 ) AFRAC( NCEL ) = 0.0
            ENDIF

            RETURN

        END IF              !  if link crosses over domain at all

    END IF              !  if link crosses cells in Y direction only


    !.....  Case: (Near)-constant Y -- 1-D problem in X only:

    IF ( ( ROW .EQ. RRR ) .OR. ( ABS( YLNK ) .LT. 1.0E-5 * YCELL ) ) THEN

        !.....  Case:  Link is outside domain
        IF ( ROW .LT. 1  .OR. ROW .GT. NROWS ) THEN
            NCEL = 0
            RETURN
        END IF

        !.....  For cases in which XBEG is not located in a grid cell to the left
        !.....  of XEND, then flip the coordinates and use the same algorithm.

        IF ( COL .GT. CCC ) THEN

            ISAV = CCC
            CCC  = COL
            COL  = ISAV

            RSAV = XEND
            XEND = XBEG
            XBEG = RSAV

        ENDIF

        !.....  Case:  Link is outside domain
        IF ( CCC .LT. 1  .OR.  COL .GT. NCOLS ) THEN
            NCEL = 0
            RETURN

            !.....  Case: Link crosses domain at some point
        ELSE

            !.....  Initialize for all cells (inside domain) intersecting link
            J  = 0
            FF = DXLNK * XCELL
            JZ = NCOLS * ( ROW - 1 )
            DO  22  IC = MAX( 1,COL ) , MIN( CCC,NCOLS )
                J  = J + 1
                ACEL ( J ) = IC + JZ
                AFRAC( J ) = FF
22          CONTINUE
            NCEL = J

            IF( NCEL .GT. NDIM ) CALL REPORT_BAD_CELL

            !.....  Reset first fraction if it starts on interior of domain
            IF( XBEG .GT. XORIG ) THEN
                XX = XORIG + XCELL * FLOAT( COL )
                AFRAC( 1 ) = DXLNK * ( XX - XBEG )
                IF( AFRAC( 1 ) < 0.0 ) AFRAC( 1 ) = 0.0
            ENDIF

            !.....  Reset last  fraction if it ends   on interior of domain
            IF( XEND .LT. XORIG + XCELL * NCOLS ) THEN
                XX = XORIG + XCELL * FLOAT( CCC - 1 )
                AFRAC( NCEL ) = DXLNK * ( XEND - XX )
                IF( AFRAC( NCEL ) < 0.0 ) AFRAC( NCEL ) = 0.0
            ENDIF

            RETURN

        END IF          !  if link crosses over domain at all

    END IF              !  if link crosses cells in X direction only

    !...................
    !.....  2-D CASE:
    !.....  (by construction, we know that | XLNK |, | YLNK | > 1.0E5*CELL)
    !.....  Find the intersections of link with the grid, starting
    !.....  from XBEG, YBEG:
    !...................

    !...................
    !.....  Set fractions of link between vertical grid lines by
    !.....  using the X fractions
    !...................

    !.....  Establish loop bounds and starting and ending coordinates
    !.....  Also, precalculate pointers (XB  XE) for adjusting fractions
    !.....  based on orientation of link.
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

    J         = 0             !  cell counter
    XFAC( J ) = 9.0E36        !  sentinel on front end-of-list

    !.....  Initialize for all cells (inside domain) intersecting link
    FF = DXLNK * XCELL

    DO  IC = IX, NX, XINC

        !.....  Skip cell ID is outside domain
        IF( IC < 1 ) CYCLE

        !.....  First X-cell of link is inside domain
        IF( IC .EQ. IX .AND. START .GE. XORIG .AND. &
            IC .GE. 1  .AND. IC    .LE. NCOLS       ) THEN
            J         = J + 1
            XX        = XORIG + XCELL * FLOAT( XB )
            XCOL( J ) = IC
            XFAC( J ) = DXLNK * FLOAT( XINC ) * ( XX - START )

        !.....  Initialize first cell when link starts outside domain
        ELSEIF( IC .EQ. 1 ) THEN
            J         = J + 1
            XCOL( J ) = IC
            XFAC( J ) = FF

        !.....  Last X-cell of link is inside domain
        ELSEIF( IC   .EQ. NX .AND.                      &
                ENDS .LE. XORIG + XCELL * NCOLS .AND.   &
                IC   .GE. 1  .AND. IC .LE. NCOLS          ) THEN
            J         = J + 1
            XX        = XORIG + XCELL * FLOAT( XE )
            XCOL( J ) = IC
            XFAC( J ) = XFAC( J-1 ) + DXLNK * FLOAT( XINC ) * ( ENDS - XX )

        !.....  Set fractions for interior of domain and non-end of link
        ELSEIF( IC .GT. 1 .AND. IC .LE. NCOLS ) THEN
            J         = J + 1
            XCOL( J ) = IC
            XFAC( J ) = XFAC( J-1 ) + FF

        ENDIF

        IF( XFAC( J ) < 0.0 ) XFAC( J ) = 0.0

    END DO

    NX           = J             !  total number of columns intersected
    XFAC( NX+1 ) = 9.0E36        !  sentinel on tail end-of-list

    !.....  Case: Link is outside domain
    IF ( NX .EQ. 0 ) THEN
        NCEL = 0
        RETURN
    END IF


    !...................
    !.....  Set fractions of link between horizontal grid lines by
    !.....  using the Y fractions
    !...................

    !.....  Establish loop bounds and starting and ending coordinates
    !.....  Also, precalculate pointers (YB  YE) for adjusting fractions
    !.....  based on orientation of link.  Note that YINC is used for
    !.....  stepping through loop _AND_ for changing sign in the YFAC calc
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

    J         = 0             !  cell counter
    YFAC( J ) = 9.0E36        !  sentinel on front end-of-list

    !.....  Calculate fractions for cells (inside domain) intersecting link
    FF = DYLNK * YCELL

    DO  IR = IY, NY, YINC

        !.....  Skip cell ID is outside domain
        IF( IR < 1 ) CYCLE

        !.....  First Y-cell of link is inside domain
        IF( IR .EQ. IY .AND. START .GE. YORIG .AND.     &
            IR .GE. 1  .AND. IR    .LE. NROWS      ) THEN
            J         = J + 1
            YY        = YORIG + YCELL * FLOAT( YB )
            YROW( J ) = IR
            YFAC( J ) = DYLNK * FLOAT( YINC ) * ( YY - START )

        !.....  Initialize first cell when link starts outside domain
        ELSEIF( IR .EQ. 1 ) THEN
            J         = J + 1
            YROW( J ) = IR
            YFAC( J ) = FF

        !.....  Last Y-cell of link is inside domain
        ELSEIF( IR   .EQ. NY .AND.                      &
                ENDS .LE. YORIG + YCELL * NROWS .AND.   &
                IR   .GE. 1 .AND. IR .LE. NROWS           ) THEN
            J         = J + 1
            YY        = YORIG + YCELL * FLOAT( YE )
            YROW( J ) = IR
            YFAC( J ) = YFAC( J-1 ) + DYLNK * FLOAT( YINC ) * ( ENDS - YY )

        !.....  Set fractions for interior of domain and non-end of link
        ELSEIF( IR .GT. 1 .AND. IR .LE. NROWS ) THEN
            J         = J + 1
            YROW( J ) = IR
            YFAC( J ) = YFAC( J-1 ) + FF

        ENDIF

        IF( YFAC( J ) < 0.0 ) YFAC( J ) = 0.0

    END DO

    NY           = J             !  total number of columns intersected
    YFAC( NY+1 ) = 9.0E36        !  sentinel on tail end-of-list

    !.....  Case: Link is outside domain
    IF ( NY .EQ. 0 ) THEN
        NCEL = 0
        RETURN
    END IF

    !...................
    !....   Now merge the two intersection lists:
    !....   This algorithm starts at one end of the link,
    !....   and continues to the other.  The XFAC and YFAC variables
    !....   contain the total fraction of the link up to the point
    !....   where the link intersects column IX (for XFAC) and/or
    !....   row IY (for YFAC).   The next intersection produces the
    !....   next AFRAC for cell (IC,IR).  If the link passes through
    !....   a cell corner, then both the IX counter and IY counter
    !....   are incremented by 1.
    !...................

    !.....  Perform the merge.  Loop terminates when both lists hit
    !....  sentinels on either end of fractions arrays

    J  = 0
    FF = 0
    IX = 1
    IY = 1
    XINC = 1
    YINC = 1
133 CONTINUE

    !.....  Intersect new column and new row at same time
    !.....  or both arrays have reached end of lists
    IF( ABS( XFAC( IX ) - YFAC( IY ) ) .LT. 1.0E-05 ) THEN
        FAC = XFAC( IX ) - FF
        FF  = XFAC( IX )
        IC  = XCOL( IX )
        IR  = YROW( IY )
        IX  = IX + XINC
        IY  = IY + YINC

    !.....  Intersect new column next
    ELSE IF ( XFAC( IX ) .LT. YFAC( IY ) ) THEN
        FAC = XFAC( IX ) - FF
        FF  = XFAC( IX )
        IC  = XCOL( IX )
        IR  = YROW( IY )
        IX  = IX + XINC

    !.....  Intersect new row next
    ELSE
        FAC = YFAC( IY ) - FF
        FF  = YFAC( IY )
        IC  = XCOL( IX )
        IR  = YROW( IY )
        IY  = IY + YINC
    END IF
    IF( IR > 0 .AND. IC > 0 ) THEN
        IF ( FF .LE. 2.0 ) THEN       ! Check for sentinel
            J = J + 1
            ACEL ( J ) = IC +  NCOLS * ( IR - 1 )
            AFRAC( J ) = MIN( 1.0, MAX( 0.0, FAC ) )
            GO TO  133                  !  to head of loop
        END IF
    END IF

    !.......   Merge complete.  Return the number of cells found:

    NCEL = J
    RETURN

    ! ********************** INTERNAL SUBPROGRAMS ****************************

CONTAINS

    SUBROUTINE REPORT_BAD_CELL

        !.....  Local variables
        CHARACTER(300) MESG

        !......................

        WRITE( MESG,94010 ) 'INTERNAL ERROR: Bad number of cells ', &
               NCEL, 'computed in ' // PROGNAME
        CALL M3MSG2( MESG )
        EFLAG = .TRUE.

        !.......   Internal buffering formats.... 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

    END SUBROUTINE REPORT_BAD_CELL

END SUBROUTINE LNK2GRD



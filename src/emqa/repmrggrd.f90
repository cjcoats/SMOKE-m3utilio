
SUBROUTINE REPMRGGRD( RCNT, NX, IX, CX, EFLAG )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !       The REPMRGGRD routine assigns grid cells and gridding factors to each
    !       of the sources still selected by the program. If there is a subgrid,
    !       it will be used to subselect the grid cells.  The final OUTREC array
    !       is created by this routine before grouping, summing, and writing
    !       emissions data.  If the report contains normalization by grid cell
    !       these factors will be included in the gridding factors.
    !
    !  PRECONDITIONS REQUIRED:
    !       Routine is only called if gridding is being used
    !       Gridding matrix is allocated and populated
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 7/2000 by M Houyoux
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
    !***********************************************************************
    USE M3UTILIO

    !.......   MODULES for public variables
    !.......   This module is the inventory arrays
    USE MODSOURC, ONLY: INDEXA

    !.......  This module contains Smkreport-specific settings
    USE MODREPRT, ONLY: NMATX, LSUBGRID, ALLRPT, NSUBGRID, SUBGNAM, VALIDCEL

    !.......  This module contains report arrays for each output bin
    USE MODREPBN, ONLY: OUTSRC, OUTBIN, OUTCELL, OUTGFAC, NOUTREC, NSRCDROP

    !.......  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: NGRID, GDTYP, XCELL, YCELL, NCOLS, YORIG

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters
    INCLUDE 'CONST3.EXT'        !  physical and mathematical constants

    !.......   SUBROUTINE ARGUMENTS
    INTEGER, INTENT (IN) :: RCNT            ! current report number
    INTEGER, INTENT (IN) :: NX( NGRID )     ! no. srcs per cell
    INTEGER, INTENT (IN) :: IX( NMATX )     ! src IDs
    REAL   , INTENT (IN) :: CX( NMATX )     ! gridding coefficients
    LOGICAL, INTENT(OUT) :: EFLAG           ! true: error found

    !.......  Local allocatable arrays
    REAL, ALLOCATABLE, SAVE :: NORMFAC( : )

    !.......  Local variables
    INTEGER         C              ! indices and counters
    INTEGER         IOS            ! i/o status
    INTEGER         ROW            ! tmp row number

    REAL            FAC            ! temporary factor
    REAL            Y              ! tmp Y cell center
    REAL            Y0             ! y origin less 0.5*dy

    LOGICAL      :: FIRSTFLAG = .TRUE.      ! True: no normalize by cell area yet found

    CHARACTER(300)  MESG           ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'REPMRGGRD'     ! program name

    !***********************************************************************
    !   begin body of subroutine REPMRGGRD

    !.......  For first time a normalized grid is encountered by this routine,
    !           compute the cell-area normalization factors for full grid.
    IF ( ALLRPT( RCNT )%NORMCELL .AND. FIRSTFLAG ) THEN

        ALLOCATE( NORMFAC( NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NORMFAC', PROGNAME )

        SELECT CASE( GDTYP )
          CASE( LAMGRD3, UTMGRD3, MERGRD3, STEGRD3, POLGRD3, EQMGRD3, TRMGRD3, ALBGRD3, LEQGRD3 )
            FAC   = 1. / ( XCELL * YCELL )
            NORMFAC = FAC       ! array

          CASE( LATGRD3 )

            Y0 = YORIG - 0.5 * YCELL
            FAC = ( PI * REARTH / 180. )**2 * XCELL * YCELL
            DO C = 1, NGRID
                ROW = 1 + INT( (C-1) / NCOLS )
                Y = Y0 + REAL( ROW ) * YCELL
                NORMFAC( C ) = 1. / ( COS( Y ) * FAC )
            END DO

          CASE DEFAULT

            WRITE( MESG,94010 ) 'INTERNAL ERROR: Grid type', GDTYP, &
                   'not recognized in ' // PROGNAME
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        END SELECT

        FIRSTFLAG = .FALSE.

    END IF

    !.......  Set report-specific local settings
    LSUBGRID = ( ALLRPT( RCNT )%SUBGNAM .NE. ' ' )

    !.......  Deallocate source list and bins before reallocating for source-cell
    !           intersections.  Valid sources are identified with INDEXA     != 0.
    IF( ALLOCATED( OUTSRC ) ) DEALLOCATE( OUTSRC, OUTBIN )

    !.......  Deallocate cell and gridding information if they have been
    !           allocated for a previous report
    IF( ALLOCATED( OUTCELL ) ) DEALLOCATE( OUTCELL, OUTGFAC )

    !.......  Count cell-source intersections...
    CALL SOURCE_CELL_X( 'COUNT', NX, IX, CX )

    !.......  Warning if no source-cell intersections
    IF( NOUTREC .EQ. 0 ) THEN
        EFLAG = .TRUE.
        MESG = BLANK5 // 'ERROR: No source-cell intersections '//&
               'for report.  Output will be empty.'
        CALL M3MSG2( MESG )

    END IF

    !.......  Allocate memory for output record arrays
    ALLOCATE( OUTSRC( NOUTREC ),    &
              OUTBIN( NOUTREC ),    &
             OUTCELL( NOUTREC ),    &
             OUTGFAC( NOUTREC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'OUTSRC...OUTGFAC', PROGNAME )

    !.......  Store sources, cells, and gridding factors
    CALL SOURCE_CELL_X( 'STORE', NX, IX, CX )

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I10, :, 1X ) )

    !******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

    !.......  This internal function counts or computes the source-cell
    !               intersections for a whole grid or subgrid, and for specifically-
    !               selected sources from the inventory.
    SUBROUTINE SOURCE_CELL_X( COMMAND, NX, IX, CX )

        !.......  Subprogram arguments (note - array dimensions from MODREPRT)
        CHARACTER(*), INTENT (IN) :: COMMAND             ! "COUNT" or "STORE"
        INTEGER     , INTENT (IN) :: NX( NGRID )         ! no. srcs per cell
        INTEGER     , INTENT (IN) :: IX( NMATX )         ! src IDs
        REAL        , INTENT (IN) :: CX( NMATX )         ! gridding coefficients

        !.......  Local variables
        INTEGER   C, IC, IG, J, K, N, S

        !----------------------------------------------------------------------

        !.......  Special case of no subgrid and all sources selected
        !.......  Have this case explicitly because it will not be slowed down
        !               by uncessary counting and conditionals.
        IF( NSRCDROP .EQ. 0 .AND. .NOT. LSUBGRID ) THEN

            NOUTREC = NMATX

            IF( COMMAND .EQ. 'COUNT' ) RETURN

            K = 0
            DO C = 1, NGRID

                DO J = 1, NX( C )

                    K = K + 1
                    OUTCELL( K ) = C
                    OUTSRC ( K ) = IX( K )
                    OUTGFAC( K ) = CX( K )

                END DO

            END DO

        !.......  Special case of no subgrid and partial sources
        ELSE IF ( .NOT. LSUBGRID ) THEN

            K = 0
            N = 0
            DO C = 1, NGRID

                DO J = 1, NX( C )

                    K = K + 1
                    S = IX( K )

                    !.......  Skip source if it has not been selected
                    IF( INDEXA( S ) .LE. 0 ) CYCLE

                    !.......  Increment count of output records
                    N = N + 1

                    !.......  If only counting, then skip the rest of the loop
                    IF( COMMAND .EQ. 'COUNT' ) CYCLE

                    !.......  Store the factors for the current cell and source
                    OUTCELL( N ) = C
                    OUTSRC ( N ) = S
                    OUTGFAC( N ) = CX( K )

                END DO

            END DO

            NOUTREC = N             ! Store output count of source-cell interstns

        !.......  Process for a subgrid and/or partial source list
        ELSE

            !.......  Determine subgrid index based on name (the validity has
            !                   already been confirmed elsewhere).
            IG = INDEX1( ALLRPT( RCNT )%SUBGNAM, NSUBGRID, SUBGNAM )

            K  = 0
            N  = 0
            IC = 1
            DO C = 1, NGRID

                !.......  If current cell is valid for the subgrid...
                IF( C .EQ. VALIDCEL( IC,IG ) ) THEN

                    !......  Loop through sources for current cell
                    DO J = 1, NX( C )

                        K = K + 1
                        S = IX( K )

                        !.......  Skip source if it has not been selected
                        IF( INDEXA( S ) .LE. 0 ) CYCLE

                        !.......  Increment count of output records
                        N = N + 1

                        !.......  If only counting, then skip the rest of the loop
                        IF( COMMAND .EQ. 'COUNT' ) CYCLE

                        !.......  Store the factors for the current cell and src
                        OUTCELL( N ) = C
                        OUTSRC ( N ) = S
                        OUTGFAC( N ) = CX( K )

                    END DO                  ! End loop over sources for current cell

                    IC = IC + 1

                !.......  Otherwise, step through gridding matrix
                ELSE

                    K = K + NX( C )

                END IF              ! Check if valid cell

            END DO                  ! End loop over cells

            NOUTREC = N             ! Store output count of source-cell interstns

        END IF                      ! If ( all sources and full grid ) or not

        !.......  If needed, update all gridding factors with division by cell
        !               area.
        IF ( ALLRPT( RCNT )%NORMCELL ) THEN

            DO N = 1, NOUTREC
                C = OUTCELL( N )
                OUTGFAC( N ) = OUTGFAC( N ) * NORMFAC( C )
            ENDDO

        END IF

        RETURN

    END SUBROUTINE SOURCE_CELL_X

END SUBROUTINE REPMRGGRD



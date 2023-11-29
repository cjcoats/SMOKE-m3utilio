
        PROGRAM SURGTOOL

C***********************************************************************
C  program body starts at line
C
C  DESCRIPTION:
C      Inputs SMOKE-formatted gridding surrogates
C      for a "fine" input grid and a grid definition for a "coarse" output
C      grid.  It produces approximate "coarse" grid surrogates file, the
C      accuracy of which depends on how fine the resolution of the input grid
C      is relative to that of the output grid.
C
C  PRECONDITIONS REQUIRED:
C      - Input & output surrogates on Lat-Lon, Lambert Conformal, or UTM map
C        projections (can perform Lambert-to-Lambert and UTM zone-to-zone
C        transformations).
C      - Input grid resolution is much finer than output grid resolution
C      - Supports SMOKE-formatted surrogate coefficient files only.
C      - Input logical file names are defined.
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Prototype  4/1998 by .Carlie J. Coats, Jr as program CROSSGRID
c       Version    4/1999 by CJC to support both EPS2 and EMS85 files
C       Version    10/00 by JMV:  source appropriated without permission as SURGTOOL
C       Version    1/02 : Updated for Models-3 surrogate format and included
C                     in SMOKE emutil module
C       Version 11/2023 by CJC:  USE M3UTILIO, MODGCTP, and related changes
C
C***********************************************************************
C
C Program copyright (C) 1998-2002 MCNC and (C) 2--3-2004 Baron Advanced
C Meteorological Systems, LLC, and released under Version 2 of the
C GNU General Public License.
C See http://www.gnu.org/copyleft/gpl.html
C
C*************************************************************************
        USE M3UTILIO
        USE MODGCTP

C...........   MODULES for public variables
C...........   This module contains the gridding surrogates tables
        USE MODSURG, ONLY: NSRGREC, IDXSRGA, SCELLA, SFIPSA, SSRGIDA,
     &                     SFRACA

C.........  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: GDTYP, P_ALP, P_BET, P_GAM, YCENT,
     &                     NCOLS, NROWS, XORIG, YORIG, XCELL, YCELL,
     &                     XCENT, YCENT, GRDNM

       IMPLICIT NONE

C...........   EXTERNAL FUNCTIONS and their descriptions:

        LOGICAL, EXTERNAL :: DSCM3GRD

C...........   LOCAL PARAMETERS
        CHARACTER(16), PARAMETER :: PROGNAME = 'SURGTOOL'  !  program name
        CHARACTER(50), PARAMETER ::
     &  CVSW = '$Name SMOKEv5.0_Jun2023$' ! CVS release tag

C............   Allocatable arrays for computing grid cell intersections
        INTEGER, ALLOCATABLE::   C2( :,: )     !  output grid cell
        REAL*8,  ALLOCATABLE::   XX( :,: )     !  output grid X
        REAL*8,  ALLOCATABLE::   YY( :,: )     !  output grid Y

C...........   Allocatable arrays for unsorted output surrogates
        INTEGER                 OUTNSRG       !  number of output entries
        INTEGER, ALLOCATABLE :: INDXA( : )    !  sorting index
        INTEGER, ALLOCATABLE :: CELLA( : )    !  cell number
        INTEGER, ALLOCATABLE :: FIPSA( : )    !  region code
        INTEGER, ALLOCATABLE :: SSCSA( : )    !  surrogate ID
        REAL   , ALLOCATABLE :: FRACA( : )    !  surrogate fraction

C...........   Logical file names and unit numbers

        INTEGER         LDEV    !  log-device
        INTEGER         SDEV    !  for  input surrogate coeff file
        INTEGER         XDEV    !  for output surrogate coeff  file

C...........    Output grid variables
        INTEGER      NCOLS2      ! number of grid columns
        INTEGER      NROWS2      ! number of grid rows
        INTEGER      NTHIK2      ! BOUNDARY:  perim thickness (cells)
        INTEGER      GDTYP2      ! grid type:  1=LAT-LON, 2=UTM, ...
        REAL(8)      P_ALP2      ! first, second, third map
        REAL(8)      P_BET2      ! projection descriptive
        REAL(8)      P_GAM2      ! parameters.
        REAL(8)      XCENT2      ! lon for coord-system X=0
        REAL(8)      YCENT2      ! lat for coord-system Y=0
        REAL(8)      XORIG2      ! X-coordinate origin of grid (map units)
        REAL(8)      YORIG2      ! Y-coordinate origin of grid
        REAL(8)      XCELL2      ! X-coordinate cell dimension
        REAL(8)      YCELL2      ! Y-coordinate cell dimension
        CHARACTER(NAMLEN3) :: GDNAM2            ! output grid name
        CHARACTER(NAMLEN3) :: COORUNIT2         !  coord sys projection units
        CHARACTER(NAMLEN3) :: COORD2            ! coord system name
        CHARACTER(MXDLEN3) :: GDESC2            ! output grid desc (if any)

C...........   Other local variables
        INTEGER          C, S, F, J, K, L, M, N, R    !  counters, subscripts
        INTEGER          L1, L2                 !  tmp lengths
        INTEGER          COL, ROW               !  tmp row and column
        INTEGER          IOS                    !  I/O status
        INTEGER          IREC                   !  input line (record) number
        INTEGER          LCEL                   !  previous cell number
        INTEGER          LFIP                   !  previous region code
        INTEGER          LSSC                   !  previous surrogate code

        REAL*8           DDX, DDY               !  inverse of cell widths
        REAL             FRAC                   !  tmp output surrogate fraction
        REAL*8           XZ, YZ

        LOGICAL       :: EFLAG = .FALSE.        !  true: error found

        CHARACTER(16)    OUTTYPE
        CHARACTER(16)    OUTUNIT
        CHARACTER(16)    SRGFMT                 !  surrogates format
        CHARACTER(256)   MESG                   !  message buffer

C***********************************************************************
C   begin body of program SURGTOOL

        LDEV = INIT3()

C.........  Write out copyright, version, web address, header info, and prompt
C           to continue running the program.

        WRITE( *, '(5X, A)' )
     &' ',
     &'Program CROSSGRID to take the SMOKE 2 format surrogate',
     &'coefficient file for a "fine" grid, a grid definition for',
     &'that "coarse" grid, and produce the approximate "coarse" ',
     &'surrogate coefficient file. ',
     &' ',
     &'Program copyright (C) 1998-2002 MCN, (C) 2003-2004 Baron',
     &'Advanced Meteorological Systems, LLC, and ',
     &'(C) 2023 Carlie J. Coats, Jr. ',
     &'Released under Version 2 of the GNU General Public License.',
     &'See URL',
     &'    http://www.gnu.org/copyleft/gpl.html',
     &' ',
     &'Comments and questions are welcome and can be sent to',
     &'    carlie@jyarborough.com.',
     &' '
        WRITE( *, '(5X, A)' )
     &'NOTE:  Current version is for Lat-Lon, Lambert, and UTM',
     &'map projections only.',
     &' '
        WRITE( *, '(5X, A)' )
     &'USAGE:',
     &'You will need to enter GRIDDESC grid names for the input and',
     &'output grids, and the logical names for the input and output',
     &'files (and to have set them prior to program launch, using ',
     &' ',
     &'    "setenv <logicalname> <pathname>").  ',
     &' ',
     &'Input files must have been sorted according to SMOKE ',
     &'conventions, prior to program execution.',
     &' ',
     &'You may use END_OF-FILE (control-D) to quit the program',
     &'during logical-name entry. Default responses are indicated',
     &'in brackets [LIKE THIS].',
     &' ',
     &'Program version $Id$',
     &  ' '

C.........  Open input surrogates file
        SDEV = PROMPTFFILE(
     &         'Enter logical name for input GRIDDING SURROGATE file',
     &         .TRUE., .TRUE., 'INFILE', PROGNAME )

C.........  Read the surrogates header and initialize the grid description
C.........  Also get the format of the surrogates file.
        CALL RDSRGHDR( .FALSE., SDEV, SRGFMT )

        MESG = 'NOTE: Input surrogates are ' // TRIM( SRGFMT ) //
     &         ' format.'
        CALL M3MSG2( MESG )

C.........  Initialize gridding module
        CALL CHKGRID( GDNAM3D, 'SURROGATES', 1, EFLAG )

C.........  Get grid parameters from grid description file for output grid.
        IF ( .NOT. DSCM3GRD( GDNAM2, GDESC2, COORD2, GDTYP2, COORUNIT2,
     &                       P_ALP2, P_BET2, P_GAM2, XCENT2,
     &                       YCENT2, XORIG2, YORIG2, XCELL2,
     &                       YCELL2, NCOLS2, NROWS2, NTHIK2 ) ) THEN

            MESG = 'Could not get Models-3 grid description.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C........  Compute X-Y coordinates of fine grid relative to coarse grid
C........  and then mapping C2 from fine cells to coarse cells:

        ALLOCATE( C2( NCOLS, NROWS ),
     &            XX( NCOLS, NROWS ),
     &            YY( NCOLS, NROWS ), STAT = IOS )
        CALL CHECKMEM( IOS, 'C2,XX,YY', PROGNAME )

        CALL GRID2XY( GDTYP2, P_ALP2, P_BET2, P_GAM2, XCENT2, YCENT2,
     &                GDTYP,  P_ALP,  P_BET,  P_GAM,  XCENT,  YCENT,
     &                NCOLS,  NROWS,  XORIG,  YORIG,  XCELL,  YCELL,
     &                XX, YY )

        DDX = 1.0d0 / XCELL2
        DDY = 1.0d0 / YCELL2

        DO  R = 1, NROWS       !  loop through cells of input grid
        DO  C = 1, NCOLS

            XZ = XX( C,R ) - XORIG2
            YZ = YY( C,R ) - YORIG2
            COL = 1 + INT( DDX * XZ )
            ROW = 1 + INT( DDY * YZ )

C.............  If center of input cell (col,row) is outside output grid,
C               then store missing intersection.
            IF ( XZ  .LT. 0.0d0   .OR.  YZ  .LT. 0.0d0  .OR.
     &           COL .GT. NCOLS2  .OR.  ROW .GT. NROWS2 ) THEN
                C2( C,R ) = IMISS3

C.............  Store which cell in output grid for current input cell (col,row)
            ELSE
                C2( C,R ) = COL + NCOLS2 * ( ROW - 1 )

            END IF

        END DO
        END DO          !  end looping through input grid

C.........  Open output file
        XDEV = PROMPTFFILE(
     &         'Enter logical name for output GRIDDING SURROGATE file',
     &           .FALSE., .TRUE., 'OUTFILE', PROGNAME )

        CALL M3MSG2( 'Reading gridding surrogates...' )

C.........  Allocate memory for and read the gridding surrogates file,
C           extracting data for a subgrid, if necessary
        CALL RDSRG( .FALSE., SDEV, SRGFMT, NROWS, NCOLS )

C.........  Allocate memory for output surrogates assuming that number of
C           records from input file will > output records.
        ALLOCATE( INDXA( NSRGREC ),
     &            CELLA( NSRGREC ),
     &            FIPSA( NSRGREC ),
     &            SSCSA( NSRGREC ),
     &            FRACA( NSRGREC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXA...FRACA', PROGNAME )

        CALL M3MSG2( 'Regridding surrogates...' )

C.........  Regrid the surrogates based on previously computed factors
        OUTNSRG = 0
        DO J = 1, NSRGREC

            K = IDXSRGA( J )
            C = SCELLA ( K )

            ROW = 1 + INT( (C-1) / NCOLS )
            COL = C - NCOLS * ( ROW-1 )

C.............  Do not store output cells if input cell does not intersect
C               any output cells

            IF ( C2( COL,ROW ) .NE. IMISS3 ) THEN

                OUTNSRG = OUTNSRG + 1

                IF ( OUTNSRG .GT. NSRGREC ) CYCLE

                INDXA( OUTNSRG ) = OUTNSRG
                CELLA( OUTNSRG ) = C2( COL,ROW )
                FIPSA( OUTNSRG ) = STR2INT( SFIPSA ( K ) )
                SSCSA( OUTNSRG ) = SSRGIDA( K )
                FRACA( OUTNSRG ) = SFRACA ( K )

            END IF

        END DO

C...........   Give error if couldn't store all of the output intersections
        IF( OUTNSRG .GT. NSRGREC ) THEN

            MESG = 'INTERNAL ERROR: Assumption about memory ' //
     &             'allocation for output surrogates is false.'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C...........   Sort and process the surrogates:
        CALL M3MSG2( 'Sorting output gridding surrogates...' )

        CALL SORTI3( OUTNSRG, INDXA, SSCSA, CELLA, FIPSA )

        CALL M3MSG2 ( 'Creating output gridding surrogate file...' )

C.........  Write header line
        MESG = '#GRID ' // GDNAM2
        WRITE( XDEV, 93980 ) TRIM( MESG ), XORIG2, YORIG2, XCELL2,
     &         YCELL2, NCOLS2, NROWS2, NTHIK2, TRIM( OUTTYPE ),
     &         TRIM( OUTUNIT ), P_ALP2, P_BET2, P_GAM2, XCENT2, YCENT2

C.........  Loop through output entries and write them
        LCEL = IMISS3
        LFIP = IMISS3
        LSSC = IMISS3
        M    = 0
        N    = 0
        FRAC = 0.0

        DO  J = 1, OUTNSRG

            K = INDXA ( J )
            C = CELLA ( K )
            F = FIPSA ( K )
            S = SSCSA ( K )

C.............  If same FIPS code, cell, and surrogate, sum surrogate fraction
            IF ( F .EQ. LFIP  .AND.
     &           C .EQ. LCEL  .AND.
     &           S .EQ. LSSC ) THEN

                FRAC = FRAC + FRACA( K )
                N    = N + 1

C.............  Otherwise...
            ELSE

                IF ( N .GT. 0 ) THEN

                    ROW = 1 + INT( (LCEL-1) / NCOLS2 )
                    COL = LCEL - NCOLS2 * ( ROW - 1 )

                    WRITE( XDEV, 93020, IOSTAT=IOS )
     &                  LSSC, LFIP, COL, ROW, FRAC

                    IF ( IOS .GT. 0 ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG, 94010 ) 'ERROR: I/O error', IOS,
     &                         'writing output file at line', M
                    END IF

                END IF

                M = M + 1
                N = 1
                LCEL = C
                LFIP = F

                LSSC = S
                FRAC = FRACA( K )

            END IF            !  if new record encountered

        END DO

C.........  Write the last entry in the output file
        IF ( N .GT. 0 ) THEN
            ROW = 1 + INT( (LCEL-1) / NCOLS2 )
            COL = LCEL - NCOLS2 * ( ROW-1 )

            FRAC = FRAC / FLOAT( N )

            WRITE( XDEV, 93020, IOSTAT=IOS )
     &              LSSC, F, COL, ROW, FRAC

C.............  Write error message, if needed
            IF ( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: I/O error', IOS,
     &                 'writing output file at line', M
            END IF

        END IF

C.........  Abort if error
        IF( EFLAG ) THEN

            MESG = 'Problem writing output file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

C..........   Normal completion
        MESG = ' '
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 0 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

93020   FORMAT( I3, 1X, I6.5, I8, I8, 1X, F11.8 )

93980   FORMAT( A, 1X, 4( E15.8, 1X ), 3( I8, 1X ), 2( A, 1X ),
     &          5( E15.8, 1X ) )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I10, :, 2X ) )

        END PROGRAM SURGTOOL

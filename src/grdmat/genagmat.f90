
SUBROUTINE GENAGMAT( GNAME, FDEV, MXSCEL, NASRC, NMATX, VFLAG,          &
                     DEFSRGID, FSGFLAG, NX, IX, CX, NCOEF, CMAX, CMIN )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created by M. Houyoux 5/99
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !************************************************************************
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
    !************************************************************************
    USE M3UTILIO

    !.......   MODULES for public variables
    !.......   This module is the source inventory arrays
    USE MODSOURC, ONLY: XLOCA, YLOCA, CIFIP, CELLID, CSOURC

    !.......   This module contains the gridding surrogates tables
    USE MODSURG, ONLY: NCELLS, FIPCELL, NSRGS, SRGLIST, NSRGFIPS,           &
                       SRGFIPS, NTSRGDSC, SRGFNAM, SRGFCOD, SRGFMT,         &
                       SRGNCOLS, SRGNROWS, NTLINES, SRGCSUM, SRGFRAC

    !.......  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: NGRID, NCOLS, NROWS, XOFF, YOFF

    !.......  This module contains the lists of unique source characteristics
    USE MODLISTS, ONLY: NINVIFIP

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: NCHARS, NSRC

    !.......   This module contains the cross-reference tables
    USE MODXREF, ONLY: ASRGID

    USE MODGRDLIB

    IMPLICIT NONE

    !.......   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: GNAME             ! gridding mtx logical name
    INTEGER     , INTENT (IN) :: FDEV              ! surg codes report file
    INTEGER     , INTENT (IN) :: MXSCEL            ! max sources per cell
    INTEGER     , INTENT (IN) :: NASRC             ! no. sources
    INTEGER     , INTENT (IN) :: NMATX             ! no. source-cell intersects
    LOGICAL     , INTENT (IN) :: VFLAG             ! true: using variable grid
    INTEGER     , INTENT (IN) :: DEFSRGID          ! default surrogate ID
    LOGICAL     , INTENT (IN) :: FSGFLAG           ! true: using default fallback surrogate
    INTEGER     , INTENT(OUT) :: NX  ( NGRID )     ! no. srcs per cell
    INTEGER     , INTENT(OUT) :: IX  ( NMATX )     ! src IDs
    REAL        , INTENT(OUT) :: CX  ( NMATX )     ! gridding coefficients
    INTEGER     , INTENT(OUT) :: NCOEF             ! no. of gridding coeffs
    INTEGER     , INTENT(OUT) :: CMAX              ! max no. of sources per cell
    INTEGER     , INTENT(OUT) :: CMIN              ! min no. of sources per cell

    !.......   EXTERNAL FUNCTIONS
    LOGICAL, EXTERNAL :: DSCM3GRD
    INTEGER, EXTERNAL :: GETFLINE
    LOGICAL, EXTERNAL :: BLKORCMT

    !.......   Local allocatable arrays...
    !.......   LOCAL VARIABLES and their descriptions:
    !.......   Local parameters
    INTEGER,       PARAMETER :: MXSEG    = 10               ! # of potential line segments
    CHARACTER(16), PARAMETER :: PROGNAME = 'GENAGMAT'     ! program name

    !.......   Other arrays
    CHARACTER(20) SEGMENT( MXSEG )                 ! Segments of parsed lines
    CHARACTER(60), ALLOCATABLE :: TMPLINE( : )     ! tmp line buffer

    !.......   Scratch Gridding Matrix (subscripted by source-within-cell, cell)
    INTEGER         IS ( MXSCEL, NGRID )     ! source IDs for each cell
    REAL            CS ( MXSCEL, NGRID )     ! factors

    !.......   Temporary array for flagging sources that are outside the
    !         domain and for flagging sources with all zero surrogate fractions
    LOGICAL             INDOMAIN( NASRC )     ! true: src is in the domain
    CHARACTER(FIPLEN3)  FIPNOSRG( NINVIFIP  )     ! cy/st/co codes w/o surrogates

    !.......   Temporary arrays for storing surrogate codes to use
    INTEGER         SURGID1( NSRC )     ! primary surrogate code
    INTEGER         SURGID2( NSRC )     ! secondary surrogate code

    !.......   Other local variables
    INTEGER         C, F, I, II, J, JJ, K, KK, N, NT, S     !  indices and counters.

    INTEGER         GDEV          !  for surrogate coeff file
    INTEGER         COL           ! tmp column
    INTEGER         TCOL          ! tmp column
    INTEGER         ID1,ID2       ! tmp primary and secondary surg codes
    INTEGER         IOS           ! i/o status
    INTEGER         ISIDX         ! tmp surrogate ID code index
    INTEGER         ISDEF         ! default surrogate ID code index
    INTEGER         JMAX          ! counter for storing correct max dimensions
    INTEGER         L2            ! string length
    INTEGER         NCEL          ! tmp number of cells
    INTEGER         NNOSRG        ! no. of cy/st/co codes with no surrogates
    INTEGER         ROW           ! tmp row
    INTEGER         TROW          ! tmp row
    INTEGER      :: NTL = 0       ! max no. of line buffers
    INTEGER         TGTSRG        ! target surrogates code
    INTEGER         SSC           ! surrogates code
    INTEGER      :: NLINES = 0    ! number of lines in input file

    REAL            FRAC        ! tmp surrogate fraction
    REAL*8          XX, YY

    LOGICAL      :: EFLAG = .FALSE.        !  true: error detected
    LOGICAL      :: LFLAG = .FALSE.        !  true: location data available
    LOGICAL      :: XYSET = .FALSE.        ! true: X/Y available for src
    LOGICAL      :: LASTIME = .FALSE.      ! true: X/Y available for src
    LOGICAL      :: CFLAG   = .FALSE.      ! true: called by sizgmat, false: called by gen[a|m]gmat
    LOGICAL      :: WFLAG   = .FALSE.      ! true: per iteration warning flag

    CHARACTER(FIPLEN3) CFIP       ! tmp country/state/county code
    CHARACTER(FIPLEN3) LFIP       ! cy/st/co code from previous iteration
    CHARACTER(200)  LINE          ! Read buffer for a line
    CHARACTER(16)   COORUNIT      !  coordinate system projection units
    CHARACTER(80)   GDESC         !  grid description
    CHARACTER(256)  MESG          !  message buffer
    CHARACTER(256)  NAMBUF        !  surrogate file name buffer
    CHARACTER(256)  NAMBUFT       !  tmp surrogate file name buffer
    CHARACTER(256)  TSRGFNAM      !  tmp surrogate file name buffer

    CHARACTER(SRCLEN3)    CSRC      ! tmp source chars string

    !***********************************************************************
    !   begin body of subroutine GENAGMAT

    !.......  Initialize number of sources per cell counter
    NX = 0       ! Array

    SURGID1 = 0       ! array
    SURGID2 = 0       ! array
    LFIP  = ' '
    NNOSRG   = 0
    JMAX  = -1

    !.......  Store the number and values of unfound cy/st/co codes
    !.......  Keep track of sources that are outside the domain
    DO I = 1, NSRC

        CFIP  = CIFIP( I )

        F = FINDC( CFIP, NSRGFIPS, SRGFIPS )

        IF ( F .LE. 0 ) THEN

            IF( CFIP .NE. LFIP ) THEN
                NNOSRG = NNOSRG + 1
                FIPNOSRG( NNOSRG ) = CFIP
                LFIP = CFIP
            END IF

            INDOMAIN( I ) = .FALSE.
            CYCLE       ! To next source

        END IF

    END DO

    !.......  Set flag to indicate that XLOCA/YLOCA are available
    LFLAG = ALLOCATED( XLOCA )

    MESG = 'Computing gridding matrix and statistics...'
    CALL M3MSG2( MESG )

    !.......  default fallback surrogate will run at last after
    !         re-assigned zero fraction surrogate

    DO II = 1, NSRGS      ! loop through only the surrogate code assigned by sources

        TGTSRG = SRGLIST( II )
        ISDEF  = FIND1( DEFSRGID, NSRGS, SRGLIST )
        KK = II

        IF( FSGFLAG ) THEN

            !.......  default fallback surrogate will run at last after
            !         re-assigned zero fraction surrogate
            IF( II >= ISDEF ) THEN
                IF( II == NSRGS ) THEN
                    TGTSRG = DEFSRGID
                    KK = ISDEF
                ELSE
                    TGTSRG = SRGLIST( II + 1 )
                    KK = II + 1
                END IF
            END IF

        END IF

        NTL = NTLINES( KK )

        !.......  Warning message when there are no surrogate available due to out of domain
        IF( NTL .EQ. 0 ) THEN
            WRITE( MESG,94010 ) 'WARNING: The surrogate', TGTSRG,           &
                ' does not exist within a domain range: emissions'//        &
                ' of sources assigned to', TGTSRG, ' will be set to zero'
            CALL M3MSG2( MESG )

            SRGFRAC = 0.0
            SRGCSUM = 0.0

        ELSE
            WRITE( MESG,94010 ) 'Looping surrogate', TGTSRG,        &
                ' through sources to generate gridding matrix'
            CALL M3MSG2( MESG )

            !......  Allocate memory for indices to surrogates tables for each source
            ALLOCATE( TMPLINE( NTL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TMPLINE', PROGNAME )

            !......  If surrogates are needed, read and store the gridding surrogates,
            !        allocate memory for the surrogate assignments, and assign
            !        surrogates to each source.
            NT = 0
            TMPLINE = ' '
            TSRGFNAM = ' '

            DO I = 1, NTSRGDSC              ! Open all surrogate files using the same srg code

                !.......  Prompt for and open I/O API output file(s)...
                CALL GETENV( 'SRGPRO_PATH', NAMBUF )
                NAMBUFT = TRIM( NAMBUF ) // SRGFNAM( I )

                IF( TGTSRG .NE. SRGFCOD( I ) ) CYCLE

                IF( NAMBUFT .NE. TSRGFNAM  ) THEN
                    CALL OPEN_SRGFILE

                    !.......  Reading surrogate files
                    DO JJ = 1, NLINES

                        READ ( GDEV, 93000, END=111, IOSTAT=IOS )LINE

                        IF ( BLKORCMT( LINE ) ) CYCLE

                        CALL PARSLINE( LINE, MXSEG, SEGMENT )
                        SSC    = STR2INT ( SEGMENT( 1 ) )
                        TCOL   = STR2INT ( SEGMENT( 3 ) )
                        TROW   = STR2INT ( SEGMENT( 4 ) )

                        !......  Check the value of the column number
                        IF( TCOL .LT. 0 .OR.  TCOL .GT. SRGNCOLS .OR.    &
                          ( TROW .EQ. 0 .AND. TCOL .NE. 0 ) ) THEN
                            WFLAG = .TRUE.
                        END IF

                        !......  Check the value of the row number
                        IF( TROW .LT. 0 .OR.  TROW .GT. SRGNROWS .OR.    &
                          ( TCOL .EQ. 0 .AND. TROW .NE. 0 ) ) THEN
                            WFLAG = .TRUE.
                            CALL M3MESG( MESG )

                        !......  Special treatment for cell (0,0) (skip for now)
                        ELSE IF( TROW .EQ. 0 .AND. TCOL .EQ. 0 ) THEN
                            CYCLE

                        END IF

                        !......  Adjust column and row for subgrid
                        TCOL = TCOL - XOFF
                        TROW = TROW - YOFF

                        !......  Skip entry after subgrid adjustment
                        IF( TCOL .LE. 0 .OR. TCOL .GT. NCOLS .OR.    &
                            TROW .LE. 0 .OR. TROW .GT. NROWS ) CYCLE

                        !......  Skip entry if rows and columns are out of range
                        IF( WFLAG ) CYCLE

                        !......  Skip entry if SSC is not in the assigned SRGLIST by source
                        IF( SSC .EQ. TGTSRG ) THEN
                            NT = NT + 1
                            TMPLINE( NT ) = LINE
                        END IF

111                 END DO

                    TSRGFNAM = NAMBUFT                ! store a previous surrogate file name buffer
                    CLOSE( GDEV )

                END IF                  ! skip if surrogate file has the same srg file

            END DO                   ! loop over all surrogate files in SRGDESC file

            !.......  Populating assigned surrogated
            CALL RDSRG4GRD( NT, TMPLINE, CFLAG )             ! populating surrogates

            DEALLOCATE( TMPLINE )

        END IF       ! end of populating surrogates

    !.......   Compute gridding matrix:
    !.......       First case:   explicit link (ILINK > 0
    !.......       Second case:  some LNKDEF entry applies
    !.......       Third case:   FIP+roadtype cross-reference match
    !.......       Fourth case:  state+roadtype cross-reference match
    !.......       fifth case:   roadtype cross-reference match
    !.......       sixth case:   fallback default

        DO S = 1, NSRC

            CFIP = CIFIP ( S )
            C    = CELLID( S )
            CSRC = CSOURC( S )
            SSC  = ASRGID( S )

            IF( SSC .NE. TGTSRG ) CYCLE

            !.......  Determine if x/y location is available
            XYSET = .FALSE.
            IF( LFLAG ) XYSET = ( XLOCA( S ) .GT. AMISS3 )

            !.......  Special case for source has an x/y location
            IF( XYSET ) THEN

                !.......  If source is in the domain, get cell number and store
                XX = XLOCA( S )
                YY = YLOCA( S )
                IF( INGRID( XX, YY, NCOLS, NROWS, COL, ROW  ) ) THEN

                    C = ( ROW-1 ) * NCOLS + COL
                    J = NX( C )

                    !.......  Check that the maximum number of sources per cell is ok
                    IF ( J .LT. MXSCEL ) THEN                      ! Not .LE.
                        J = J + 1
                        IS ( J,C ) = S
                        CS ( J,C ) = 1.0
                    !.......  Keep track of the maximum sources per cell for err mesg
                    ELSE
                        IF( J+1 .GT. JMAX ) JMAX = J+1

                    END IF

                    NX( C ) = J

                    !......  Otherwise, mark source as being outside domain
                    !         ELSE
                    INDOMAIN( S ) = .FALSE.
                END IF

                CYCLE                           ! To head of loop over sources

            !.......  Special case for source already assigned a grid cell
            ELSE IF( C .GT. 0 ) THEN

                J = NX( C )

                !.......  Check that the maximum number of sources per cell is ok
                IF( J .LT. MXSCEL ) THEN
                    J = J + 1
                    IS( 1,C ) = S
                    CS( 1,C ) = 1.0

                ELSE
                    IF( J+1 .GT. JMAX ) JMAX = J+1
                END IF

                NX( C ) = J
                CYCLE                       ! To head of loop over sources

            END IF

            !.......  For non-cell sources...

            !.......  Retrieve the indices to the surrogates tables
            ISIDX = 1
            F   = FINDC( CFIP, NSRGFIPS, SRGFIPS )

            !.......  Store the number and values of unfound cy/st/co codes
            !.......  Keep track of sources that are outside the domain
            IF ( F .LE. 0 ) THEN

            !.......  Re-assigning org assigned srg to default fallback srg
                IF( FSGFLAG ) ASRGID( S ) = DEFSRGID
                CYCLE               ! To next source

            END IF

            !.......  Loop through all of the cells intersecting this FIPS code.
            DO K = 1, NCELLS( F )

                C    = FIPCELL( K,F )              ! Retrieve cell number

                !.......  Set the surrogate fraction
                CALL SETFRAC( S, ISIDX, TGTSRG, K, F, NCHARS,               &
                             INDOMAIN( S ), CSRC, DEFSRGID, FSGFLAG,        &
                             ID1, ID2, FRAC, CFLAG )

                !.......  Re-assigning org assigned srg to default fallback srg
                IF( ID2 .EQ. DEFSRGID .AND. FSGFLAG ) THEN
                    ASRGID( S ) = DEFSRGID
                    CYCLE
                END IF

                !.......  Store surg IDs for reporting
                SURGID1( S ) = ID1
                SURGID2( S ) = ID2

                IF( FRAC .GT. 0 ) THEN

                    J    = NX( C )
                    !.......  Check that the maximum number of sources per cell is ok
                    !.......  Note that this J comparison to MXSCEL is not the typical
                    !         .LE. on purpose.
                    IF ( J .LT. MXSCEL .AND. FRAC .NE. 0. ) THEN
                        J = J + 1
                        IS ( J,C ) = S
                        CS ( J,C ) = FRAC

                    !.......  Keep track of the maximum sources per cell for err mesg
                    ELSE
                        IF( J+1 .GT. JMAX ) JMAX = J+1
                    END IF

                    !.......  Store the count of sources for current cell
                    NX( C )   = J

                END IF              ! if surrogate fraction > 0.

            END DO        !  end of loop on cells K for this FIP

        END DO            !  end loop on sources S, computing gridding matrix.

    END DO            !  end loop on assigned surrogate list

    !.......  For first time routine is called in all cases,
    IF( LASTIME ) THEN
        MESG = 'WARNING: Some surrogates renormalized when ' //         &
               'total of surrogates by county were ' // CRLF() //       &
                BLANK10 // 'greater than 1.'
        CALL M3MSG2( MESG )
    ENDIF

    !.......  Abort if overflow occurred
    IF ( JMAX .GT. MXSCEL ) THEN

        WRITE( MESG, 94010 )                                            &
         'INTERNAL ERROR: Gridding matrix not written.' //              &
         CRLF() // BLANK10 // 'Arrays would have overflowed.' //        &
         CRLF() // BLANK10 // 'Current maximum sources per cell (MXSCEL) =', MXSCEL, '.' // &
         CRLF() // BLANK10 // 'Actual  maximum sources per cell          =', JMAX  , '.'
        CALL M3MSG2( MESG )
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

    END IF

    !.......  Compress matrix into I/O representation from scratch
    !.......  representation and compute statistics.

    K    = 0
    CMAX = 0
    CMIN = 99999999
    DO C = 1, NGRID      ! Loop through cells

        J = NX( C )

        IF (      J .GT. CMAX ) THEN
            CMAX = J
        ELSE IF ( J .GT. 0 .AND. J .LT. CMIN ) THEN
            CMIN = J
        END IF

        DO N = 1, J      ! Loop through sources in this cell
            K = K + 1
            IF ( K .LE. NMATX ) THEN
                S       = IS( N,C )
                IX( K ) = S
                CX( K ) = CS( N,C )
            END IF
        END DO

    END DO        !  end of loop on cells C for this FIP

    NCOEF = K

    !.......  Write gridding matrix
    MESG = 'Writing out GRIDDING MATRIX file...'
    CALL M3MSG2( MESG )

    IF( .NOT. WRITE3( GNAME, 'ALL', 0, 0, NX ) ) THEN
        MESG = 'Error writing GRIDDING MATRIX file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Write output surrogates codes
    DO S = 1, NSRC
        WRITE( FDEV,93360 ) S, SURGID1( S ), SURGID2( S )
    END DO

    !.......  Report FIPS that don't have surrogate data
    !.......  Report links that are outside the grid
    !        CALL RPSRCOUT( NNOSRG, 0, FIPNOSRG, ' ' )

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Formatted file I/O formats...... 93xxx
93000 FORMAT( A )

93360 FORMAT( I8, 1X, I8, 1X, I8 )

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I9, :, 1X ) )

    !******************  INTERNAL SUBPROGRAMS  *****************************
CONTAINS

    !.......  This internal subprogram opens individual surrogate file

    SUBROUTINE OPEN_SRGFILE

        !.......  Set logical file name
        IF( .NOT. SETENVVAR( 'SRGPRO_PATH', NAMBUFT )) THEN
            MESG = 'Could not set logical name of file ' // TRIM( NAMBUFT )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        !.......  Get the number of lines in the surrogate description file desription file
        GDEV = GETEFILE( 'SRGPRO_PATH',.TRUE., .TRUE., PROGNAME )

        IF( GDEV .LT. 0 ) THEN
            MESG = 'Could not open input surrogate file' // TRIM( NAMBUFT )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        REWIND( GDEV )

        NLINES = GETFLINE( GDEV, 'Reading surrogate files' )

        IF( .NOT. SETENVVAR( 'SRGPRO_PATH', NAMBUF )) THEN
            MESG = 'Could not set logical file name of ' // TRIM( NAMBUF )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

    END SUBROUTINE OPEN_SRGFILE

END SUBROUTINE GENAGMAT


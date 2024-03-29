
SUBROUTINE SIZGMAT( CATEGORY, NSRC, VFLAG, DEFSRGID, SRGFLAG,       &
                    MXSCEL, MXCSRC, MXCCL, NMATX, NMATXU )

    !***********************************************************************
    !  subroutine body starts at line 102
    !
    !  DESCRIPTION:
    !      This subroutine determines the sizes needed for creating the for the
    !      area and mobile gridding matrices
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version ??/???? by ???
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !**************************************************************************
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

    !.......   MODULES for public variables
    !.......   This module is the source inventory arrays
    USE MODSOURC, ONLY: XLOCA, YLOCA, CIFIP, CELLID, CLINK,         &
                        XLOC1, YLOC1, XLOC2, YLOC2

    !.......   This module contains the cross-reference tables
    USE MODXREF, ONLY: ASRGID

    !.......  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: NGRID, NCOLS, NROWS, XOFF, YOFF

    !.......   This module contains the gridding surrogates tables
    USE MODSURG, ONLY: NCELLS, FIPCELL, NSRGS, SRGLIST, NSRGFIPS,   &
                       SRGFIPS, NTSRGDSC, SRGFNAM, SRGFCOD, NTLINES,&
                       MXCFIP, SRGNCOLS, SRGNROWS, SRGCSUM, SRGFRAC

    USE MODGRDLIB

    IMPLICIT NONE

    !.......   INCLUDES:

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: CATEGORY      ! source category
    INTEGER     , INTENT (IN) :: NSRC          ! local number of sources
    LOGICAL     , INTENT (IN) :: VFLAG         ! true: using variable grid
    INTEGER     , INTENT (IN) :: DEFSRGID      ! default surrogate code
    LOGICAL     , INTENT (IN) :: SRGFLAG       ! true: using default fallback surrogate
    INTEGER     , INTENT(OUT) :: MXSCEL        ! max sources per cell
    INTEGER     , INTENT(OUT) :: MXCSRC        ! max cells per source
    INTEGER     , INTENT(OUT) :: MXCCL         ! max cells per county or link
    INTEGER     , INTENT(OUT) :: NMATX         ! no. src-cell intersections
    INTEGER     , INTENT(OUT) :: NMATXU        ! no. county-cell intrsctns for all sources

    !.......   EXTERNAL FUNCTIONS and their descriptions:
    LOGICAL, EXTERNAL :: BLKORCMT
    INTEGER, EXTERNAL :: GETFLINE

    INTEGER,       PARAMETER :: MXSEG = 10               ! # of potential line segments
    CHARACTER(16), PARAMETER :: PROGNAME = 'SIZGMAT'     ! program name

    !.......   Local arrays dimensioned by subroutine arguments
    !.......   Note that the NGRID dimension could conceivably be too small if
    !              a link winds through the whole domain, but this is a case that
    !              is not worth going to extra trouble for since it is not realistic
    INTEGER         NX  ( NGRID )        ! number of srcs per cell
    INTEGER         ACEL( NGRID )        ! number of cell intersections per src
    REAL            AFAC( NGRID )        ! fraction of link in cell

    !.......   Other local variables
    INTEGER         C, F, J, JJ, K, KK, I, II, N, NT, S       ! counters and indices

    INTEGER         CCNT                 ! counters for no. non-zero-surg cells
    INTEGER      :: CELLSRC = 0          ! cell number as source char
    INTEGER         COL                  ! tmp column
    INTEGER         TCOL                 ! tmp column
    INTEGER         GDEV                 !  for surrogate coeff file
    INTEGER         ID1, ID2             ! primary and 2ndary surg codes
    INTEGER         LC, LR               ! length of COLRANGE  ROWRANGE
    INTEGER         ISIDX                ! tmp surrogate ID code index
    INTEGER         IOS                  ! i/o status
    INTEGER         IREC                 ! Record counter
    INTEGER         ISDEF                !  default surrogate ID code index
    INTEGER         NWARN                ! tmp number of WARNING msg
    INTEGER         NCEL                 ! tmp number of cells
    INTEGER         ROW                  ! tmp row
    INTEGER         TROW                 ! tmp row
    INTEGER         NTL                  ! max no. of line buffers
    INTEGER         TGTSRG               ! target surrogates code
    INTEGER         SSC                  ! surrogates code
    INTEGER      :: NLINES = 0           ! number of lines in input file

    !.......   Other arrays
    CHARACTER(20) SEGMENT( MXSEG )                 ! Segments of parsed lines
    CHARACTER(60), ALLOCATABLE :: TMPLINE( : )       ! tmp line buffer

    REAL            ALEN            ! link length
    REAL*8          XX, YY

    LOGICAL      :: EFLAG = .FALSE.     ! true: error flag
    LOGICAL      :: LFLAG = .FALSE.     ! true: location data available
    LOGICAL      :: XYSET = .FALSE.     ! true: X/Y available for src
    LOGICAL      :: CFLAG = .TRUE.      ! true: called by sizgmat, false: called by gen[a|m]gmat
    LOGICAL      :: WFLAG = .FALSE.         ! true: per iteration warning flag

    CHARACTER(FIPLEN3)  CFIP                     ! country/state/county code
    CHARACTER(FIPLEN3)  LFIP                     ! tmp country/state/county code
    CHARACTER(200)      LINE                     ! Read buffer for a line
    CHARACTER(300)      MESG                     !  message buffer
    CHARACTER(256)      NAMBUF                   !  surrogate file name buffer
    CHARACTER(256)      NAMBUFT                  !  tmp surrogate file name buffer
    CHARACTER(256)      TSRGFNAM                 !  tmp surrogate file name buffer
    CHARACTER(20)       COLRANGE                 ! buffer w/ column range
    CHARACTER(20)       ROWRANGE                 ! buffer w/ row range

    CHARACTER(LNKLEN3) :: CLNK = ' '           ! tmp link ID

    !***********************************************************************
    !   begin body of subroutine SIZGMAT

    !.......  Print status message
    MESG = 'Computing gridding matrix size...'
    CALL M3MSG2( MESG )

    !......  Set flag to indicate that XLOCA/YLOCA are available
    LFLAG = ALLOCATED( XLOCA )

    !.......  Initialize the count of sources per cell
    NX = 0       ! array

    !.......  Loop through sources
    MXCSRC  = 0
    MXCCL   = 0
    NMATXU  = 0
    CELLSRC = 0

    !.......  Allocate memory for indices to surrogates tables for each source
    ALLOCATE( NTLINES( NSRGS ), STAT=IOS )
    CALL CHECKMEM( IOS, 'NTLINES', PROGNAME )

    !.......  Create message fields for errors
    WRITE( COLRANGE, '( "( ", I1, " to ", I4, " )" )' ) 1, SRGNCOLS
    WRITE( ROWRANGE, '( "( ", I1, " to ", I4, " )" )' ) 1, SRGNROWS

    LC = LEN_TRIM( COLRANGE )
    LR = LEN_TRIM( ROWRANGE )

    !.......  Count total line buffers to define memory size
 
    DO II = 1, NSRGS            ! loop through only the surrogate code assigned by sources

        NTL = 0
        TSRGFNAM = ' '

        TGTSRG = SRGLIST( II )
        KK = II

        ISDEF  = FIND1( DEFSRGID, NSRGS, SRGLIST )

        !......  default fallback surrogate will run for the last for gap filling.
        IF( II >= ISDEF ) THEN

            !.......  default fallback surrogate will run at last after
            !         re-assigned zero fraction surrogate
            IF( II == NSRGS ) THEN
                TGTSRG = DEFSRGID
                KK = ISDEF
            ELSE
                TGTSRG = SRGLIST( II + 1 )
                KK = II + 1
            END IF

        END IF

        !.......   Count total line buffers to define memory size of TMPLINE
        DO I = 1, NTSRGDSC          ! Open all surrogate files using the same srg code

            !......  Prompt for and open I/O API output file(s)...
            IF( TGTSRG .NE. SRGFCOD( I ) ) CYCLE

            CALL GETENV( 'SRGPRO_PATH', NAMBUF )
            WRITE( NAMBUFT, '( 2A )' ) TRIM( NAMBUF ), SRGFNAM( I )

            IF( NAMBUFT .NE. TSRGFNAM  ) THEN
                CALL OPEN_SRGFILE

                IREC = 0
                !.......  Reading surrogate files
                DO JJ = 1, NLINES

                    READ ( GDEV, 93000, END=111, IOSTAT=IOS ) LINE
                    IREC = IREC + 1

                    IF ( IOS .GT. 0 ) THEN
                        WRITE( MESG, 94010) 'I/O error', IOS,       &
                            'reading gridding surrogates file at line', IREC
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF

                    IF ( BLKORCMT( LINE ) ) CYCLE

                    CALL PARSLINE( LINE, MXSEG, SEGMENT )
                    SSC    = STR2INT ( SEGMENT( 1 ) )
                    TCOL   = STR2INT ( SEGMENT( 3 ) )
                    TROW   = STR2INT ( SEGMENT( 4 ) )

                    !......  Check the value of the column number
                    IF( TCOL .LT. 0 .OR.  TCOL .GT. SRGNCOLS  .OR.&
                      ( TROW .EQ. 0 .AND. TCOL .NE. 0 ) ) THEN
                        WFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'WARNING: Column value', TCOL,&
                               'is outside range ' // COLRANGE( 1:LC ) //&
                               ' from FIPS ' // CFIP // ' and surrogate', SSC
                        CALL M3MESG( MESG )
                    END IF

                    !......  Check the value of the row number
                    IF( TROW .LT. 0 .OR.  TROW .GT. SRGNROWS  .OR.&
                      ( TCOL .EQ. 0 .AND. TROW .NE. 0 ) ) THEN
                        WFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'WARNING: Row value ', TROW,    &
                               'is outside range ' // ROWRANGE( 1:LR ) //   &
                               ' from FIPS ' // CFIP // ' and surrogate', SSC
                        CALL M3MESG( MESG )

                    !......  Special treatment for cell (0,0) (skip for now)
                    ELSE IF( TROW .EQ. 0 .AND. TCOL .EQ. 0 ) THEN
                        CYCLE

                    END IF

                    !......  Adjust column and row for subgrid
                    TCOL = TCOL - XOFF
                    TROW = TROW - YOFF

                    !.......  Skip entry after subgrid adjustment
                    IF( TCOL .LE. 0 .OR. TCOL .GT. NCOLS .OR.   &
                        TROW .LE. 0 .OR. TROW .GT. NROWS ) CYCLE

                    !...... Skip entry if rows and columns are out of range
                    IF( WFLAG ) CYCLE

                    !.......  Skip entry if SSC is not in the assigned SRGLIST by source
                    IF( SSC .EQ. TGTSRG ) NTL = NTL + 1

111             END DO

                TSRGFNAM = NAMBUFT                ! store a previous surrogate file name buffer
                CLOSE( GDEV )

            END IF                  ! skip if surrogate file has the same srg file

        END DO                   ! loop over all surrogate files in SRGDESC file

        !......  Store no of line buffers of each surrogate
        NTLINES( KK ) = NTL

        !.......  Write the status of reading surrogate files.
        WRITE( MESG,94010 ) 'Reading surrogate', TGTSRG,        &
                            'to define the size of gridding matrix'
        CALL M3MSG2( MESG )

        !.......  Warning message when there are no surrogate available.
        IF( NTL .EQ. 0 ) THEN
            WRITE( MESG,94010 ) 'WARNING: The surrogate', TGTSRG,   &
               ' does not exist within a modeling domain'
            CALL M3MESG( MESG )
        END IF

        !.......  Allocate memory for indices to surrogates tables for each source
        ALLOCATE( TMPLINE( NTL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMPLINE', PROGNAME )

        !......  If surrogates are needed, read and store the gridding surrogates,
        !        allocate memory for the surrogate assignments, and assign
        !        surrogates to each source.
        NT = 0
        TMPLINE = ' '
        TSRGFNAM = ' '

        DO I = 1, NTSRGDSC  ! Open all surrogate files using the same srg code

            !.......  Prompt for and open I/O API output file(s)...

            CALL GETENV( 'SRGPRO_PATH', NAMBUF )
            WRITE( NAMBUFT, '( 2A )' ) TRIM( NAMBUF ), SRGFNAM( I )

            IF( TGTSRG .NE. SRGFCOD( I ) ) CYCLE

            IF( NAMBUFT .NE. TSRGFNAM  ) THEN
                CALL OPEN_SRGFILE

                    !.......  Reading surrogate files
                DO JJ = 1, NLINES

                    READ ( GDEV, 93000, END=113, IOSTAT=IOS )LINE

                    IF ( BLKORCMT( LINE ) ) CYCLE

                    CALL PARSLINE( LINE, MXSEG, SEGMENT )
                    SSC    = STR2INT ( SEGMENT( 1 ) )
                    TCOL   = STR2INT ( SEGMENT( 3 ) )
                    TROW   = STR2INT ( SEGMENT( 4 ) )

                        !......  Check the value of the column number
                    IF( TCOL .LT. 0 .OR.  TCOL .GT. SRGNCOLS .OR.   &
                      ( TROW .EQ. 0 .AND. TCOL .NE. 0 ) ) THEN
                        WFLAG = .TRUE.
                    END IF

                        !......  Check the value of the row number
                    IF( TROW .LT. 0 .OR.  TROW .GT. SRGNROWS .OR.   &
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
                    IF( TCOL .LE. 0 .OR. TCOL .GT. NCOLS .OR.   &
                        TROW .LE. 0 .OR. TROW .GT. NROWS ) CYCLE

                        !......  Skip entry if rows and columns are out of range
                    IF( WFLAG ) CYCLE

                        !......  Skip entry if SSC is not in the assigned SRGLIST by source
                    IF( SSC .EQ. TGTSRG ) THEN
                        NT = NT + 1
                        TMPLINE( NT ) = LINE
                    END IF

113             END DO

                TSRGFNAM = NAMBUFT            ! store a previous surrogate file name buffer
                CLOSE( GDEV )

            END IF              ! skip if surrogate file has the same srg file

        END DO               ! loop over all surrogate files in SRGDESC file

        CALL RDSRG4GRD( NT, TMPLINE, CFLAG )         ! populating surrogates

        DEALLOCATE( TMPLINE )

        LFIP = ' '
        NWARN = 0

            !.......  Loop over sources per each assigned surrogate
        DO S = 1, NSRC

            CFIP = CIFIP( S )
            SSC  = ASRGID( S )
            IF( CATEGORY .EQ. 'AREA'   ) CELLSRC = CELLID( S )
            IF( CATEGORY .EQ. 'MOBILE' ) CLNK    = CLINK( S )

            IF( SSC .NE. TGTSRG ) CYCLE

                !.......  Determine if x/y location is available
            XYSET = .FALSE.
            IF( LFLAG ) XYSET = ( XLOCA( S ) .GT. AMISS3 )

                !.......  If cell-specific source...
            IF ( CELLSRC .GT. 0 ) THEN
                NCEL = 1
                ACEL( 1 ) = CELLID( S )
                AFAC( 1 ) = 1.

                !.......  Check if source has been converted to point src
            ELSE IF( XYSET ) THEN

                    !......  If source is in the domain....
                XX = XLOCA( S )
                YY = YLOCA( S )
                IF( INGRID( XX, YY, NCOLS, NROWS, COL, ROW  ) ) THEN

                    !......  Set as 1 cell and get the cell number
                    NCEL = 1
                    ACEL( 1 ) = ( ROW-1 ) * NCOLS + COL
                    AFAC( 1 ) = 1.

                    !.......  Otherwise, skip this source because it's outside the grid
                ELSE
                    NCEL = 0

                END IF

                !......  If area/non-link source...
            ELSE IF( CLNK .EQ. ' ' ) THEN

                    !.......  Retrieve the index to the surrogates cy/st/co list
                ISIDX = 1
                F     = FINDC( CFIP, NSRGFIPS, SRGFIPS )

                !.......  Retrieve the cell intersection info from the
                !         surrogates tables from MODSURG
                IF ( F .GT. 0 ) THEN

                    NCEL = NCELLS( F )
                    ACEL( 1:NCEL ) = FIPCELL( 1:NCEL, F )                       ! arrays

                    DO K = 1, NCEL
                        CALL SETFRAC( S, ISIDX, TGTSRG, K, F, 1,        &
                                     .FALSE.,' ', DEFSRGID, SRGFLAG,     &
                                      ID1, ID2, AFAC( K ), CFLAG )

                    !.......  Re-assigning org assigned srg to default fallback srg
                        IF( ID2 .EQ. DEFSRGID .AND. SRGFLAG ) THEN
                            ASRGID( S ) = DEFSRGID
                            CYCLE
                        END IF

                    END DO

                    !.......  Otherwise, skip this source because it's outside the grid
                ELSE
                    IF( SRGFLAG ) THEN
                        ASRGID( S ) = DEFSRGID
                        IF( II < NSRGS ) CYCLE
                    END IF

                    !......  Critical warning message about zeroing emission
                    !                            due to no surrogates for this co/st/cy
                    IF( LFIP .NE. CFIP ) THEN
                        WRITE( MESG,94010 ) 'WARNING: Causing ' //      &
                           'zeroing emissions due to missing '//        &
                           'surrogate', TGTSRG, ' for co/st/cy :: '//   &
                           CFIP
                        NWARN = NWARN + 1
                        IF( NWARN < 100 ) CALL M3MESG( MESG )
                        LFIP = CFIP
                    END IF

                    CYCLE

                END IF


                !.......  If link source, determine the number of cells for this source
            ELSE
                CALL LNK2GRD( NGRID, XLOC1( S ), YLOC1( S ),        &
                              XLOC2( S ), YLOC2( S ), NCEL, ACEL,   &
                              AFAC, ALEN, EFLAG)

                !.......  Make sure that there was enough storage
                IF( EFLAG ) THEN
                    WRITE( MESG,94010 ) 'INTERNAL ERROR: Overflow for source', S
                    CALL M3MSG2( MESG )
                    CYCLE
                END IF

            END IF
            !
            !.......  Loop through the cells for this source and increment the number
            !                   of sources per cell.
            CCNT = 0
            DO N = 1, NCEL

                IF( AFAC( N ) .GT. 0.0 ) THEN
                    C = ACEL( N )
                    NX( C ) = NX( C ) + 1
                    CCNT = CCNT + 1
                END IF

            END DO                ! End loop on cells for this source

            !.......  Update the maximum number of cells per source
            IF( CCNT .GT. MXCSRC ) MXCSRC = CCNT

            !.......  Update the maximum number of cells per county or link
            IF( NCEL .GT. MXCCL ) MXCCL = NCEL

        END DO                ! End loop on sources

    END DO                ! End loop on assigned surrogates

        !.......  Abort if error
    IF( EFLAG ) THEN
        MESG = 'Problem determining memory for gridding matrix.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Determine maximum number of sources per cell
    !.......  And determine the total number of source-cell intersections
    MXSCEL = NX( 1 )
    NMATX  = NX( 1 )
    DO C = 2, NGRID

        J = NX( C )
        IF( J .GT. MXSCEL ) THEN
            MXSCEL = J
        END IF

        NMATX = NMATX + J

    END DO                ! End loop on cells

    !..........
    !.......  Estimating ungridding matrices sizes ......
    !..........

        !...... Loop over sources per each assigned surrogate
    DO S = 1, NSRC

        CFIP = CIFIP( S )
        SSC  = ASRGID( S )
        IF( CATEGORY .EQ. 'AREA' ) CELLSRC = CELLID( S )
        IF( CATEGORY .EQ. 'MOBILE' ) CLNK = CLINK( S )

        !.......  Determine if x/y location is available
        XYSET = .FALSE.
        IF( LFLAG ) XYSET = ( XLOCA( S ) .GT. AMISS3 )

        !.......  If cell-specific source...
        IF ( CELLSRC .GT. 0 ) THEN
            NCEL = 1
            ACEL( 1 ) = CELLID( S )
            AFAC( 1 ) = 1.

        !.......  Check if source has been converted to point src
        ELSE IF( XYSET ) THEN

            !......  If source is in the domain....
            XX = XLOCA( S )
            YY = YLOCA( S )
            IF( INGRID( XX, YY, NCOLS, NROWS, COL, ROW  ) ) THEN

                !......  Set as 1 cell and get the cell number
                NCEL = 1
                ACEL( 1 ) = ( ROW-1 ) * NCOLS + COL
                AFAC( 1 ) = 1.

            !.......  Otherwise, skip this source because it's outside the grid
            ELSE
                NCEL = 0

            END IF

        !......  If area/non-link source...
        ELSE IF( CLNK .EQ. ' ' ) THEN

            !.......  Retrieve the index to the surrogates cy/st/co list
            ISIDX = 1
            F     = FINDC( CFIP, NSRGFIPS, SRGFIPS )

            !.......  Retrieve the cell intersection info from the
            !         surrogates tables from MODSURG
            IF ( F .GT. 0 ) THEN

                NCEL = NCELLS( F )

            !.......  Otherwise, skip this source because it's outside the grid
            ELSE
                NCEL = 1

            END IF

        !.......  If link source, determine the number of cells for this source
        ELSE
            CALL LNK2GRD( NGRID, XLOC1( S ), YLOC1( S ),        &
                         XLOC2( S ), YLOC2( S ), NCEL, ACEL,    &
                         AFAC, ALEN, EFLAG)

            !.......  Make sure that there was enough storage
            IF( EFLAG ) THEN
                WRITE( MESG,94010 ) 'INTERNAL ERROR: Overflow for source', S
                CALL M3MSG2( MESG )
                CYCLE
            END IF

        END IF
        !
        !.......  Count all county/cell intersections for all sources.
        !         This is needed for ungridding matrix.
        DO N = 1, NCEL

            NMATXU = NMATXU + 1

        END DO            ! End loop on cells for this source

    END DO                ! End loop on sources

    RETURN

        !******************  FORMAT  STATEMENTS   ******************************

        !.......   Internal buffering formats...... 94xxx

93000 FORMAT( A )

94010 FORMAT( 10( A, :, I8, :, 1X ) )


CONTAINS        !******************  INTERNAL SUBPROGRAMS  *****************************


    !.......  This internal subprogram opens individual surrogate file

    SUBROUTINE OPEN_SRGFILE

        !.......  Set logical file name
        IF( .NOT. SETENVVAR( 'SRGPRO_PATH', NAMBUFT )) THEN
            MESG = 'Could not set logical  name of file ' // TRIM( NAMBUFT )
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
            MESG = 'Could not set logical  name of file ' // TRIM( NAMBUF )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

    END SUBROUTINE OPEN_SRGFILE

END SUBROUTINE SIZGMAT

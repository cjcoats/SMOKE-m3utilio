
SUBROUTINE RDSRG( VFLAG, FDEV, SRGFMT, SRGNROWS, SRGNCOLS )

    !***********************************************************************
    !  subroutine body starts at line 117
    !
    !  DESCRIPTION:
    !      This subroutine allocates memory for the spatial surrogate
    !      arrays, reads the the spatial surrogates file, and then fills
    !      the spatial surrogate arrays.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version ??/???? by ???
    !       Version 11/2023 by CJC:  USE M3UTILIO and related changes
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

    !...........   Modules for public variables
    !...........   This module contains the gridding surrogates tables
    USE MODSURG, ONLY: NSRGREC, IDXSRGA, SFIPSA, SCELLA, SSRGIDA,&
                       IDXSRGB, NSRGFIPS, NSRGS, SRGFIPS, SRGLIST,&
                       FIPCELL, SRGFRAC, SRGCSUM, NCELLS, SFRACA

    !.........  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: XOFF, YOFF, NCOLS, NROWS

    IMPLICIT NONE

    !...........   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constat parameters

    !...........   Subroutine arguments
    LOGICAL      , INTENT  (IN) :: VFLAG          ! true: using variable grid
    INTEGER      , INTENT  (IN) :: FDEV           ! File unit number
    CHARACTER(*) , INTENT  (IN) :: SRGFMT         ! Format of surrogates file
    INTEGER      , INTENT  (IN) :: SRGNROWS       ! number of grid rows from hdr
    INTEGER      , INTENT  (IN) :: SRGNCOLS       ! number of grid cols from hdr

    !...........   EXTERNAL FUNCTIONS and their descriptions:

    LOGICAL, EXTERNAL :: BLKORCMT

    !...........   Local parameters

    INTEGER, PARAMETER :: MXSEG = 10              ! # of potential line segments

    !...........   Other arrays

    CHARACTER(20) SEGMENT( MXSEG )                 ! Segments of parsed lines

    CHARACTER(200), ALLOCATABLE :: SORTBUF( : )     ! concatenated info for sorting

    !...........   Local variables
    INTEGER         I, J, K, L                ! indices and counters

    INTEGER         C                         ! tmp cell number
    INTEGER         CELCNT                    ! cell counter
    INTEGER         CNTCHK                    ! check for message overflow count
    INTEGER         COL                       ! Temp grid column number (x-axis)
    INTEGER         IOS                       ! i/o status
    INTEGER         IREC                      ! Record counter
    INTEGER         LC, LR                    ! length of COLRANGE  ROWRANGE
    INTEGER         LCEL                      ! cell ID from previous iteration
    INTEGER         LSSC                      ! srg ID from previous iteration
    INTEGER         NSRGALL                   ! No. entries in surrgoates file
    INTEGER         MXCFIP                    ! Max cells per county code
    INTEGER         ROW                       ! Temp grid row number (y-axis)
    INTEGER         SSC                       ! Temp spatial surrogate code

    REAL            RATIO                     ! Temp spatial surrogate Ratio

    LOGICAL      :: GFLAG = .FALSE.           ! true: per iteration surg > 1
    LOGICAL      :: HFLAG = .FALSE.           ! true: header found
    LOGICAL      :: OFLAG = .FALSE.           ! true: overall warning flag
    LOGICAL      :: RFLAG = .FALSE.           ! true: renormalized surrogates
    LOGICAL         WFLAG                     ! true: per iteration warning flag

    CHARACTER(FIPLEN3) CFIP                   ! tmp country/state/county code
    CHARACTER(FIPLEN3) LFIP                   ! county code from prev iteration
    CHARACTER(20)   COLRANGE                  ! buffer w/ column range
    CHARACTER(20)   ROWRANGE                  ! buffer w/ row range
    CHARACTER(200)  LINE                      ! Read buffer for a line
    CHARACTER(600)  MESG                      ! Message buffer

    CHARACTER(16) :: PROGNAME = 'RDSRG'        !  program name

    !***********************************************************************
    !   Begin body of subroutine RDSRG

    IREC    = 0
    J       = 0
    NSRGREC = 0

    !......... Determine the number surrogate file entries

    REWIND( FDEV )

    SELECT CASE( SRGFMT )

      CASE( 'MODELS3' )

        HFLAG = .FALSE.
        DO

            READ ( FDEV, 93000, END=12, IOSTAT=IOS ) LINE

            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010)&
                  'I/O error', IOS, 'reading gridding surrogates '//&
                  'file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. HFLAG ) THEN

                CALL UPCASE( LINE )

                IF( VFLAG ) THEN              ! check current line for header
                    I = INDEX( LINE, '#VARIABLE_GRID' )
                ELSE
                    I = INDEX( LINE, '#GRID' )
                END IF

                IF( I .GT. 0 ) THEN           ! skip if current line is header
                    HFLAG = .TRUE.
                    CYCLE
                END IF

            ELSE IF ( BLKORCMT( LINE ) ) THEN     ! skip if current line is blank OR comment
                CYCLE

            ELSE
                NSRGREC = NSRGREC + 1
            END IF

        END DO

12      CONTINUE        ! end of read on input file

    END SELECT

    !......... Allocate memory for surrogate arrays
    IF( ALLOCATED( IDXSRGA ) ) DEALLOCATE( IDXSRGA )
    ALLOCATE( IDXSRGA( NSRGREC ),                  &
              IDXSRGB( NSRGREC ),                  &
              SORTBUF( NSRGREC ),                  &
               SCELLA( NSRGREC ),                  &
               SFIPSA( NSRGREC ),                  &
              SSRGIDA( NSRGREC ),                  &
               SFRACA( NSRGREC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'IDXSRGA...SFRACA', PROGNAME )

    !.........  Create message fields for errors
    WRITE( COLRANGE, '( "( ", I1, " to ", I4, " )" )' ) 1, SRGNCOLS
    WRITE( ROWRANGE, '( "( ", I1, " to ", I4, " )" )' ) 1, SRGNROWS

    LC = LEN_TRIM( COLRANGE )
    LR = LEN_TRIM( ROWRANGE )

    !.........  Fill surrogate arrays

    REWIND( FDEV )
    IREC    = 0

    SELECT CASE( SRGFMT )

      CASE( 'MODELS3' )

        HFLAG = .FALSE.
        DO

            READ ( FDEV, 93000, END=24, IOSTAT=IOS ) LINE

            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010)                                     &
                  'I/O error', IOS, 'reading gridding surrogates '//    &
                  'file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. HFLAG ) THEN

                CALL UPCASE( LINE )

                IF( VFLAG ) THEN              ! check current line for header
                    I = INDEX( LINE, '#VARIABLE_GRID' )
                ELSE
                    I = INDEX( LINE, '#GRID' )
                END IF

                IF( I .GT. 0 ) THEN       ! skip if current line is header
                    HFLAG = .TRUE.
                    CYCLE
                END IF

            ELSE IF ( BLKORCMT( LINE ) ) THEN     ! skip if current line is blank OR comment
                CYCLE

            END IF

            !.................  Parse the line of data into segments based on the rules
            !                   for "list-formatted" in fortran, but not requiring
            !                   quotes around the text strings

            CALL PARSLINE( LINE, MXSEG, SEGMENT )

            SSC    = STR2INT ( SEGMENT( 1 ) )
            CFIP   = SEGMENT( 2 )
            CALL PADZERO( CFIP )
            COL    = STR2INT ( SEGMENT( 3 ) )
            ROW    = STR2INT ( SEGMENT( 4 ) )
            RATIO  = STR2REAL( SEGMENT( 5 ) )

            WFLAG = .FALSE.

            !.................  Check the value of the column number
            IF( COL .LT. 0 .OR.  COL .GT. SRGNCOLS  .OR.            &
              ( ROW .EQ. 0 .AND. COL .NE. 0     )    ) THEN
                WFLAG = .TRUE.
                OFLAG = .TRUE.
                WRITE( MESG,94010 ) 'WARNING: Column value ', COL,  &
                       'is outside range ' // COLRANGE( 1:LC ) //   &
                       ' at line', IREC
                CALL M3MESG( MESG )
            END IF

            !.................  Check the value of the row number
            IF( ROW .LT. 0 .OR.  ROW .GT. SRGNROWS  .OR.            &
              ( COL .EQ. 0 .AND. ROW .NE. 0     )    ) THEN
                WFLAG = .TRUE.
                OFLAG = .TRUE.
                WRITE( MESG,94010 ) 'WARNING: Row value ', ROW,     &
                       'is outside range ' // ROWRANGE( 1:LR ) //   &
                       ' at line', IREC
                CALL M3MESG( MESG )

            !.................  Special treatment for cell (0,0) (skip for now)
            ELSE IF( ROW .EQ. 0 .AND. COL .EQ. 0 ) THEN
                CYCLE

            END IF

            !.................  Adjust column and row for subgrid
            COL = COL - XOFF
            ROW = ROW - YOFF

            !.................  Skip entry after subgrid adjustment
            IF( COL .LE. 0 .OR. COL .GT. NCOLS .OR.                 &
                ROW .LE. 0 .OR. ROW .GT. NROWS  ) CYCLE

            !.................  Check the value of the ratio value
            IF( RATIO .GT. 1. ) THEN
                WRITE( MESG,94020 )                                 &
                       'WARNING: resetting surrogates ratio at ' // &
                       'line', IREC, 'from', RATIO, 'to 1.'
                CALL M3MESG( MESG )
                RATIO = 1.
            END IF

    !.................  Skip entry if rows and columns are out of range
            IF( WFLAG ) CYCLE

            J = J + 1

            IDXSRGA( J ) = J
            IDXSRGB( J ) = J
            SCELLA ( J ) = (ROW-1)*NCOLS + COL
            SFIPSA ( J ) = CFIP
            SSRGIDA( J ) = SSC
            SFRACA ( J ) = RATIO
            WRITE( SORTBUF( J ), '(A,I8,I8)' ) SFIPSA( J ), SCELLA( J ), SSRGIDA( J )

        END DO

24      CONTINUE        ! end of read on input file

        NSRGALL = J

        !.............  Write out final warning for row/col out of range
        IF( OFLAG ) THEN
            MESG = 'WARNING: Lines skipped in surrogate file ' //   &
                   'because rows or columns were out of range.'
            CALL M3MSG2( MESG )
        END IF

      CASE DEFAULT
        MESG = 'INTERNAL ERROR: Surrogate format "' // SRGFMT //    &
               '" not understood'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END SELECT

    !.........  Now create the derived surrogates tables from the original data...

    !.........  Sort surrogates by county code  cell  surrogate code
    CALL SORTIC( NSRGALL, IDXSRGA, SORTBUF )

    !.........  Sort surrogates by surrogate code
    CALL SORTI1( NSRGALL, IDXSRGB, SSRGIDA )

    !.........  Count county codes in surrogates file and maximum number of cells
    !           per cy/st/co code.
    LFIP     = ' '
    LCEL     = -1
    MXCFIP   = 0
    NSRGFIPS = 0
    CELCNT   = 0
    DO I = 1, NSRGALL

        J    = IDXSRGA( I )
        CFIP = SFIPSA ( J )
        C    = SCELLA ( J )

        IF( CFIP .NE. LFIP ) THEN

            IF( CELCNT .GT. MXCFIP ) MXCFIP = CELCNT

            NSRGFIPS = NSRGFIPS + 1      ! incrmt cntr for county in srg file
            CELCNT   = 0                 ! init cell counter per county
            LFIP     = CFIP

        END IF

        IF( C   .NE. LCEL .OR. CELCNT .EQ. 0 ) THEN
            CELCNT = CELCNT + 1          ! increment cell counter per county
            LCEL = C
        END IF

    END DO

    IF( CELCNT .GT. MXCFIP ) MXCFIP = CELCNT

    !.........  Count surrogate codes in surrogates file
    LSSC  = -1
    NSRGS = 0
    DO I = 1, NSRGALL

        J   = IDXSRGB( I )
        SSC = SSRGIDA( J )

        IF( SSC .NE. LSSC ) THEN
            NSRGS = NSRGS + 1
            LSSC = SSC
        END IF

    END DO

    !.........  Allocate memory for derived surrogates tables
    IF( ALLOCATED( SRGLIST ) ) DEALLOCATE( SRGLIST )

    ALLOCATE(  NCELLS( NSRGFIPS ),                  &
              SRGFIPS( NSRGFIPS ),                  &
              SRGLIST( NSRGS ),                     &
              FIPCELL( MXCFIP, NSRGFIPS ),          &
              SRGFRAC( NSRGS, MXCFIP, NSRGFIPS ),   &
              SRGCSUM( NSRGS, NSRGFIPS ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SRGCSUM', PROGNAME )

    !.........  Initialize arrays
    NCELLS  = 0      ! array
    SRGFIPS = ' '    ! array
    SRGLIST = 0      ! array
    FIPCELL = 0      ! array
    SRGFRAC = 0.     ! array
    SRGCSUM = 0.     ! array

    !.........  Store derived surrogates tables...

    !.........  Store the surrogate ID list
    LSSC  = -1
    NSRGS = 0
    DO I = 1, NSRGALL

        J    = IDXSRGB( I )
        SSC = SSRGIDA( J )

        IF( SSC .NE. LSSC ) THEN
            NSRGS = NSRGS + 1
            SRGLIST( NSRGS ) = SSC
            LSSC = SSC
        END IF

    END DO

    !.........  Initialize arrays that might not be totally populated
    FIPCELL = 0       ! array
    SRGFRAC = 0       ! array

    !.........  Store the surrogate fractions, FIPS codes, and cell numbers...
    LFIP     = ' '
    LCEL     = -1
    NSRGFIPS = 0
    DO I = 1, NSRGALL

        J     = IDXSRGA( I )
        CFIP  = SFIPSA ( J )
        C     = SCELLA ( J )
        SSC   = SSRGIDA( J )
        RATIO = SFRACA ( J )

        IF( CFIP .NE. LFIP ) THEN

            NSRGFIPS = NSRGFIPS + 1          ! incrmt cntr for county in srg file
            SRGFIPS( NSRGFIPS ) = CFIP
            CELCNT   = 0                     ! init cell counter per county
            LFIP     = CFIP

        END IF

        IF( C   .NE. LCEL .OR. CELCNT .EQ. 0 ) THEN
            CELCNT = CELCNT + 1              ! increment cell counter per county
            LCEL = C
        END IF

        !.............  Store cell count at every iteration, and it will get
        !               overwritten until the maximum value for each county is stored
        NCELLS ( NSRGFIPS ) = CELCNT

        !.............  Store cell number  at every iteration, but it will be the same
        !               each time
        FIPCELL( CELCNT, NSRGFIPS ) = C

        !.............  Find surrogate code in sorted list and use position to store
        !               the surrogates fraction
        K = FIND1( SSC, NSRGS, SRGLIST )

        SRGFRAC( K, CELCNT, NSRGFIPS ) = RATIO
        SRGCSUM( K, NSRGFIPS ) = SRGCSUM( K, NSRGFIPS ) + RATIO

    END DO

    !.........  Now check to make sure the surrogate totals for each county are
    !           less than or equal to 1.
    DO I = 1, NSRGFIPS

        GFLAG = .FALSE.
        CNTCHK = 0
        DO K = 1, NSRGS

            !.................  Check if county total surrogates greater than 1
            IF( SRGCSUM( K,I ) .GT. 1.001 ) THEN

            !.....................  If first problem on this line
                IF( .NOT. GFLAG ) THEN
                    WRITE( MESG,94030 ) 'WARNING: County ' //           &
                      'surrogate total greater than 1. for '//          &
                      'county' // CRLF() // BLANK10 // SRGFIPS( I ) //  &
                      ', SSC(', SRGLIST( K ), '):', SRGCSUM( K,I )
                    GFLAG = .TRUE.
                    CNTCHK = CNTCHK + 1

                !.....................  If multiple problems on this line, but fewer than
                !                       the MESG length will permit, add to message
                ELSE IF( CNTCHK .LE. 27 ) THEN
                    L = LEN_TRIM( MESG )
                    WRITE( MESG,94031 ) MESG( 1:L )// ', SSC(',         &
                           SRGLIST( K ), '):', SRGCSUM( K,I )
                    CNTCHK = CNTCHK + 1

                END IF

            END IF

            !.................  Renormalize all surrogates with totals > 1
            IF( SRGCSUM( K,I ) .GT. 1. ) THEN
                RFLAG = .TRUE.              ! Set to give global warning
                DO C = 1, NCELLS( I )
                    SRGFRAC( K,C,I ) = SRGFRAC( K,C,I ) / SRGCSUM( K,I )
                END DO
            END IF

        END DO

            !.............  Give a warning message for significant counties
        IF( GFLAG ) CALL M3MESG( MESG )

    END DO

    !.........  Reset number of surrogate records stored in the module with
    !           the correct number after reading file and removing records that
    !           are outside the subgrid (if any)
    NSRGREC = NSRGALL

    IF( RFLAG ) THEN
        MESG = 'WARNING: Some surrogates renormalized when ' //&
               'total of surrogates by county were ' // CRLF() //&
               BLANK10 // 'greater than 1.'
        CALL M3MSG2( MESG )
    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )


    !...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

94020 FORMAT( A, 1X, I8, 1X, A, 1X, F10.6, 1X, A )

94030 FORMAT( A, I3.2, A, F8.4 )

94031 FORMAT( A, 1X, I3.2, A, F8.4 )

END SUBROUTINE RDSRG

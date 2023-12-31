
SUBROUTINE RDSRG4GRD( NT, TMPLINE, CFLAG )

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
    !       Created by B.H. Baek 02/06
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

    !...........   Modules for public variables
    !...........   This module contains the gridding surrogates tables
    USE MODSURG, ONLY: IDXSRGA, IDXSRGB, SCELLA, SFIPSA, SSRGIDA,       &
                       SFRACA, NSRGREC, NSRGFIPS, SRGFIPS, FIPCELL,     &
                       NCELLS, SRGFRAC, SRGCSUM, SRGFMT, SRGNCOLS,      &
                       SRGNROWS, MXCFIP

    !.........  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: XOFF, YOFF, NCOLS, NROWS

    IMPLICIT NONE

    !...........   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constat parameters

    !...........   Subroutine arguments
    INTEGER      , INTENT  (IN) :: NT                ! no of surrogate entry
    CHARACTER(*) , INTENT  (IN) :: TMPLINE( NT )     ! tmp line buffer
    LOGICAL      , INTENT  (IN) :: CFLAG             ! true: called by sizgmat, false: called by gen[a|m]gmat

    !...........   EXTERNAL FUNCTIONS and their descriptions:

    LOGICAL, EXTERNAL :: BLKORCMT

    !...........   Local parameters

    INTEGER,        PARAMETER :: MXSEG = 10              ! # of potential line segments
    CHARACTER(16) , PARAMETER :: PROGNAME = 'RDSRG4GRD'        !  program name

    !...........   Other arrays

    CHARACTER(20) SEGMENT( MXSEG )                 ! Segments of parsed lines

    !...........   Local variables
    INTEGER         I, F, J, K, L, N             ! indices and counters

    INTEGER         C                         ! tmp cell number
    INTEGER         CELCNT                    ! cell counter
    INTEGER         CNTCHK                    ! check for message overflow count
    INTEGER         COL                       ! Temp grid column number (x-axis)
    INTEGER         IOS                       ! i/o status
    INTEGER         LCEL                      ! cell ID from previous iteration
    INTEGER         NSRGALL                   ! No. entries in surrgoates file
    INTEGER         ROW                       ! Tmp grid row number (y-axis)
    INTEGER         SSC                       ! Tmp spatial surrogate code
    INTEGER         LCC                       ! previous tmp cell ID

    REAL            RATIO                     ! Tmp spatial surrogate Ratio
    REAL            LRATIO                    ! previous tmp spatial surrogate Ratio

    LOGICAL      :: EFLAG = .FALSE.           ! true: errir found
    LOGICAL      :: GFLAG = .FALSE.           ! true: per iteration surg > 1
    LOGICAL      :: HFLAG = .FALSE.           ! true: header found
    LOGICAL      :: OFLAG = .FALSE.           ! true: overall warning flag
    LOGICAL         WFLAG                     ! true: per iteration warning flag

    CHARACTER(200), ALLOCATABLE :: SORTBUF( : )     ! concatenated info for sorting
    CHARACTER(FIPLEN3) CFIP                   ! tmp country/state/county code
    CHARACTER(FIPLEN3) LFIP                   ! county code from prev iteration
    CHARACTER(200)  LINE                      ! Read buffer for a line
    CHARACTER(300)  MESG                      ! Message buffer

    !***********************************************************************
    !   Begin body of subroutine RDSRG4GRD

    NSRGREC = NT

    !......... Allocate memory for surrogate arrays
    IF (ALLOCATED(IDXSRGA)) DEALLOCATE (IDXSRGA)
    IF (ALLOCATED(IDXSRGB)) DEALLOCATE (IDXSRGB)
    IF (ALLOCATED(SCELLA))  DEALLOCATE (SCELLA)
    IF (ALLOCATED(SFIPSA))  DEALLOCATE (SFIPSA)
    IF (ALLOCATED(SSRGIDA)) DEALLOCATE (SSRGIDA)
    IF (ALLOCATED(SFRACA))  DEALLOCATE (SFRACA)
    IF (ALLOCATED(NCELLS))  DEALLOCATE (NCELLS)
    IF (ALLOCATED(SRGFIPS)) DEALLOCATE (SRGFIPS)
    IF (ALLOCATED(SRGFRAC)) DEALLOCATE (SRGFRAC)
    IF (ALLOCATED(SRGCSUM)) DEALLOCATE (SRGCSUM)
    IF (ALLOCATED(FIPCELL)) DEALLOCATE (FIPCELL)

    ALLOCATE( IDXSRGA( NSRGREC ),   &
              IDXSRGB( NSRGREC ),   &
              SORTBUF( NSRGREC ),   &
               SCELLA( NSRGREC ),   &
               SFIPSA( NSRGREC ),   &
              SSRGIDA( NSRGREC ),   &
               SFRACA( NSRGREC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'IDXSRGA...SFRACA', PROGNAME )

    !.........  Initialize arrays
    IDXSRGA = 0      ! array
    IDXSRGB = 0      ! array
    SCELLA  = 0      ! array
    SFIPSA  = ' '    ! array
    SSRGIDA = 0.     ! array
    SFRACA  = 0.     ! array

    !.........  Fill surrogate arrays

    SELECT CASE( SRGFMT )

      CASE( 'MODELS3' )

        HFLAG = .FALSE.
        J = 0
        DO I = 1, NT

            LINE = TMPLINE( I )

            IF ( .NOT. HFLAG ) CALL UPCASE( LINE )

            !........  Parse the line of data into segments based on the rules
            !          for "list-formatted" in fortran, but not requiring
            !          quotes around the text strings

            CALL PARSLINE( LINE, MXSEG, SEGMENT )

            SSC    = STR2INT ( SEGMENT( 1 ) )
            CFIP   = SEGMENT( 2 )
            CALL PADZERO( CFIP )
            COL    = STR2INT ( SEGMENT( 3 ) )
            ROW    = STR2INT ( SEGMENT( 4 ) )
            RATIO  = STR2REAL( SEGMENT( 5 ) )

            WFLAG = .FALSE.

            !........  Adjust column and row for subgrid
            COL = COL - XOFF
            ROW = ROW - YOFF

            !........  Skip entry after subgrid adjustment
            IF( COL .LE. 0 .OR. COL .GT. NCOLS .OR.&
                ROW .LE. 0 .OR. ROW .GT. NROWS       ) CYCLE

            !........  Check the value of the ratio value
            IF( RATIO .GT. 1. .AND. CFLAG ) THEN
                WRITE( MESG,94020 )&
                   'WARNING: resetting surrogates ratio ' //&
                   ' of Co/St/Ct (FIPS): ' // CFIP // ' and surrogate :',&
                   SSC, ' from', RATIO, 'to 1.0'
                CALL M3MESG( MESG )

                RATIO = 1.
            END IF

            !........  Skip entry if rows and columns are out of range
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

        CONTINUE

        NSRGALL = J

        !.............  Write out final warning for row/col out of range
        IF( OFLAG ) THEN
            MESG = 'WARNING: Lines skipped in surrogate file ' //&
                   'because rows or columns were out of range.'
            CALL M3MSG2( MESG )
        END IF

      CASE DEFAULT
        MESG = 'INTERNAL ERROR: Surrogate format "' // SRGFMT //&
               '" not understood'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END SELECT

    !.........  Now create the derived surrogates tables from the original data...

    !.........  Sort surrogates by county code  cell  surrogate code
    CALL SORTIC( NSRGALL, IDXSRGA, SORTBUF )

    !.........  Count county codes in surrogates file and maximum number of cells
    !       per cy/st/co code.
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

    K = 1                ! single srg code per run

    !.........  Allocate memory for derived surrogates tables
    ALLOCATE(  NCELLS( NSRGFIPS ),              &
              SRGFIPS( NSRGFIPS ),              &
              FIPCELL( MXCFIP, NSRGFIPS ),      &
              SRGFRAC( K, MXCFIP, NSRGFIPS ),   &
              SRGCSUM( K, NSRGFIPS ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SRGCSUM', PROGNAME )

    !.........  Initialize arrays
    NCELLS  = 0      ! array
    SRGFIPS = ' '    ! array
    FIPCELL = 0      ! array
    SRGFRAC = 0.     ! array
    SRGCSUM = 0.     ! array

    !........  Store the surrogate fractions, FIPS codes, and cell numbers...
    LFIP     = ' '
    LCEL     = -1
    LCC      = -1
    LRATIO   = -1.
    NSRGFIPS = 0
    DO I = 1, NSRGALL

        J     = IDXSRGA( I )
        CFIP  = SFIPSA ( J )
        C     = SCELLA ( J )
        SSC   = SSRGIDA( J )
        RATIO = SFRACA ( J )

        IF( CFIP .NE. LFIP ) THEN

            NSRGFIPS = NSRGFIPS + 1      ! incrmt cntr for county in srg file
            SRGFIPS( NSRGFIPS ) = CFIP
            CELCNT   = 0                 ! init cell counter per county
            LFIP     = CFIP
            LRATIO = -1.
            LCC    = -1

        ELSE
            IF( LCC .EQ. C .AND. LRATIO .EQ. RATIO .AND. CFLAG) THEN
                WRITE( MESG,94011 ) 'ERROR: Duplicate entries' //       &
                     ' with same FIPS' // CFIP // 'surrogate ', SSC,    &
                     ' and fraction value', RATIO
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
            END IF

            IF( LCC .EQ. C .AND. LRATIO .NE. RATIO .AND. CFLAG) THEN
                WRITE( MESG,94012 ) 'WARNING: Duplicate entries' //     &
                     ' with same FIPS' // CFIP // 'and surrogate ', SSC,&
                     ' with different fraction values', RATIO,          &
                     ' and ', LRATIO
                CALL M3MSG2( MESG )
            END IF

            LCC    = C
            LRATIO = RATIO

        END IF

        IF( C   .NE. LCEL .OR. CELCNT .EQ. 0 ) THEN

            CELCNT = CELCNT + 1              ! increment cell counter per county
            LCEL = C

        END IF

        !.....  Store cell count at every iteration, and it will get
        !       overwritten until the maximum value for each county is stored
        NCELLS ( NSRGFIPS ) = CELCNT

        !.....  Store cell number  at every iteration, but it will be the same
        !       each time
        FIPCELL( CELCNT, NSRGFIPS ) = C

        SRGFRAC( K, CELCNT, NSRGFIPS ) = RATIO
        SRGCSUM( K, NSRGFIPS ) = SRGCSUM( K, NSRGFIPS ) + RATIO

    END DO

    !.........  Now check to make sure the surrogate totals for each county are
    !           less than or equal to 1.
    IF( .NOT. CFLAG ) THEN

        DO I = 1, NSRGFIPS

            GFLAG = .FALSE.
            CNTCHK = 0

            !........  Check if county total surrogates greater than 1
            IF( SRGCSUM( K,I ) .GT. 1.0 ) THEN
                WRITE( MESG,94030 ) 'WARNING: Re-normalizing ' //       &
                   'county total surrogate fractions of Co/St/Ct: ' //  &
                   SRGFIPS( I ) // ' :: Surrogate:', SSC ,              &
                    ' greater than 1 :: ', SRGCSUM( K,I )
                CALL M3MESG( MESG )

                !........  Renormalize all surrogates with totals > 1
                DO C = 1, NCELLS( I )
                    SRGFRAC( K,C,I ) = SRGFRAC( K,C,I ) / SRGCSUM( K,I )
                END DO
            END IF


        END DO

    END IF

    !.........  Deallocate local variables
    DEALLOCATE( IDXSRGA, IDXSRGB, SORTBUF, SCELLA, SFIPSA, SSRGIDA, SFRACA )

    IF( EFLAG ) THEN
        MESG = ' Problem with defining the size of gridding matrix.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

94011 FORMAT( A, I5, 1X, A, F8.4 )

94012 FORMAT( A, I5, 1X, A, F8.4, A, F8.4 )

94020 FORMAT( A, 1X, I5.5, 1X, A, 1X, F11.8, 1X, A )

94030 FORMAT( A, 1X, I5.5, A, 1X, F11.8 )

END SUBROUTINE RDSRG4GRD


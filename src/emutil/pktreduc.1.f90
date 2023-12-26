
PROGRAM PKTREDUC

!***********************************************************************
!  program body starts at line
!
!  DESCRIPTION:
!       This program reduces the size of control and projection packets, but
!       keeps the same information. This is done by itentification and
!       implementation of state defaults where possible.
!
!       Currently, the program can be used for projection packets only.  It
!       has been updated to handle state defaults already being in the
!       file. When there are state defaults already (county=000), pktreduc
!       does not reduce that state.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!       Models-3 I/O
!
!  REVISION  HISTORY:
!       Created 12/00 by M Houyoux. Updates 12/2001 for SMOKE 1.4
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
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
!*************************************************************************
    USE M3UTILIO

    IMPLICIT NONE

!...........   INCLUDES:

    INCLUDE 'EMCNST3.EXT'     ! emissions constant parameters

!...........   PARAMETERS and their descriptions:

    CHARACTER(16), PARAMETER :: PROGNAME = 'PKTREDUC'   !  program name
    CHARACTER(50), PARAMETER :: SCCSW = '%W%'

!.........  EXTERNAL FUNCTIONS and their descriptions:

    INTEGER, EXTERNAL :: GETNLIST
    LOGICAL, EXTERNAL :: CHKINT

!.........  Allocatable arrays
    INTEGER,      ALLOCATABLE :: INDXA  ( : )
    INTEGER,      ALLOCATABLE :: FIPSA  ( : )
    INTEGER,      ALLOCATABLE :: RCOUNT ( : )
    INTEGER,      ALLOCATABLE :: RINDX  ( : )
    INTEGER,      ALLOCATABLE :: RMAX   ( : )

    LOGICAL,      ALLOCATABLE :: DEFAULT( : )
    LOGICAL,      ALLOCATABLE :: FIXED  ( : )

    CHARACTER(66), ALLOCATABLE :: CFSCA  ( : )
    CHARACTER(66), ALLOCATABLE :: CSSCA  ( : )
    CHARACTER(66), ALLOCATABLE :: CRECA  ( : )

!.........  Fixed-dimension arrays
!.........  NOTE: The per-packet arrays should eventually be allocatable
    INTEGER        NRECS ( 5 )
    INTEGER        PSTLIN( 5 )

    CHARACTER(20)  PTYPE ( 5 )
    CHARACTER(50)  SEGMENT( 4 )
    CHARACTER(256) PENBUF( 5 )
    CHARACTER(256) PSTBUF( 5 )

!.........  Logical file names and unit numbers
    INTEGER   IDEV         ! input control packets file
    INTEGER   LDEV         ! log file unit number
    INTEGER   ODEV         ! output control packets file

!.........  Local variables

    INTEGER   I, J, K, K1, K2, K3, K4, L, N         ! indices and counters

    INTEGER   CNY          ! tmp county code
    INTEGER   FIP          ! tmp FIPS code
    INTEGER   STA          ! tmp state ID
    INTEGER   IDXS1        ! saved index (for start of FIP/SCC/factor)
    INTEGER   IDXS2        ! saved index (for start of FIP/SCC)
    INTEGER   IOS          ! i/o status
    INTEGER   IREC         ! record counter
    INTEGER   NPKT         ! number of packets in file
    INTEGER   NSEG         ! number of segments in auto-detect format
    INTEGER   SCCCOL       ! column for the SCC data

    REAL      RBUF         ! real values tmp buffer

    LOGICAL     :: EFLAG    = .FALSE. ! true: error found
    LOGICAL     :: COUNTFLG = .FALSE. ! true: count packet records

    CHARACTER(6)        COL1DESC ! description of column 1
    CHARACTER(12)       CBUF    ! tmp factor
    CHARACTER(10)       FIPFMT  ! format for int->str FIPS code conversion
    CHARACTER(66)       LCFSC   ! previous CFSCA
    CHARACTER(66)       LCREC   ! previous CRECA
    CHARACTER(66)       LCSSC   ! previous CSSCA
    CHARACTER(256)      LINE    ! line buffer
    CHARACTER(300)      MESG    ! temporary message array
    CHARACTER(FIPLEN3)  CFIP    ! char FIPS code
    CHARACTER(FIPLEN3)  CSTA    ! char state code
    CHARACTER(SCCLEN3)  TSCC    ! tmp source category code

!***********************************************************************
!   begin body of program PKTREDUC

    LDEV = INIT3()

!.........  Setup formats
    WRITE( FIPFMT, '("(I",I2.2,".",I2.2,")")' ) FIPLEN3, FIPLEN3

!.........  Write out copyright, version, web address, header info, and prompt
!           to continue running the program.
    CALL INITEM( LDEV, SCCSW, PROGNAME )

!.........  Set source category based on environment variable setting
    CALL GETCTGRY

!.........  Prompt for name of control packets file file
    MESG = 'Enter logical name of the INPUT PACKETS file'
    IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'GCNTL', PROGNAME )

!.........  Prompt for names of output file
    MESG = 'Enter logical name for the OUTPUT PACKETS file'
    ODEV= PROMPTFFILE( MESG,.FALSE.,.TRUE.,'GCNTL_OUT',PROGNAME )

    CALL M3MSG2( 'Scanning input file...' )

!.........  Determine location, type, and size of packets in file
! NOTE: This section could share code with ALOCPKTS (and its subprogram
!    N: CHECK_PACKETS) for reading and counting packets.  For now, this
!    N: section will just be hardcoded for the projection packet.
    NPKT = 0
    IREC = 0
    COUNTFLG = .FALSE.
    DO

        READ( IDEV,93000,END=555,IOSTAT=IOS ) LINE
        IREC = IREC + 1

        LINE = ADJUSTL( LINE )
        IF( LINE( 1:11 ) .EQ. '/PROJECTION' ) THEN
            NPKT = NPKT + 1
            PSTLIN( NPKT ) = IREC
            PSTBUF( NPKT ) = LINE
            PTYPE ( NPKT ) = 'PROJECTION'
            COUNTFLG = .TRUE.

        ELSE IF( LINE( 1:5 ) .EQ. '/END/' ) THEN
            PENBUF( NPKT ) = LINE
            COUNTFLG = .FALSE.

        ELSE IF( COUNTFLG ) THEN
            NRECS( NPKT ) = NRECS( NPKT ) + 1

        END IF

    END DO

555 CONTINUE  ! End of read loop

    CALL M3MSG2( 'Reading input file...' )

!.........  Loop through packets in file
    DO N = 1, NPKT

!.............  Allocated memory for unsorted arrays
        ALLOCATE( INDXA( NRECS( N ) ),&
        &          FIPSA( NRECS( N ) ),&
        &          CFSCA( NRECS( N ) ),&
        &          CSSCA( NRECS( N ) ),&
        &          CRECA( NRECS( N ) ),&
        &         RCOUNT( NRECS( N ) ),&
        &          RINDX( NRECS( N ) ),&
        &           RMAX( NRECS( N ) ),&
        &        DEFAULT( NRECS( N ) ),&
        &          FIXED( NRECS( N ) ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INDXA...FIXED', PROGNAME )

        INDXA   = 0         ! Array
        FIPSA   = 0         ! Array
        CFSCA   = ' '       ! Array
        CSSCA   = ' '       ! Array
        CRECA   = ' '       ! Array
        RCOUNT  = 0         ! Array
        RMAX    = 0         ! Array
        DEFAULT = .FALSE.   ! Array
        FIXED   = .FALSE.   ! Array

!.............  Rewind file
        REWIND( IDEV )

!.............  Based on type of packet...
        SELECT CASE( PTYPE( N ) )
          CASE( 'PROJECTION' )

!.................  Skip lines to get to correct spot in file
            CALL SKIPL( IDEV, PSTLIN( N ) )

!.................  Loop through entries in packet
            IREC = PSTLIN( N )
            K = 0
            DO I = 1, NRECS( N )

!.....................  Read unsorted packet records
                READ( IDEV, 93000, END=1011, IOSTAT=IOS ) LINE

                IF( IOS .NE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )&
                    &   'I/O error', IOS,&
                    &   'reading speciation x-ref file at line', IREC
                    CALL M3MESG( MESG )
                    CYCLE

                END IF

!.....................  For first line, get format by counting the number of
!                       columns
                IF ( I .EQ. 1 ) THEN
                    L = LEN_TRIM( LINE )
                    NSEG = GETNLIST( L, LINE )

                    IF( NSEG .LT. 3 ) THEN
                        MESG = 'First line of packet has < 3 '//&
                        &       'columns'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                    ELSE IF ( NSEG .EQ. 3 ) THEN
                        MESG = 'NOTE: Assuming region codes are' //&
                        &       ' single values (no spaces)'
                        CALL M3MSG2( MESG )
                        COL1DESC = 'Region'
                        SCCCOL = 2

                    ELSE IF ( NSEG .EQ. 4 ) THEN
                        MESG = 'NOTE: Assuming region codes are ' //&
                        &       'separated into country/state and '//&
                        &       'county sections (one space)'
                        CALL M3MSG2( MESG )
                        COL1DESC = 'Co/Sta'
                        SCCCOL = 3

                    ELSE IF ( NSEG .GT. 4 ) THEN
                        MESG = 'First line of packet has > 4 '//&
                        &       'columns'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                    END IF

                END IF

!.....................  Decompose and reassemble lines of file for sorting needs
                CALL PARSLINE( LINE, NSEG, SEGMENT )

                IF( .NOT. CHKINT( SEGMENT( 1 ) ) ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: ' // COL1DESC //&
                    &       ' code is non-integer at entry', I
                    CALL M3MESG( MESG )
                END IF

!.....................  For country/state and county in separate columns
                IF ( NSEG .GT. 3 ) THEN

                    IF( .NOT. CHKINT( SEGMENT( 2 ) ) ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: County' //&
                        &   ' code is non-integer at entry', I
                        CALL M3MESG( MESG )
                    END IF
                    STA = STR2INT( SEGMENT( 1 ) ) * 1000
                    CNY = STR2INT( SEGMENT( 2 ) )
                    FIP = STA + CNY

!.....................  For region code with no spaces (CSSYYY)
                ELSE
                    FIP  = STR2INT( SEGMENT( 1 ) )
                    STA  = ( FIP / 1000 ) * 1000

                END IF

                RBUF = STR2REAL( SEGMENT( NSEG ) )

                WRITE( CFIP,FIPFMT ) FIP
                WRITE( CSTA,FIPFMT ) STA
                WRITE( CBUF,'(F12.6)' ) RBUF

                TSCC = SEGMENT( SCCCOL )
                TSCC = ADJUSTR( TSCC )
                L = LEN_TRIM( TSCC )

                K = K + 1
                INDXA( K ) = K
                FIPSA( K ) = FIP
                CFSCA( K ) = CFIP// TSCC( 1:L )
                CSSCA( K ) = CSTA// TSCC( 1:L )
                CRECA( K ) = CSTA// TSCC( 1:L )// CBUF

            END DO
            NRECS( N ) = K

            CALL M3MSG2( '    Sorting packet info...' )

!.................  Sort packet information by state, scc, factor
            CALL SORTIC( NRECS( N ), INDXA, CRECA )

            CALL M3MSG2( '    Reducing packet info...' )

!.................  Count duplicate state-SCC records and store first position
!                   of each with different packet entries
            LCREC = ' '
            LCFSC = ' '
            DO I = 1, NRECS( N )

                J = INDXA( I )
                FIP = FIPSA( J )
                CNY = FIP - (FIP/1000) * 1000    !  integer math

!.....................  If record is different from previous
                IF( CRECA( J ) .NE. LCREC ) THEN
                    RCOUNT( J ) = 1
                    RINDX ( J ) = J
                    IDXS1 = J

!.....................  Otherwise
                ELSE
                    RCOUNT( IDXS1 ) = RCOUNT( IDXS1 ) + 1
                    RINDX ( J ) = IDXS1

                END IF

!.....................  If default is already in file, mark record as
!                       not needing a new default.
                IF ( CNY .EQ. 0 ) FIXED( J ) = .TRUE.

!.....................  If State-SCC is different from the previous. Note that
!                       state-SCC will be different when RCOUNT( J ) has just
!                       been initialized to 1
                IF( CSSCA( J ) .NE. LCSSC ) THEN
                    IDXS2 = J
                    RMAX   ( J ) = RCOUNT( J )
                    DEFAULT( J ) = .TRUE.

!.....................  Otherwise, compare with current max for this state-scc
!                       combo and reset index if needed
!.....................  Do not reset counter or position of default if
!                       the state default was in the original inputs (county
!                       code is zero).
!.....................  Cannot assume that the original state default will
!                       be encountered first in the sorted list.
                ELSE IF ( CNY .EQ. 0 .OR.&
                &        ( .NOT. FIXED( IDXS2 )   .AND.&
                &         RCOUNT( IDXS1 ) .GT. RMAX( IDXS2 ) ) ) THEN

                    DEFAULT( IDXS2 ) = .FALSE.
                    IDXS2 = IDXS1
                    RMAX   ( IDXS2 ) = RCOUNT( IDXS1 )
                    DEFAULT( IDXS1 ) = .TRUE.

                END IF

                LCSSC = CSSCA( J )
                LCREC = CRECA( J )

            END DO

!.................  Resort so that the output order is easier to read.
            DO I = 1, NRECS( N )
                INDXA( I ) = I
            END DO

            CALL SORTIC( NRECS( N ), INDXA, CFSCA )

            CALL M3MSG2( '    Writing packet info...' )

!.................  Write packet header
            L = LEN_TRIM( PSTBUF( N ) )
            WRITE( ODEV,93000 ) PSTBUF( N )( 1:L )

!.................  Write defaults and all other records
            K1 = FIPLEN3 + 1
            K2 = FIPLEN3 + SCCLEN3
            K3 = K2 + 1
            DO I = 1, NRECS( N )

                J    = INDXA( I )
                K    = RINDX( J )
                K4   = LEN_TRIM( CRECA( J ) )
                FIP  = FIPSA( J )
                CSTA = CRECA( J )( 1:FIPLEN3 )
                TSCC = CRECA( J )( K1:K2 )
                CBUF = ADJUSTL( CRECA( J )( K3:K4 ) )

!.....................  If record is a state default, write out as such
                IF( DEFAULT( J ) ) THEN
                    WRITE( ODEV,93480 ) CSTA, TSCC, CBUF

                ELSE IF( .NOT. DEFAULT( K )        .AND.&
                &                RCOUNT( K ) .GT. 0       ) THEN
                    WRITE( ODEV,93490 ) FIP, TSCC, CBUF

                END IF

            END DO

!.................  Write packet tail
            L = LEN_TRIM( PENBUF( N ) )
            WRITE( ODEV,93000 ) PENBUF( N )( 1:L )

          CASE DEFAULT
            EFLAG = .TRUE.
!                WRITE( )

        END SELECT ! end packet type selection

    END DO         ! end loop on packets

!.........  Normal completion of program
    CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

!.........  Error completion

1011 MESG = 'Unexpected end of file'
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )


!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

93480 FORMAT( A, 1X, A, 1X, A )

93490 FORMAT( I6.6, 1X, A, 1X, A )

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10 ( A, :, I5, :, 2X ) )


END PROGRAM PKTREDUC


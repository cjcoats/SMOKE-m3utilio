
SUBROUTINE RDINVSRCS( FDEV, FNAME, NRAWBP, NRAWSRCS, ORLFLG )

    !***********************************************************************
    !  subroutine body starts at line 133
    !
    !  DESCRIPTION:
    !      This subroutine controls reading an ASCII inventory file for any source
    !      category from one of many formats.  It determines the format and
    !      calls the appropriate reader subroutines. It controls the looping
    !      through multiple files when a list-formatted file is used as input.
    !      This routine only reads the unique source characteristics from the inventories.
    !
    !  PRECONDITIONS REQUIRED:
    !      Input file unit FDEV opened
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !      Subroutines:
    !      Functions:
    !
    !  REVISION  HISTORY:
    !       Created 1/03 by C. Seppanen (based on rdinven.f)
    !       Updated 2/06 by B. Baek (adding wildfire format)
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
    !.......   This module is the inventory arrays
    USE MODSOURC, ONLY: CSOURCA, SRCIDA, NSTRECS, SRCSBYREC, RECIDX

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: NSRC, CATEGORY

    !.......  This module contains the lists of unique inventory information
    USE MODLISTS, ONLY: FILFMT, LSTSTR, FIREFLAG, FF10FLAG

    !.......  This module is for mobile-specific data
    USE MODMOBIL, ONLY: NSCCTBL, SCCTBL, SCCRVC, SCCMAPFLAG,    &
                        NSCCMAP, SCCMAPLIST, EXCLSCCFLAG

    IMPLICIT NONE

    !.......   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    INTEGER,      INTENT (IN) :: FDEV             ! unit no. of inv file
    CHARACTER(*), INTENT (IN) :: FNAME            ! logical name of file
    INTEGER,      INTENT(OUT) :: NRAWBP           ! no. of sources with pols/acts
    INTEGER,      INTENT(OUT) :: NRAWSRCs         ! no. of raw sources
    LOGICAL,      INTENT(OUT) :: ORLFLG           ! true: read ORL inventory

    !.......   EXTERNAL FUNCTIONS and their descriptions:

    INTEGER, EXTERNAL :: GETPID
    INTEGER, EXTERNAL :: GETFLINE
    INTEGER, EXTERNAL :: GETFORMT
    INTEGER, EXTERNAL :: GETINVYR
    INTEGER, EXTERNAL :: FIND1FIRST
    LOGICAL, EXTERNAL :: BLKORCMT
    LOGICAL, EXTERNAL :: CHKINT
    LOGICAL, EXTERNAL :: USEEXPGEO

    !.......   Local parameters
    INTEGER      , PARAMETER :: MXRECS   = 3000000      ! maximum records per iteration
    INTEGER      , PARAMETER :: NSCSEG   = 8            ! num. segments in scratch file
    INTEGER      , PARAMETER :: NSEG     = 70           ! maximum no of segments
    CHARACTER(16), PARAMETER :: PROGNAME =  'RDINVSRCS'     ! program name

    !.......   Local arrays
    CHARACTER(ALLLEN3+OBRLEN3) :: TMPCSOURC( MXRECS )       ! source information from inventory file(s)
    INTEGER                    :: TCSRCIDX ( MXRECS )       ! index for sorting source info
    INTEGER                    :: FRSNUMS  ( MXRECS,3 )     ! triplets of file, record, and source number
    CHARACTER(SRCLEN3) SCSEGMENT( NSCSEG )       ! segments from scratch file
    CHARACTER( 40 )    SEGMENT  ( NSEG )         ! segments of line

    INTEGER,            ALLOCATABLE:: CSRCIDX  ( : )        ! index for sorting CSOURCA

    !.......   File units and logical/physical names
    INTEGER         CDEV            !  scratch file

    !.......   Other local variables
    INTEGER         I, J, JJ, K, KK, K1, K2, L, NP, S     !  counters and indices
    INTEGER         L0, L1, L2, L3, L4, L5, L6, L7, L8, L9

    INTEGER         CSRC_LEN         !  length of source characteristics
    INTEGER         CURFMT           !  format of current inventory file
    INTEGER         CURFIL           !  current file from list formatted inventory
    INTEGER         IOS              !  i/o status
    INTEGER         INVFMT           !  inventory format code
    INTEGER         IREC             !  no. of records read
    INTEGER         ISTREC           !  no. of records stored
    INTEGER         IVT              !  vehicle type code
    INTEGER         LDEV             !  device no. for log file
    INTEGER         NSCC             !  tmp no of reference SCCs
    INTEGER         MXWARN           !  maximum number of warnings
    INTEGER         NLINE            !  number of lines in list format file
    INTEGER         NPOLPERLN        !  no. of pollutants per line of inventory file
    INTEGER         NRECPERLN        !  no. of records per line
    INTEGER      :: NWARN0= 0        !  current number of warnings
    INTEGER      :: NWARN1= 0        !  current number of warnings 1
    INTEGER      :: NWRLINE = 0      !  no. of lines in file writting to log
    INTEGER*4       PID              !  UNIX process ID at runtime
    INTEGER         ROAD             !  road class number
    INTEGER         RWT              !  roadway type
    INTEGER      :: TOTSRCS = 0      !  total number of sources
    INTEGER      :: TOTRECS = 0      !  total number of records
    INTEGER      :: NORSID = 0       !  no of Oris IDs under same EIS unit

    LOGICAL      :: ORSFLAG = .FALSE.     ! true: process multiple ORIS units under the same EIS unit
    LOGICAL      :: EFLAG   = .FALSE.     ! true: error occured
    LOGICAL      :: HDRFLAG               ! true: current line is part of header
    LOGICAL      :: LSTTIME = .FALSE.     ! true: last time through

    CHARACTER(FIPLEN3) CFIP        ! fips code
    CHARACTER(LNKLEN3) CLNK        ! link ID
    CHARACTER(VIDLEN3) CIVT        ! vehicle type ID
    CHARACTER(RWTLEN3) CROAD       ! road class no.
    CHARACTER(RWTLEN3) CRWT        ! roadway type
    CHARACTER(VTPLEN3) VTYPE       ! tmp vehicle type
    CHARACTER(RWTLEN3+VTPLEN3) CRVC        ! tmp roadway // vehicle type

    CHARACTER(PLTLEN3) FCID        ! facility/plant ID
    CHARACTER(CHRLEN3) PTID, NPTID        ! point ID
    CHARACTER(CHRLEN3) SKID        ! stack ID
    CHARACTER(CHRLEN3) SGID        ! segment ID
    CHARACTER(CHRLEN3) DVID        ! device ID
    CHARACTER(CHRLEN3) PRID        ! process ID

    CHARACTER(ORSLEN3) :: CORS = ' '      ! DOE plant ID
    CHARACTER(BLRLEN3) :: BLID = ' '      ! boiler ID

    CHARACTER(SCCLEN3) TSCC        ! scc code
    CHARACTER(ALLLEN3) TCSOURC     ! concatenated src (minus pollutant)
    CHARACTER(ALLLEN3+OBRLEN3) FCSOURC     ! concatenated src (minus pollutant)

    CHARACTER(10)      CREC        ! record number
    CHARACTER(4)       CFIL        ! file number
    CHARACTER(300)     BUFFER      ! tmp line buffer
    CHARACTER(300)     OUTLINE     ! line to write to scratch file
    CHARACTER(300)     INFILE      ! input file line buffer
    CHARACTER(500)     LINE        ! input file line buffer
    CHARACTER(300)     MESG        ! message buffer
    CHARACTER(20)      VIDFMT      ! vehicle type ID format
    CHARACTER(20)      RWTFMT      ! roadway type number format
    CHARACTER(1024)    TMPFILNAM      ! File name of tmp file

    CHARACTER(512)     PATHNM               ! path name for tmp file
    CHARACTER(300)     TENLINES( 10 )       ! first ten lines of inventory file

    !***********************************************************************
    !   begin body of subroutine RDINVSRCS

    !.......  Get log file number for reports
    LDEV = INIT3()

    TMPCSOURC = ' '     ! array
    TCSRCIDX  = 0       ! array
    FRSNUMS   = 0       ! array

    !.......  Get maximum number of warnings
    MXWARN = ENVINT( WARNSET, ' ', 100, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble WARNSET', 2 )
    END IF

    !.......  Handle multiple ORIS units under same EIS Unit
    MESG = 'Process multiple ORIS/Boiler units under the same source unit'
    ORSFLAG = ENVYN( 'PROCESS_MULT_ORIS_UNITS_YN', MESG, .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "PROCESS_MULT_ORIS_UNITS_YN"', 2 )
    END IF

    !.......  Get temporary directory location
    MESG = 'Path where temporary import file will be written'
    CALL ENVSTR( 'SMK_TMPDIR', MESG, '.', PATHNM, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK_TMPDIR"', 2 )
    END IF

    IF( IOS /= 0 ) THEN
        IF( NWARN0 < MXWARN ) THEN
            MESG = 'WARNING: Temporary input file will be ' //    &
                   'placed in executable directory because ' //    &
                   CRLF() // BLANK10 // 'SMK_TMPDIR environment '//    &
                   'variable is not set properly'
            CALL M3MSG2( MESG )
            NWARN0 = NWARN0 + 1
        END IF
    END IF

    !.......  Get process ID for using in tmp file name
    PID = GETPID()

    !.......  Build tmp file name and set environment variable to its value,
    !           so calling script can clean up tmp file if desired.
    WRITE( TMPFILNAM, '(A,I8)') TRIM( PATHNM )// '/import_tmp_', PID
    IF ( .NOT. SETENVVAR( 'SMKINVEN_TMPFILE', TMPFILNAM )) THEN
        MESG = 'Could not set environment variable for Smkinven '//    &
              'temporary file name:'// CRLF()// BLANK10// TMPFILNAM
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Initialize ORL (fire) flag to false
    ORLFLG = .FALSE.
    FIREFLAG = .FALSE.

    !.......  Create formats for mobile data
    IF( CATEGORY == 'MOBILE' ) THEN
        WRITE( VIDFMT, '("(I",I2.2,")")' ) VIDLEN3
        WRITE( RWTFMT, '("(I",I2.2,")")' ) RWTLEN3
    END IF

    !.......  Determine file format of inventory file
    INVFMT = GETFORMT( FDEV, -1 )

    !.......  If SMOKE list format, read file and check file for formats.
    !           NOTE- LSTFMT defined in EMCNST3.h90
    IF( INVFMT == LSTFMT ) THEN

        !.......  Generate message for GETFLINE and RDLINES calls
        MESG = TRIM( CATEGORY ) // ' inventory file, ' // TRIM( FNAME ) // ', in list format'

        !.......  Get number of lines of inventory files in list format
        NLINE = GETFLINE( FDEV, MESG )

        !.......  Allocate memory for storing contents of list-formatted file
        ALLOCATE( FILFMT( NLINE ),    &
                  LSTSTR( NLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LSTSTR', PROGNAME )

        FILFMT = -1          ! array
        LSTSTR = ' '         ! array

        !.......  Store lines of PTINV file
        CALL RDLINES( FDEV, MESG, NLINE, LSTSTR )

        !.......  Reset number of lines to remove blanks
        !               (RDLINES does not store blank lines)
        DO I = 1, NLINE
            IF( LSTSTR( I ) == ' ' ) THEN
                NLINE = I - 1
                EXIT
            END IF
        END DO

        !.......  Check the format of the list-formatted inventory file and
        !               return the code for the type of files it contains
        CALL CHKLSTFL( NLINE, FNAME, LSTSTR, FILFMT )

        !.......  Close original inventory file (will reuse device number for individual files)
        CLOSE( FDEV )

    !.......  If not list format, then set FILFMT to the type of file (IDA,EPS)
    ELSE

        NLINE = 1
        ALLOCATE( FILFMT( NLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FILFMT', PROGNAME )
        FILFMT = INVFMT

    END IF

    !.......  Set default inventory characteristics (declared in MODINFO) used
    !           by the IDA and EPS formats, including NPPOL
    DO I = 1, NLINE
        IF( FILFMT( I ) > 0 ) THEN
            CALL INITINFO( FILFMT( I ) )
            EXIT
        END IF
    END DO

    CURFIL = 1

    !.......  If file is list format, open first file
    IF( INVFMT == LSTFMT ) THEN

        LINE = LSTSTR( CURFIL )

        !.......  Skip #LIST lines (must be first)
        IF( INDEX( LINE, 'LIST' ) > 0 ) THEN
            CURFIL = CURFIL + 1
            LINE = LSTSTR( CURFIL )
        END IF

        !.......  Check for inventory year packet
        IF( GETINVYR( LINE ) > 0 ) THEN
            CURFIL = CURFIL + 1          ! move to next file in list
        END IF

        !.......  Make sure there are more files
        IF( CURFIL > NLINE ) THEN
            MESG = 'No individual inventory files in list file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        !.......  Store path of file name
        INFILE = LSTSTR( CURFIL )

        !.......  Open current file
        OPEN( FDEV, FILE=INFILE, STATUS='OLD', IOSTAT=IOS )

        !.......  Check for errors while opening file
        IF( IOS /= 0 ) THEN

            WRITE( MESG,94010 ) 'Problem at line ', CURFIL, 'of ' //    &
               TRIM( FNAME ) // '.' // ' Could not open file:' //    &
               CRLF() // BLANK5 // TRIM( INFILE )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSE
            WRITE( MESG,94010 ) 'Successful OPEN for ' //    &
               'inventory file:' // CRLF() // BLANK5 //    &
               TRIM( INFILE )
            CALL M3MSG2( MESG )

        END IF

        CURFMT = FILFMT( CURFIL )

        !.......  Set default inventory characteristics that depend on file format
        CALL INITINFO( CURFMT )

    ELSE

        !.......  If not list format, set current format to inventory format
        CURFMT = INVFMT

    END IF

    IREC = 0       ! current record number

    !.......  Open scratch file for writing record numbers
    CDEV = JUNIT()
    OPEN( CDEV, FILE=TMPFILNAM, IOSTAT=IOS, STATUS='NEW' )

    IF( IOS /= 0 ) THEN
        MESG = 'Could not open temporary import file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Loop over files and multiples of MXRECS
    DO

        !.......  Reset counters
        S = 0               ! source number
        ISTREC= 0           ! number of records stored

        !.......  Loop through records in current file
        DO
            IF( ISTREC == MXRECS ) EXIT

            READ( FDEV, 93000, IOSTAT=IOS ) LINE

            IREC = IREC + 1
            IF( IOS > 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'I/O error', IOS,    &
                   'reading inventory file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

            !.......  Check if we've reached the end of the file
            IF( IOS < 0 ) THEN

                !.......  If list format, try to open next file
                IF( INVFMT == LSTFMT ) THEN

                    !.......  Close current file and reset counter
                    CLOSE( FDEV )
                    IREC = 0

                    !.......  Advance to next file
                    CURFIL = CURFIL + 1

                    !.......  Check if there are more files to read
                    IF( CURFIL <= NLINE ) THEN
                        LINE = LSTSTR( CURFIL )

                        !.........  Check for #LIST line
                        IF( INDEX( LINE, 'LIST' ) > 0 ) THEN
                            CURFIL = CURFIL + 1                          ! move to next file in list
                        END IF

                        !.........  Make sure current line is not INVYEAR packet
                        IF( GETINVYR( LINE ) > 0 ) THEN
                            CURFIL = CURFIL + 1                          ! move to next file in list
                        END IF

                        !.........  Make sure there are still files to read
                        IF( CURFIL > NLINE ) THEN
                            LSTTIME = .TRUE.
                            EXIT
                        END IF

                        INFILE = LSTSTR( CURFIL )

                        OPEN( FDEV, FILE=INFILE, STATUS='OLD',    &
                              IOSTAT=IOS )

                        !.........  Check for errors while opening file
                        IF( IOS /= 0 ) THEN

                            WRITE( MESG,94010 ) 'Problem at line ',    &
                               CURFIL, 'of ' // TRIM( FNAME ) //    &
                               '.' // ' Could not open file:' //    &
                               CRLF() // BLANK5 // TRIM( INFILE )
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                        ELSE
                            WRITE( MESG,94010 )    &
                              'Successful OPEN for ' //    &
                              'inventory file(s):' // CRLF() //    &
                              BLANK5 // TRIM( INFILE )
                            CALL M3MSG2( MESG )

                        END IF

                        !.........  Set default inventory characteristics that depend on file format
                        CALL INITINFO( FILFMT( CURFIL ) )
                        CURFMT = FILFMT( CURFIL )
                        NWRLINE = 0

                        !.........  Skip back to the beginning of the loop
                        CYCLE

                    !.......  Otherwise, no more files to read, so exit
                    ELSE
                        LSTTIME = .TRUE.
                        EXIT
                    END IF

                !.......  Otherwise, not a list file, so exit
                ELSE
                    LSTTIME = .TRUE.
                    EXIT
                END IF

            END IF               ! end check for end of file

            !.......  Skip blank lines
            IF( LINE == ' ' ) CYCLE

            !.......  Check if format works with expanded geographic codes
            IF( USEEXPGEO() ) THEN
                IF( CURFMT == ORLFMT .OR. CURFMT == ORLNPFMT .OR. CURFMT == ORLFIREFMT ) THEN
                    MESG = 'ERROR: Expanded geographic codes are only ' //    &
                           'supported for inventories in FF10 format.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
            END IF

            !.......  Process line depending on file format and source category
            SELECT CASE( CURFMT )

              CASE( MEDSFMT )
                SELECT CASE( CATEGORY )
                  CASE( 'POINT' )           ! used for pregridded MEDS format inv
                    CALL RDSRCMEDSPT( LINE, CFIP, FCID, PTID, SKID,    &
                                     SGID, TSCC, NPOLPERLN, HDRFLAG,   &
                                     EFLAG )
                END SELECT

              CASE( FF10FMT )
                ORLFLG = .TRUE.
                FF10FLAG = .TRUE.

                SELECT CASE( CATEGORY )
                  CASE( 'AREA' )           ! used for nonroad only
                    CALL RDSRCFF10AR( LINE, CFIP, TSCC, NPOLPERLN, HDRFLAG, EFLAG )
                  CASE( 'MOBILE' )
                    CALL RDSRCFF10MB( LINE, CFIP, CLNK, TSCC, NPOLPERLN, HDRFLAG, EFLAG )
                  CASE( 'POINT' )
                    CALL RDSRCFF10PT( LINE, CFIP, FCID, PTID, SKID,         &
                                     SGID, TSCC, CORS, BLID, NPOLPERLN,     &
                                     HDRFLAG, EFLAG )
                END SELECT

              CASE( ORLFMT )
                ORLFLG = .TRUE.

                SELECT CASE( CATEGORY )
                  CASE( 'AREA' )           ! used for nonroad only
                    CALL RDSRCORLAR( LINE, CFIP, TSCC, NPOLPERLN, HDRFLAG, EFLAG )
                  CASE( 'MOBILE' )
                    CALL RDSRCORLMB( LINE, CFIP, CLNK, TSCC, NPOLPERLN, HDRFLAG, EFLAG )
                  CASE( 'POINT' )
                    CALL RDSRCORLPT( LINE, CFIP, FCID, PTID, SKID,    &
                                     SGID, TSCC, CORS, BLID, NPOLPERLN,    &
                                     HDRFLAG, EFLAG )
                END SELECT

              CASE( ORLNPFMT )
                ORLFLG = .TRUE.

                CALL RDSRCORLNP( LINE, CFIP, TSCC, NPOLPERLN, HDRFLAG, EFLAG )

              CASE( ORLFIREFMT )
                ORLFLG   = .TRUE.
                FIREFLAG = .TRUE.

                CALL RDSRCORLFR( LINE, CFIP, FCID, PTID, SKID,    &
                                 SGID, TSCC, NPOLPERLN, HDRFLAG,  &
                                 EFLAG )

              CASE DEFAULT
                WRITE( MESG, 94010 ) 'Routine rdinvsrc.f not '//  &
                       'expecting to read file of format code',   &
                        CURFMT
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END SELECT

            !.......  Check for header lines
            IF( HDRFLAG ) THEN
                CYCLE
            END IF

            !.......  SCC mapping loop : Mobile activity inventory use only.
            KK = 0
            NSCC = 0
            IF( CATEGORY == 'MOBILE' ) THEN
                IF( SCCMAPFLAG ) THEN
                    CALL PADZERO( TSCC )
                    KK   = INDEX1( TSCC, NSCCMAP, SCCMAPLIST( :,1 ) )
                    IF( KK > 0 ) THEN
                        NSCC = STR2INT( SCCMAPLIST( KK,3 ) )
                    ELSE
                        IF( EXCLSCCFLAG ) THEN
                            MESG = 'WARNING: Dropping SCC "' // TRIM( TSCC ) //    &
                                '" not listed in SCCXREF file'
                            CALL M3MESG( MESG )
                            CYCLE      ! skip SCC not found in SCCXREF file
                        END IF
                    END IF
                END IF
            END IF

            !.......  loop over mapped SCC
            DO JJ = 0, NSCC

                IF( JJ > 0 .AND. KK > 0 ) IREC = IREC + 1                 ! increment no of records by reference SCCs
                IF( SCCMAPFLAG .AND. KK > 0 ) TSCC = SCCMAPLIST( KK+JJ,2 )

                !.......  Write first ten lines of inventory to log file
                IF( NWRLINE < 10 .AND. .NOT. FIREFLAG ) THEN
                    NWRLINE = NWRLINE + 1
                    TENLINES( NWRLINE ) = BLANK10 // TRIM( LINE )

                    IF( NWRLINE == 10 ) THEN
                        MESG = BLANK10 // 'First 10 lines of current inventory:'
                        WRITE( LDEV, '(A)' ) TRIM( MESG )

                        DO I = 1,NWRLINE
                            WRITE( LDEV, '(A)' ) TRIM( TENLINES( I ) )
                        END DO
                    END IF
                END IF

                !.......  Check that source characteristics are correct
                !.......  Make sure some emissions are kept for this source
                IF( NPOLPERLN == 0 ) THEN
                    IF( NWARN0 < MXWARN ) THEN
                        WRITE( MESG,94010 ) 'WARNING: No kept '//    &
                               'pollutants found at line', IREC, '. ' //    &
                               'The source will be dropped.'
                        CALL M3MESG( MESG )
                        NWARN0 = NWARN0 + 1
                    END IF
                    CYCLE
                END IF

                IF( .NOT. USEEXPGEO() .AND. .NOT. CHKINT( CFIP ) ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: State and/or ' //    &
                           'county code is non-integer at line', IREC
                    IF( NWARN1 < MXWARN ) CALL M3MESG( MESG )
                    NWARN1 = NWARN1 + 1
                END IF

                IF( .NOT. USEEXPGEO() .AND.    &
                    ( CFIP( FIPEXPLEN3+2:FIPEXPLEN3+3 ) == '00' .OR.    &
                      CFIP( FIPEXPLEN3+4:FIPEXPLEN3+6 ) == '000'     ) ) THEN
                    WRITE( MESG,94010 ) 'WARNING: State and/or ' //    &
                           'county code is zero (missing) at line', IREC
                    IF( NWARN1 < MXWARN ) CALL M3MESG( MESG )
                    NWARN1 = NWARN1 + 1
                END IF

                !.......  Check source specific characteristics
                IF( CATEGORY == 'AREA' ) THEN

                    !.......  Make sure SCC is at least 8 characters long
                    IF( LEN_TRIM( TSCC ) < 8 ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: SCC code must ' //    &
                            'be at least 8 characters long at line',    &
                            IREC
                        CALL M3MESG( MESG )
                    END IF

                END IF

                IF( CATEGORY == 'MOBILE' ) THEN

                !.......  Set vehicle type and road class
                    CALL PADZERO( TSCC )

                    IVT = STR2INT( TSCC( 13:14 ) )
                    RWT = STR2INT( TSCC( 15:16 ) )

                ELSE IF( CATEGORY == 'POINT' ) THEN

                    IF( CURFMT == ORLFMT .OR. CURFMT == MEDSFMT .OR.    &
                        CURFMT == ORLFIREFMT .OR. CURFMT == FF10FMT ) THEN

                        !.......  Make sure SCC is at least 8 characters long
                        IF( LEN_TRIM( TSCC ) < 8 ) THEN
                            IF( NWARN0 < MXWARN ) THEN
                                WRITE( MESG,94010 ) 'WARNING: SCC ' //    &
                                   'code is less than 8 characters ' //    &
                                   'long at line', IREC, '. Adding ' //    &
                                   'trailing zeros.'
                                CALL M3MESG( MESG )
                                NWARN0 = NWARN0 + 1
                            END IF

                            DO I = LEN_TRIM( TSCC )+1, 8
                                TSCC( I:I ) = '0'
                            END DO
                        END IF

                        !.......  Make sure we have a facility/plant ID
                        IF( FCID == ' ' ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 'ERROR: Missing ' //    &
                                'plant/facility ID code at line', IREC
                            CALL M3MESG( MESG )
                        ELSE IF( LEN_TRIM( FCID ) > PLTLEN3 ) THEN
                            IF( NWARN0 < MXWARN ) THEN
                                WRITE( MESG,94010 ) 'WARNING: Facility ' //    &
                                   'ID is longer than maximum 20 characters ' //    &
                                   'long at line', IREC
                                CALL M3MESG( MESG )
                                NWARN0 = NWARN0 + 1
                            END IF
                        END IF

                    END IF

                END IF

                !.......  Skip rest of loop if an error has occured
                IF( EFLAG ) CYCLE

                !.......  Build concatenated source information
                SELECT CASE( CATEGORY )
                  CASE( 'AREA' )
                    CALL PADZERO( TSCC )

                    CALL BLDCSRC( CFIP, TSCC, CHRBLNK3, CHRBLNK3,    &
                                  CHRBLNK3, CHRBLNK3, CHRBLNK3,      &
                                  CHRBLNK3, TCSOURC )
                  CASE( 'MOBILE' )
                    CALL FLTRNEG( CLNK )
                    WRITE( CRWT,RWTFMT ) RWT
                    WRITE( CIVT,VIDFMT ) IVT

                    CALL BLDCSRC( CFIP, CRWT, CLNK, CIVT, TSCC,    &
                                  CHRBLNK3, CHRBLNK3, CHRBLNK3,    &
                                  TCSOURC )

                  CASE( 'POINT' )
                    CALL PADZERO( TSCC )

                    CALL BLDCSRC( CFIP, FCID, PTID, SKID, SGID,    &
                                  TSCC, CHRBLNK3, CHRBLNK3,        &
                                  TCSOURC )
                END SELECT

                CSRC_LEN = LEN_TRIM( TCSOURC )

                !.......  Append DOE plant ORIS and Boiler IDs to define src
                FCSOURC = TCSOURC // CORS // BLID

                !.......  Store source info on first time through
                IF( S == 0 ) THEN
                    S = S + 1
                    TCSRCIDX ( S ) = S
                    TMPCSOURC( S ) = FCSOURC
                END IF

                !.......  On subsequent passes, only store source info
                !                   if it does not match previous source
                IF( FCSOURC /= TMPCSOURC( S ) ) THEN
                    S = S + 1
                    TCSRCIDX ( S ) = S
                    TMPCSOURC( S ) = FCSOURC
                END IF

                !.......  Store current source number for this record
                ISTREC = ISTREC + 1
                FRSNUMS( ISTREC,1 ) = CURFIL
                FRSNUMS( ISTREC,2 ) = IREC
                FRSNUMS( ISTREC,3 ) = S

                !.......  Update total number of sources with pollutants
                NRAWBP = NRAWBP + NPOLPERLN

            END DO      ! loop through SCCMAPLIST (NSCC)

        END DO          ! loop through MXRECS lines

        !.......  Abort if there was a reading error
        IF( EFLAG ) THEN
            MESG = 'Error reading raw inventory file ' // TRIM( FNAME )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        !.......  Update total number of stored records
        TOTRECS = TOTRECS + ISTREC

        !.......  Sort source info
        CALL SORTIC( S, TCSRCIDX, TMPCSOURC )

        !.......  Write source info and record numbers to file
        DO I = 1, S

            J = TCSRCIDX( I )

            !.......  On first time through, set up first line and counter
            IF( I == 1 ) THEN
                OUTLINE = TMPCSOURC( J )( 1:CSRC_LEN )
                TCSOURC = TMPCSOURC( J )( 1:CSRC_LEN )
                FCSOURC = TMPCSOURC( J )
                NORSID = 0
                NRECPERLN = 0
                TOTSRCS = TOTSRCS + 1
            END IF

            !.......  If current source does not match previous, write old output line
            !                   and start new line
            IF( TMPCSOURC( J )( 1:CSRC_LEN ) /= FCSOURC( 1:CSRC_LEN) ) THEN
                WRITE( CDEV, '(A)' ) TRIM( OUTLINE )
                OUTLINE = TMPCSOURC( J )( 1:CSRC_LEN )
                TCSOURC = TMPCSOURC( J )( 1:CSRC_LEN )
                FCSOURC = TMPCSOURC( J )
                NORSID = 0                             ! reset number of oris/boiler IDs
                NRECPERLN = 0                          ! reset number of read records
                TOTSRCS = TOTSRCS + 1                  ! increment total number of sources

            !.......  Added new source when there are more than one
            !                       ORIS and Boilers under same plant ID
            ELSE IF( ORSFLAG ) THEN
                IF( TMPCSOURC( J ) /= FCSOURC ) THEN
                    NORSID = NORSID + 1
                    CFIP = TMPCSOURC( J )( PTBEGL3( 1 ):PTENDL3( 1 ) )
                    FCID = TMPCSOURC( J )( PTBEGL3( 2 ):PTENDL3( 2 ) )
                    PTID = TMPCSOURC( J )( PTBEGL3( 3 ):PTENDL3( 3 ) )
                    SKID = TMPCSOURC( J )( PTBEGL3( 4 ):PTENDL3( 4 ) )
                    SGID = TMPCSOURC( J )( PTBEGL3( 5 ):PTENDL3( 5 ) )
                    TSCC = TMPCSOURC( J )( PTBEGL3( 6 ):PTENDL3( 6 ) )

                    !.......  Update Unit ID with ## when multiple oris IDs
                    PTID = ADJUSTL( PTID )
                    L1 = LEN_TRIM( PTID )
                    IF( L1 > CHRLEN3-3 ) L1 = CHRLEN3 - 3
                    WRITE( NPTID, '( A,I2.2)' ) PTID( 1:L1 )//'_',NORSID

                    !.......  Warning msg for duplicate sources
                    CALL FMTCSRC( TCSOURC, 7, BUFFER, L2 )
                    MESG = 'WARNING: Multiple ORIS/Boiler IDs under same point source'//    &
                           CRLF() //BLANK10// BUFFER( 1:L2 ) //    &
                           CRLF() //BLANK10// 'Renamed Release Point ID "'//TRIM(PTID) //    &
                           '" to "'// TRIM(NPTID)//'"'
                    CALL M3MSG2( MESG )

                    CALL BLDCSRC( CFIP, FCID, NPTID, SKID, SGID,    &
                                  TSCC, CHRBLNK3, CHRBLNK3,    &
                                  TCSOURC )

                    WRITE( CDEV, '(A)' ) TRIM( OUTLINE )
                    OUTLINE = TCSOURC( 1:CSRC_LEN )
                    FCSOURC = TMPCSOURC( J )               ! preserve FCSOURC before it gest updated/modified
                    NRECPERLN = 0                          ! reset number of read records
                    TOTSRCS = TOTSRCS + 1                  ! increment total number of sources
                END IF

            END IF

            !.......  Find source number in records array
            K = FIND1FIRST( J, ISTREC, FRSNUMS( :,3 ) )

            IF( K <= 0 ) THEN
                WRITE( MESG,94010 ) 'INTERNAL ERROR: Could not ' //    &
                   'find source number', J, 'in file and record ' //    &
                    'number array'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            !.......  Loop over records corresponding to current source
            DO

                !.......  Make sure not to go outside the array
                IF( K > ISTREC ) EXIT

                IF( FRSNUMS( K,3 ) == J ) THEN

                    !.......  If already read NSCSEG-1 records, write line with continuation
                    !                           character and start new line
                    IF( NRECPERLN == NSCSEG-1 ) THEN
                        OUTLINE = TRIM( OUTLINE ) // ' \'
                        WRITE( CDEV, '(A)' ) TRIM( OUTLINE )
                        OUTLINE = TCSOURC( 1:CSRC_LEN )
                        NRECPERLN = 0
                    END IF

                    WRITE( CFIL, '(I4)' ) FRSNUMS( K,1 )
                    WRITE( CREC, '(I10)' ) FRSNUMS( K,2 )

                    !.......  If writing the first record, can't use trim otherwise
                    !                           will lose any blank source characteristics
                    IF( NRECPERLN == 0 ) THEN
                        OUTLINE = OUTLINE( 1:CSRC_LEN ) // ' ' //    &
                                  TRIM( ADJUSTL( CFIL ) ) // '/' //    &
                                  TRIM( ADJUSTL( CREC ) )
                    ELSE
                        OUTLINE = TRIM( OUTLINE ) // ' ' //    &
                                  TRIM( ADJUSTL( CFIL ) ) // '/' //    &
                                  TRIM( ADJUSTL( CREC ) )
                    END IF

                    NRECPERLN = NRECPERLN + 1
                    K = K + 1
                ELSE
                    EXIT
                END IF
            END DO                ! loop over records

            !.......  If last source, write final line
            IF( I == S ) THEN
                WRITE( CDEV, '(A)' ) TRIM( OUTLINE )
            END IF
        END DO           ! loop to write sources to file

        !.......  Check if this is last time
        IF( LSTTIME ) EXIT

    END DO

    !.......  Rewind scratch file
    REWIND( CDEV )

    !.......  Allocate memory to read complete scratch file
    ALLOCATE( CSRCIDX( TOTSRCS ),    &
              CSOURCA( TOTSRCS ),    &
            SRCSBYREC( TOTRECS,3 ),  &
               RECIDX( TOTRECS ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CSRCIDX...RECIDX', PROGNAME )

    S = 0
    ISTREC = 0
    TCSOURC = ' '

    !.......  Read source info from scratch file
    DO

        READ( CDEV, 93000, IOSTAT=IOS ) LINE

        !.......  Check for I/O errors
        IF( IOS > 0 ) THEN
            WRITE( MESG, 94010 ) 'I/O error', IOS, 'reading scratch file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        !.......  Check for end of file
        IF( IOS < 0 ) EXIT

        !.......  Parse line after source info into segments
        BUFFER = LINE( CSRC_LEN+1:LEN_TRIM( LINE ) )
        CALL PARSLINE( BUFFER, NSCSEG, SCSEGMENT )

        !.......  Check if this is a continuation of a previous line
        IF( LINE( 1:CSRC_LEN ) /= TCSOURC ) THEN
            S = S + 1

        !.......  Store information from line
            CSRCIDX( S ) = S
            CSOURCA( S ) = LINE( 1:CSRC_LEN )
            TCSOURC = CSOURCA( S )
        END IF

        !.......  Loop through segments (file and record numbers)
        DO I = 1, NSCSEG-1

            !.......  Exit if segment is blank (reached end of line)
            IF( SCSEGMENT( I ) == ' ' ) EXIT

            !.......  Increment record counter and initialize sorting array
            ISTREC = ISTREC + 1
            RECIDX( ISTREC ) = ISTREC

            !.......  Find location of separator
            K = INDEX( SCSEGMENT( I ), '/' )
            L = LEN_TRIM( SCSEGMENT( I ) )

            !.......  Store file number
            SRCSBYREC( ISTREC,1 ) = STR2INT( SCSEGMENT( I )( 1:K-1 ) )
            SRCSBYREC( ISTREC,2 ) = STR2INT( SCSEGMENT( I )( K+1:L ) )
            SRCSBYREC( ISTREC,3 ) = S

        END DO

    END DO

    NRAWSRCS = S
    NSTRECS  = ISTREC

    !.......  Sort sources by record array by file number then record number
    CALL M3MESG( 'Sorting sources by file and line number...' )

    CALL SORTI2( NSTRECS, RECIDX, SRCSBYREC( :,1 ),    &
                 SRCSBYREC( :,2 ) )

    !.......  Sort inventory sources if needed (only if have more than MXRECS values)
    IF( NSTRECS > MXRECS ) THEN
        CALL M3MESG( 'Sorting sources by characteristics...' )

        CALL SORTIC( NRAWSRCS, CSRCIDX, CSOURCA )
    END IF

    !.......  Allocate memory for source numbers
    ALLOCATE( SRCIDA( NRAWSRCS ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SRCIDA', PROGNAME )

    !.......  Loop through sources to determine source IDs
    TCSOURC = ' '
    NSRC = 0

    DO I = 1, NRAWSRCS
        J = CSRCIDX( I )

        IF( CSOURCA( J ) /= TCSOURC ) THEN
            NSRC = NSRC + 1
            TCSOURC = CSOURCA( J )
        END IF

        SRCIDA( J ) = NSRC
    END DO

    !.......  Deallocate arrays that are no longer needed
    DEALLOCATE( CSRCIDX )

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

93000 FORMAT( A )

94010 FORMAT( 10( A, :, I8, :, 1X ) )

94060 FORMAT( 10( A, :, E10.3, :, 1X ) )

94070 FORMAT( I3, A1, I8 )

END SUBROUTINE RDINVSRCS


SUBROUTINE RDREPIN( NSLIN, NSSIN, RDEV, SDEV, GDEV, PDEV, TDEV,     &
                    EDEV, YDEV, NDEV, NIDEV, NPDEV, NMDEV, NNDEV,   &
                    NODEV, ADEV, ENAME, CUNAME, GNAME, LNAME,       &
                    PRNAME, SLNAME, SSNAME, NX, IX, CX, SSMAT, SLMAT )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      The RDREPIN routine reads in the SMOKE intermediate files and other
    !      files needed for generating the reports.
    !
    !  PRECONDITIONS REQUIRED:
    !    REPCONFIG file is opened
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
    !***********************************************************************
    USE M3UTILIO

    !.......   MODULES for public variables
    !.......   This module is the inventory arrays
    USE MODSOURC, ONLY: SRGID, CSOURC, CDOM, CWEK, CMON,            &
                        CMND, CTUE, CWED, CTHU, CFRI, CSAT,         &
                        CSUN, CMET, CIFIP, STKHT, STKDM, STKTK,     &
                        STKVE, CINTGR

    !.......  This module contains Smkreport-specific settings
    USE MODREPRT, ONLY: NMATX, AFLAG, NREPORT, GFLAG, GSFLAG,       &
                        CUFLAG, SLFLAG, SSFLAG, TSFLAG, PSFLAG,     &
                        YFLAG, VFLAG, LFLAG, NFLAG, SDATE, STIME,   &
                        NSTEPS, TSTEP, LOC_BEGP, LOC_ENDP,          &
                        ASCREC, PRRPTFLG, PRFLAG, MINC,             &
                        NSPCPOL, SPCPOL, NMAJOR, NPING, ALLRPT,     &
                        STKX, STKY, LSPCPOL, NIFLAG, NMFLAG, NNFLAG,&
                        NOFLAG, SDFLAG

    !.......  This module contains report arrays for each output bin
    USE MODREPBN, ONLY: NSVARS, SPCOUT

    !.......  This module contains the control packet data and control matrices
    USE MODCNTRL, ONLY: NVPROJ, NVCMULT, PRMAT, ACUMATX,    &
                        PNAMPROJ, PNAMMULT

    !.......  This module contains arrays for plume-in-grid and major sources
    USE MODELEV, ONLY: LMAJOR, LPING, GROUPID

    !.......  This module contains the lists of unique source characteristics
    USE MODLISTS, ONLY: NINVIFIP, INVCFIP

    !.......  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: NGRID

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: NSRC, NCHARS, CATEGORY, CRL, JSCC, ATTRUNIT, &
                       SC_BEGP, SC_ENDP, NIPPA

    !.......  This module is required for the FileSetAPI
    USE MODFILESET

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters
    INCLUDE 'SETDECL.h90'       !  FileSetAPI function declarations

    !.......   SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN) :: NSLIN      ! no. mass spec input vars
    INTEGER     , INTENT (IN) :: NSSIN      ! no. mass spec input vars
    INTEGER     , INTENT (IN) :: RDEV(3)    ! control report files
    INTEGER     , INTENT (IN) :: SDEV       ! unit no.: ASCII inven file
    INTEGER     , INTENT (IN) :: GDEV       ! unit no.: gridding supplemental
    INTEGER     , INTENT (IN) :: PDEV       ! unit no.: speciation supplemental
    INTEGER     , INTENT (IN) :: TDEV       ! unit no.: temporal supplemental
    INTEGER     , INTENT (IN) :: EDEV       ! unit no.: elevated ID file (PELV)
    INTEGER     , INTENT (IN) :: YDEV       ! unit no.: cy/st/co file
    INTEGER     , INTENT (IN) :: NDEV       ! unit no.: SCC descriptions
    INTEGER     , INTENT (IN) :: NIDEV      ! unit no.: SIC descriptions
    INTEGER     , INTENT (IN) :: NPDEV      ! unit no.: speciation GSPRO desc.
    INTEGER     , INTENT (IN) :: NMDEV      ! unit no.: MACT descriptions
    INTEGER     , INTENT (IN) :: NNDEV      ! unit no.: NAICS descriptions
    INTEGER     , INTENT (IN) :: NODEV      ! unit no.: ORIS descriptions
    INTEGER     , INTENT (IN) :: ADEV       ! unit no.: ASCII elevated file
    CHARACTER(*), INTENT (IN) :: ENAME      ! name for I/O API inven input
    CHARACTER(*), INTENT (IN) :: CUNAME     ! mulitplicative control matrix name
    CHARACTER(*), INTENT (IN) :: GNAME      ! gridding matrix name
    CHARACTER(*), INTENT (IN) :: LNAME      ! layer fractions file name
    CHARACTER(*), INTENT (IN) :: PRNAME     ! projection matrix name
    CHARACTER(*), INTENT (IN) :: SLNAME     ! speciation matrix name
    CHARACTER(*), INTENT (IN) :: SSNAME     ! speciation matrix name
    INTEGER     , INTENT(OUT) :: NX( NGRID )     ! no. srcs per cell
    INTEGER     , INTENT(OUT) :: IX( NMATX )     ! src IDs
    REAL        , INTENT(OUT) :: CX( NMATX )     ! gridding coefficients
    REAL        , INTENT(OUT) :: SSMAT( NSRC, NSSIN )     ! mass spec coefs
    REAL        , INTENT(OUT) :: SLMAT( NSRC, NSLIN )     ! mole spec coefs

    !....... EXTERNAL FUNCTIONS and their descriptions:
    LOGICAL,EXTERNAL :: CHKINT
    LOGICAL,EXTERNAL :: CHKREAL
    INTEGER,EXTERNAL :: GETFLINE
    INTEGER,EXTERNAL :: GETNLIST
    LOGICAL,EXTERNAL :: USEEXPGEO

    !.......  Local allocatable arrays
    REAL        LFRAC1L( NSRC )      ! 1st-layer fraction

    !.......  Array that contains the names of the inventory variables needed for
    !         this program
    CHARACTER(NAMLEN3) IVARNAMS( MXINVARR )

    !.......  For parsing lines
    CHARACTER(64)         SEGMENT( 10 )
    CHARACTER(CHRLEN3) :: CHARS  ( 5 )       ! tmp plant characteristics

    !.......   Local variables that depend on module variables
    INTEGER    SWIDTH( NCHARS )

    !.......   Other local variables
    INTEGER          I, J, K, L, L1, L2, N, V, S, T     ! counters and indices

    INTEGER          DIU                    ! tmp diurnal profile number
    INTEGER          IOS                    ! i/o status
    INTEGER          IREC                   ! tmp record number
    INTEGER       :: JDATE = 0              ! Julian date
    INTEGER       :: JTIME = 0              ! time (HHMMSS)
    INTEGER          MON                    ! tmp monthly profile number
    INTEGER       :: NINVARR = 0            ! no. actual inventory inputs
    INTEGER          NREPLIN                ! no. lines in input report
    INTEGER          NS                     ! tmp no. strings on line
    INTEGER          NV                     ! tmp no. variables in temporal suplm
    INTEGER          NSTK                   ! no. of stacks in ASCII elevated file
    INTEGER          COLS                   ! no. of cols in grid
    INTEGER          ROWS                   ! no. of rows in grid
    INTEGER       :: SRGID1                 ! tmp primary surrogate IDs
    INTEGER       :: SRGID2                 ! tmp fallback surrogate IDs
    INTEGER          WEK                    ! tmp weekly profile number

    REAL             XO, YO                 ! x and y origins
    REAL             XC, YC                 ! x and y cell widths
    REAL             XL, YL                 ! x and y stack locations

    LOGICAL       :: LRDREGN = .FALSE.      !  true: read region code
    LOGICAL       :: EFLAG   = .FALSE.      !  true: error found
    LOGICAL       :: MSGFLAG = .FALSE.      !  true: don't repeat message
    LOGICAL       :: LTMP    = .FALSE.      !  true: temporary logical
    LOGICAL, SAVE :: FIRSTIME =.FALSE.      !  true: first time

    CHARACTER(TMPLEN3) TTYP, TBUF       !  temporal profile entry type
    CHARACTER(TMPLEN3) FIRSTTYP         !  1st temporal profile entry type
    CHARACTER(16) ::   BNAME = ' '      !  name buffer
    CHARACTER(50)      BUFFER           !  string buffer
    CHARACTER(256)     LINE             !  input line
    CHARACTER(256)     MESG             !  message buffer
    CHARACTER(NAMLEN3) CBUF             !  tmp pollutant name
    CHARACTER(MXDLEN3) DBUF             !  tmp variable name
    CHARACTER(FIPLEN3) CFIP             !  tmp ASCII FIPS code
    CHARACTER(LNKLEN3) CLNK             !  tmp link code
    CHARACTER(SRCLEN3) CSRC             !  tmp source chars
    CHARACTER(PLTLEN3) FCID             !  tmp facility code
    CHARACTER(PLTLEN3) PLT              !  tmp plant code
    CHARACTER(SRCLEN3) SRCBUF           !  tmp source chars
    CHARACTER(SCCLEN3) TSCC             !  tmp SCC code

    CHARACTER(16), PARAMETER :: PROGNAME = 'RDREPIN'     ! program name

    !***********************************************************************
    !   begin body of subroutine RDREPIN

    !.......  If not using the ASCII elevated file
    IF( .NOT. AFLAG ) THEN
        !.......  Set local variables for determining input inventory variables
        LRDREGN = ( ANY( ALLRPT%BYGEO1 )   .OR.     &
                    ANY( ALLRPT%BYCNRY )   .OR.     &
                    ANY( ALLRPT%BYSTAT )   .OR.     &
                    ANY( ALLRPT%BYCNTY )   .OR.     &
                    ANY( ALLRPT%BYPLANT )  .OR.     &
                    ANY( ALLRPT%BYORIS )   .OR.     &
                    ANY( ALLRPT%BYBOILER ) .OR.     &
                    ANY( ALLRPT%BYSPC )    .OR.     &
                    ANY_CVAL( NREPORT, ALLRPT%REGNNAM ) .OR. YFLAG )

        !.......  Build array of inventory variable names based on report settings
        !.......  Region code
        IF( LRDREGN ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'CIFIP'
        END IF

        !.......  Road class code
        IF( ANY( ALLRPT%BYRCL ) ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'IRCLAS'
        END IF

        !.......  SIC code
        IF( ANY( ALLRPT%BYSIC ) ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'CISIC'
        END IF

        !.......  SCC code
        IF( LRDREGN .OR. ANY( ALLRPT%BYSCC ) ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'CSCC'
        END IF

        !.......  INTEGRATE flag
        IF( ANY( ALLRPT%BYINTGR ) ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'CINTGR'
        END IF

        !.......  MACT code
        IF( ANY( ALLRPT%BYMACT ) ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'CMACT'
        END IF

        !.......  NAICS code
        IF( ANY( ALLRPT%BYNAICS ) ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'CNAICS'
        END IF

        !.......  ORIS code
        IF( ANY( ALLRPT%BYORIS ) ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'CORIS'
        END IF

        !.......  Boiler code
        IF( ANY( ALLRPT%BYBOILER ) ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'CBLRID'
        END IF

        !.......  Unit ID code
        IF( ANY( ALLRPT%BYUNIT ) ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'CNEIUID'
        END IF

        !.......  Source type code
        IF( ANY( ALLRPT%BYSRCTYP ) ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'CSRCTYP'
        END IF

        !.......  Source description
        IF( ANY( ALLRPT%BYSRC ) .OR.&
            ANY( ALLRPT%BYPLANT ) ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'CSOURC'
        END IF

        !.......  Stack parameters
        IF( ANY( ALLRPT%STKPARM ) ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'STKHT'
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'STKDM'
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'STKTK'
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'STKVE'
        END IF

        !.......  Fugutive parameters
        IF( ANY( ALLRPT%FUGPARM ) ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'FUGHGT'
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'FUGWID'
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'FUGLEN'
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'FUGANG'
        END IF

        !.......  Point-source lat-lon
        IF( ANY( ALLRPT%LATLON ) ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'XLOCA'
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'YLOCA'
        END IF

        !.......  Plant name
        IF( ANY( ALLRPT%SRCNAM ) ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'CPDESC'
        END IF

        !.......  Emissions release point type
        IF( ANY( ALLRPT%BYERPTYP ) ) THEN
            NINVARR = NINVARR + 1
            IVARNAMS( NINVARR ) = 'CERPTYP'
        END IF

        !.......  Allocate memory for and read in required inventory characteristics
        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

    ELSE

        !.......  Read ASCII elevated file
        ALLOCATE( ATTRUNIT( 9 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ATTRUNIT', PROGNAME )
        ATTRUNIT = ''

        ATTRUNIT( 6 ) = 'm'
        ATTRUNIT( 7 ) = 'm'
        ATTRUNIT( 8 ) = 'deg K'
        ATTRUNIT( 9 ) = 'm/s'

        DO I = 1, 3
            ASCREC = ASCREC + 1
            READ( ADEV, '(A)' ) LINE
        END DO

        !......  Read in x, y origin and grid cell size
        READ( ADEV, '(F10.0,F10.0)' ) XO, YO
        ASCREC = ASCREC + 1
        READ( ADEV, '(F10.0,F10.0)' ) XC, YC
        ASCREC = ASCREC + 1
        READ( ADEV, '(I10,I10)' ) COLS, ROWS
        ASCREC = ASCREC + 1

        !......  Read in point source characteristics
        DO I = 1, 3
            ASCREC = ASCREC + 1
            READ( ADEV, '(A)' ) LINE
        END DO

        ALLOCATE( CIFIP( NSRC ),&
                 CSOURC( NSRC ),&
                  STKHT( NSRC ),&
                  STKDM( NSRC ),&
                  STKTK( NSRC ),&
                  STKVE( NSRC ),&
                   STKX( NSRC ),&
                   STKY( NSRC ),&
                  LPING( NSRC ),&
                 LMAJOR( NSRC ),&
                GROUPID( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CIFIP...GROUPID', PROGNAME )
        LPING  = .FALSE.
        LMAJOR = .TRUE.
        GROUPID= 0

        DO I = 1, NSRC

            CHARS = ''        ! array

            ASCREC = ASCREC + 1
            READ( ADEV, 93010 ) N, XL, YL, FCID,&
                        CHARS( 1 ), CIFIP( I )

            CALL BLDCSRC( CIFIP( I ), FCID, CHARS(1), CHARS(2),&
                          CHARS(3), CHARS(4), CHARS(5),&
                          POLBLNK3, CSOURC( I ) )

            READ( ADEV, 93020 ) STKHT( I ), STKDM( I ),&
                        STKTK( I ), STKVE( I )
            ASCREC = ASCREC + 1

            IF( STKDM( I ) .LT. 0. ) THEN
                LPING( I ) = .TRUE.
                STKDM( I ) = -STKDM( I )
            END IF

            STKVE( I ) = STKVE( I ) / 3600.

            STKX( I ) = INT( ( XL - XO ) / XC ) + 1
            STKY( I ) = INT( ( YL - YO ) / YC ) + 1

            IF( STKX( I ) .LT. 0 .OR. STKX( I ) .GT. COLS&
                .OR. STKY( I ) .LT. 0 .OR. STKY( I ) .GT. ROWS&
                ) THEN
                WRITE( MESG, 94010 ) 'ERROR: Source at ' //&
                       'line', IREC, 'is outside of the grid.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        END DO

        !......  Leave read pointer at beginning of emissions data
        DO I = 1, 11
            ASCREC = ASCREC + 1
            READ( ADEV, '(A)' ) LINE
        END DO

    END IF

    !.......  Create unique source characteristic lists
    CALL GENUSLST

    !.......  If needed, read in gridding matrix
    !.......  Initialize all to 1 for point sources
    IF( GFLAG ) THEN

        MESG = 'Reading gridding matrix...'
        CALL M3MSG2( MESG )

        CALL RDGMAT( GNAME, NGRID, NMATX, NMATX, NX, IX, CX )

        !.......  Initialize part of gridding matrix array for point sources
        IF( CATEGORY .EQ. 'POINT' ) THEN
            CX = 1       ! array
        END IF

    END IF

    !.......  If needed, read in gridding supplementation fike
    IF( GSFLAG ) THEN

        MESG = 'Reading supplemental gridding file...'
        CALL M3MSG2( MESG )

        ALLOCATE( SRGID( NSRC,2 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRGID', PROGNAME )
        SRGID = -9

        MESG = 'Supplemental gridding file'
        N = GETFLINE( GDEV, MESG )

        IREC = 0
        DO I = 1, N

            READ( GDEV, *, END=999, IOSTAT=IOS ) S, SRGID1, SRGID2
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'I/O error', IOS,&
                  'reading supplemental gridding file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

            SRGID( S,1 ) = SRGID1
            SRGID( S,2 ) = SRGID2

        END DO

    END IF

    !......  Allocate memory for projection factors
    !......  Ensure that NVPROJ = 1 if no controls are being run,
    !        because genrprt.f will still use array. First column
    !        will always be an array of ones.
    I = NVPROJ
    IF( PRRPTFLG ) I = 2 * I
    ALLOCATE( PRMAT( NSRC,1+I ), STAT=IOS )
    CALL CHECKMEM( IOS, 'PRMAT', PROGNAME )
    PRMAT = 1.

    !......  If needed, read in projection matrix
    IF( PRFLAG ) THEN

        MESG = 'Reading projection matrix...'
        CALL M3MSG2( MESG )

        !.......  Read in projection factors for each projection variable
        !.......  Note that openrepin.f contrains the no. of vars to 1,
        !         since that is how Cntlmat currently works.
        DO V = 1, NVPROJ
            IF( .NOT. READSET( PRNAME, PNAMPROJ( V ), ALLAYS3,&
                             ALLFILES, 0, 0, PRMAT( 1,1+V ) ) ) THEN

                MESG = 'ERROR: Could not read "' // TRIM( PNAMPROJ( V ) ) //&
                      '" from "' // TRIM( PRNAME ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END DO
    END IF

    !......  Read projection factors report...
    IF( PRRPTFLG ) THEN

        MESG = 'Reading projection report...'
        CALL M3MSG2( MESG )

        !......  Get the number of lines in the file
        MESG = 'projection report'
        NREPLIN = GETFLINE( RDEV( 1 ), MESG )

        !......  Loop through all lines in the file
        DO I = 1, NREPLIN

            !.......  Read line
            READ( RDEV( 1 ), '(A)', END = 9001 ) LINE

            !.......  Skip header lines (not very robust way)
            IF( I .LE. 6 ) CYCLE

            L = LEN_TRIM( LINE )
            NS = GETNLIST( L, LINE )

            !.......  Parse line into parts
            SEGMENT = ' '              ! array
            CALL PARSLINE( LINE, NS, SEGMENT )
            CFIP = SEGMENT( 1 )( 1:FIPLEN3 )

            !.......  Build source string
            SELECT CASE( CATEGORY )
              CASE ( 'AREA' )
                TSCC = SEGMENT( 2 )( 1:SCCLEN3 )
                CALL BLDCSRC( CFIP, TSCC, CHRBLNK3, CHRBLNK3,&
                              CHRBLNK3, CHRBLNK3, CHRBLNK3,&
                              POLBLNK3, CSRC )

              CASE ( 'MOBILE' )
                TSCC = SEGMENT( 2 )( 1:SCCLEN3 )
                CLNK = SEGMENT( 3 )( 1:LNKLEN3 )
                CALL BLDCSRC( CFIP, TSCC, CLNK, CHRBLNK3,&
                              CHRBLNK3, CHRBLNK3, CHRBLNK3,&
                              POLBLNK3, CSRC )

              CASE ( 'POINT' )
                PLT = SEGMENT( 2 )( 1:PLTLEN3 )
                DO J = 3, MAX( 8,NS-1 )
                    CHARS( J-2 ) = SEGMENT( J )( 1:CHRLEN3 )
                END DO

                IF( JSCC .GT. 0 ) THEN
                    TSCC = SEGMENT( JSCC )
                ELSE
                    TSCC = SEGMENT( NS-1 )
                END IF

                CALL BLDCSRC( CFIP, PLT, CHARS( 1 ), CHARS( 2 ),&
                              CHARS( 3 ), CHARS( 4 ), CHARS( 5 ),&
                              POLBLNK3, SRCBUF )
                CSRC = SRCBUF( 1:SRCLEN3 ) // TSCC

            END SELECT

            !.......  Search for source string in source list
            S = FINDC( CSRC, NSRC, CSOURC )

            !.......  If not found, give error
            IF( S .LE. 0 ) THEN

                !.......  Check if the first segment is an integer, and
                !         if not, then there is garbage at the end of
                !         the report (or an old report with multiple
                !         reports in one file).  If so, end read loop.
                IF( .NOT. CHKINT( SEGMENT( 1 ) ) ) THEN
                    EXIT
                END IF

                !......  Otherwise, there is an error because the report
                !        file is not for the inventory used in this run.
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Projection report ' //&
                       'entry at line', I, 'could not' //&
                       CRLF() // BLANK10 // 'be matched to the inventory.'
                CALL M3MESG( MESG )
                CYCLE
            END IF

            !.......  If found store factor
            MSGFLAG = .FALSE.
            J = NS - NVPROJ
            DO V = 1, NVPROJ
                J = J + 1

                !......  Check to ensure that segment is a real first.
                IF( .NOT. CHKREAL( SEGMENT( J ) ) ) THEN
                    EFLAG = .TRUE.

                    !......  If message has not been written for this line...
                    IF( .NOT. MSGFLAG ) THEN
                        WRITE( MESG,94010 ) 'ERROR: Bad format or value at line', IREC,     &
                                            'of projection report.'
                        CALL M3MSG2( MESG )
                    END IF
                    MSGFLAG = .TRUE.
                    CYCLE

                !......  If field is a real, store it
                ELSE
                    PRMAT( S,1+NVPROJ+V ) = STR2REAL( SEGMENT( J ) )
                END IF

            END DO

        END DO
    END IF         ! end if projection report or not

    !......  Allocate memory for multiplicative control factors
    !......  Ensure that NVCMULT = 1 if no controls are being run,
    !        because genrprt.f will still use array.  
    !        First column will always be an array of ones.
    ALLOCATE( ACUMATX( NSRC,1+NVCMULT ), STAT=IOS )
    CALL CHECKMEM( IOS, 'ACUMATX', PROGNAME )
    ACUMATX = 1.

    !.......  If needed, read in multiplicative control matrix
    IF( CUFLAG ) THEN

        MESG = 'Reading multiplicative control matrix...'
        CALL M3MSG2( MESG )

        DO V = 1, NVCMULT
            IF( .NOT. READSET( CUNAME, PNAMMULT( V ), ALLAYS3,&
                               ALLFILES, 0, 0, ACUMATX( 1,1+V ) )) THEN

                MESG = 'ERROR: Could not read "' // TRIM( PNAMMULT( V ) ) // &
                       '" from file "' // TRIM( CUNAME ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF
        END DO
    END IF

    !.......  If needed, read in speciation matrices
    !.......  NOTE that only the variables that are needed are read in
    IF( SLFLAG .OR. SSFLAG ) THEN

        MESG = 'Reading speciation matrices...'
        CALL M3MSG2( MESG )

        IF( SLNAME .NE. ' ' ) BNAME = SLNAME
        IF( SSNAME .NE. ' ' ) BNAME = SSNAME

        !.......  Get file header for variable names
        IF ( .NOT. DESCSET( BNAME, ALLFILES ) ) THEN

            MESG = 'Could not get description of file "' // TRIM( BNAME ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

        K = 0
        DO V = 1, NSVARS

            IF( SPCOUT( V ) ) THEN

                K = K + 1
                DBUF = VDESCSET( V )
                IF( SLFLAG ) THEN
                    CALL RDSMAT( SLNAME, DBUF, SLMAT( 1,K ) )
                END IF

                IF( SSFLAG ) THEN
                    CALL RDSMAT( SSNAME, DBUF, SSMAT( 1,K ) )
                END IF

            END IF

        END DO

    END IF

    !.......  If needed, read in speciation supplementation file
    IF( PSFLAG ) CALL RDSSUP( PDEV )

    !.......  If needed, read in temporal supplementation matrix
    IF( TSFLAG ) THEN

        MESG = 'Reading supplemental temporal file...'
        CALL M3MSG2( MESG )

        ALLOCATE( CMON( NSRC ),&
                  CWEK( NSRC ),&
                  CDOM( NSRC ),&
                  CMND( NSRC ),&
                  CTUE( NSRC ),&
                  CWED( NSRC ),&
                  CTHU( NSRC ),&
                  CFRI( NSRC ),&
                  CSAT( NSRC ),&
                  CSUN( NSRC ),&
                  CMET( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CMON...CMET', PROGNAME )
        CDOM = ' '
        CWEK = ' '
        CMON = ' '
        CMND = ' '
        CTUE = ' '
        CWED = ' '
        CTHU = ' '
        CFRI = ' '
        CSAT = ' '
        CSUN = ' '
        CMET = ' '
        CBUF = ' '

        MESG = 'Supplemental temporal file'
        N = GETFLINE( TDEV, MESG )

        !.......  Skip file header
        READ( TDEV, * )

        S = 0
        IREC = 1
        DO I = 2, N

            READ( TDEV, '(A)', END=1003, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )&
                  'I/O error', IOS,&
                  'reading supplemental temporal file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

            !.......  Parse line into parts
            SEGMENT = ' '              ! array
            CALL PARSLINE( LINE, 20, SEGMENT )

            TTYP = TRIM   ( SEGMENT( 1 ) )
            NV   = STR2INT( SEGMENT( 2 ) )
            TBUF = TRIM   ( SEGMENT( 4 ) )

            !.......  NV > 1 is not supported
            IF ( NV .NE. 1 ) THEN
                MESG = 'Number of variables in TSUP file is not supported.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF( .NOT. FIRSTIME ) THEN
                FIRSTTYP = TTYP
                FIRSTIME = .TRUE.
            END IF

            IF( FIRSTTYP == TTYP .OR. TTYP == 'MET' ) S = S + 1

            SELECT CASE( TTYP )
              CASE ( 'MTH' )
                CMON( S ) = TBUF
              CASE ( 'WEK' )
                CWEK( S ) = TBUF
              CASE ( 'DOM' )
                CDOM( S ) = TBUF
              CASE ( 'MON' )
                CMND( S ) = TBUF
              CASE ( 'TUE' )
                CTUE( S ) = TBUF
              CASE ( 'WED' )
                CWED( S ) = TBUF
              CASE ( 'THU' )
                CTHU( S ) = TBUF
              CASE ( 'FRI' )
                CFRI( S ) = TBUF
              CASE ( 'SAT' )
                CSAT( S ) = TBUF
              CASE ( 'SUN' )
                CSUN( S ) = TBUF
              CASE ( 'MET' )
                CMET( S ) = TBUF
            END SELECT

        END DO

        !.......  Check file format, assuming that pollutants weren't processed
        !         in > 1 groups.  (this routine doesn't handle grouped processing)
        IF ( S .NE. NSRC ) THEN
            MESG = 'INTERNAL ERROR: ' // CRL// 'TSUP file has '//&
                   'inconsistent number of lines with NSRC'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
        END IF

    END IF

    !.......  If needed, read in country, state, county file
    IF( YFLAG ) THEN
        IF( USEEXPGEO() ) THEN
            CALL RDGEOCODES( NINVIFIP, INVCFIP )
        ELSE
            CALL RDSTCY( YDEV, NINVIFIP, INVCFIP )
        END IF
    END IF

    !.......  If needed, read in elevated source indentification file
    IF( VFLAG ) THEN
        LTMP = ( .NOT. LFLAG )
        CALL RDPELV( EDEV, NSRC, LTMP, NMAJOR, NPING )
    END IF

    !.......  If needed, read in SCC descriptions file
    IF( NFLAG ) CALL RDSCCDSC( NDEV )

    !.......  If needed, read in SIC descriptions file
    IF( NIFLAG ) CALL RDSICDSC( NIDEV )

    !.......  If needed, read in GSPRO descriptions file
    IF( SDFLAG ) CALL RDSPROFDSC( NPDEV )

    !.......  If needed, read in MACT descriptions file
    IF( NMFLAG ) CALL RDMACTDSC( NMDEV )

    !.......  If needed, read in NAICS descriptions file
    IF( NNFLAG ) CALL RDNAICSDSC( NNDEV )

    !.......  If needed, read in ORIS descriptions file
    IF( NOFLAG ) CALL RDORSDSC( NODEV )

    !.......  If needed, read in layer fractions file to identify elevated sources
    IF( LFLAG ) THEN

        IF( .NOT. ALLOCATED( LMAJOR ) ) THEN
            ALLOCATE( LMAJOR( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LMAJOR', PROGNAME )
            LMAJOR = .FALSE.       ! array
        END IF

        IF( .NOT. ALLOCATED( LPING ) ) THEN
            ALLOCATE( LPING( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LPING', PROGNAME )
            LPING  = .FALSE.       ! array
        END IF

        IF( .NOT. ALLOCATED( GROUPID ) ) THEN
            ALLOCATE( GROUPID( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GROUPID', PROGNAME )
            GROUPID = 0       ! array
        END IF

        MESG = 'Reading layer fractions...'
        CALL M3MSG2( MESG )

        JDATE = SDATE
        JTIME = STIME
        DO T = 1, NSTEPS

            IF( READ3( LNAME, 'LFRAC', 1, JDATE, JTIME, LFRAC1L ) ) THEN
                DO S = 1, NSRC
                    IF( LFRAC1L( S ) .LT. 1. ) LMAJOR( S ) = .TRUE.
                END DO

            ELSE      !  Read failed
                MESG = 'Could not read "LFRAC" from '// LNAME
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

            END IF

            CALL NEXTIME( JDATE, JTIME, TSTEP )

        END DO

    END IF

    !.......  Reformat source characteristics and set widths.  Do this once
    !         for the entire run of the program, so that it doesn't have to be
    !         done for each report (it is slow)

    IF( ANY( ALLRPT%BYSRC  ) .OR.&
        ANY( ALLRPT%BYPLANT ) ) THEN

        !.......  Determine width of source chararactistic columns over the
        !         whole inventory
        SWIDTH = 0           ! initialize array
        DO S = 1, NSRC

            K = 0
            DO J = MINC, NCHARS

                K  = K + 1
                L1 = SC_BEGP( J )
                L2 = SC_ENDP( J )
                BUFFER = ADJUSTL( CSOURC( S )( L1:L2 ) )
                SWIDTH( K ) = MAX( SWIDTH( K ), LEN_TRIM( BUFFER ) )

            END DO

        END DO

        !.......  Reset CSOURC based on these widths
        !.......  Also remove SCC from the source characteristics
        DO S = 1, NSRC

            CSRC = CSOURC( S )

            L = SC_ENDP( MINC-1 )
            K  = 0
            DO J = MINC, NCHARS

                K = K + 1
                IF( J .NE. JSCC ) THEN
                    L1 = SC_BEGP( J )
                    L2 = SC_ENDP( J )
                    CSOURC( S ) = CSOURC( S )( 1:L ) //&
                                  ADJUSTL( CSRC( L1:L2 ) )
                    L = L + SWIDTH( K )
                END IF

            END DO

        END DO

        !.......  Allocate and initialize arrays for storing new field lengths
        ALLOCATE( LOC_BEGP( NCHARS ),&
                  LOC_ENDP( NCHARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LOC_BEGP,LOC_ENDP', PROGNAME )
        LOC_BEGP = SC_BEGP           ! array
        LOC_ENDP = SC_ENDP           ! array

        !.......  Set local start and end fields based on new widths
        K = 0
        DO J = MINC, NCHARS
            K = K + 1
            IF( J .EQ. JSCC ) CYCLE
            LOC_BEGP( J ) = LOC_ENDP( J-1 ) + 1
            LOC_ENDP( J ) = LOC_BEGP( J ) + SWIDTH( K ) - 1
        END DO

        IF( JSCC .GT. 0 ) NCHARS = NCHARS - 1

    END IF

    !.......  Exit if any errors encountered
    IF( EFLAG ) THEN
        MESG = 'Problem(s) reading input files.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    RETURN

999 MESG = 'Unexpected end of file reached while reading supplementary gridding file.'
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

1001 MESG = 'Unexpected end of file reached while reading supplementary speciation file.'
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

1003 MESG = 'Unexpected end of file reached while reading supplementary temporal file.'
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

9001 MESG = 'Unexpected end of file reached while reading projection report file.'
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Formatted file I/O formats...... 93xxx

93010 FORMAT( I10, 10X, 2F10.0, 2A10, I10.5 )

93020 FORMAT( F10.1, F10.2, F10.1, F10.0 )

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I10, :, 1X ) )

    !******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------

    !.......  This internal function scans a character array for any
    !         non-blank values, and if it finds one, returns true.
    LOGICAL FUNCTION ANY_CVAL( NDIM, CHARARR )

        !.......  Subprogram arguments
        INTEGER     , INTENT (IN) :: NDIM
        CHARACTER(*), INTENT (IN) :: CHARARR( NDIM )

        !.......  Local variables
        INTEGER   I

        !----------------------------------------------------------------------

        ANY_CVAL = .FALSE.

        DO I = 1, NDIM

            IF( CHARARR( I ) .NE. ' ' ) THEN
                ANY_CVAL = .TRUE.
                RETURN
            END IF

        END DO

        RETURN

    END FUNCTION ANY_CVAL

END SUBROUTINE RDREPIN

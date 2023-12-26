
PROGRAM MRGGRID

    !***********************************************************************
    !  program body starts at line
    !
    !  DESCRIPTION:
    !    Program MRGGRID reads 2-D and 3-D I/O API files and merges them
    !    into a single 2-D or 3-D file (depending on the inputs)
    !    The time period merged is adjusted based on the latest
    !    starting file and earliest ending file, unless MRG_DIFF_DAY is
    !    set in which case the time period is based on the standard
    !    environment variables. All variables are merged, even if different
    !    variables are in each file.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Original by M. Houyoux 4/98
    !       Modified by B.H. Baek  4/08
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***********************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !       System
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

    !.....   MODULES for public variables
    !.....  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: NGRID, NCOLS, NROWS, NLAYS, VGLVS, VGTYP, VGTOP

    !.....  This module is required for the FileSetAPI
    USE MODFILESET, ONLY : FILE_INFO, RNAMES, NVARSET, VNAMESET, VUNITSET, VDESCSET

    IMPLICIT NONE

    !.....   INCLUDES:
    INCLUDE 'CONST3.EXT'        !  physical constants from I/O API

    INCLUDE 'EMCNST3.h90'
    INCLUDE 'SETDECL.h90'       !  FileSetAPI variables and functions

    !.....   EXTERNAL FUNCTIONS
    LOGICAL, EXTERNAL :: BLKORCMT
    INTEGER, EXTERNAL :: GETFLINE
    LOGICAL, EXTERNAL :: CHKREAL

    !.....  LOCAL PARAMETERS and their descriptions:

    CHARACTER(16), PARAMETER :: PROGNAME = 'MRGGRID'     ! program name
    CHARACTER(50), PARAMETER ::  CVSW    = '$Name SMOKEv5.0_Jun2023$'     ! CVS release tag

    !.....   LOCAL VARIABLES and their descriptions:
    CHARACTER(16) :: SEGMENT( 3 )

    !.....   Emissions arrays
    REAL, ALLOCATABLE :: E2D ( : )        ! 2-d emissions
    REAL, ALLOCATABLE :: E3D( :,: )       ! 3-d emissions
    REAL, ALLOCATABLE :: EOUT( :,: )      ! output emissions
    REAL, ALLOCATABLE :: BEFORE_ADJ( : )      ! emissions before factors applied
    REAL, ALLOCATABLE :: AFTER_ADJ ( : )      ! emissions after factors applied
    REAL, ALLOCATABLE :: BEFORE_SPC( : )      ! emissions before factors applied
    REAL, ALLOCATABLE :: AFTER_SPC ( : )      ! emissions after factors applied

    !.....   Input file descriptors
    INTEGER,       ALLOCATABLE :: DURATA( : )     ! no. time steps
    INTEGER,       ALLOCATABLE :: NCOLSA( : )     ! no. columns
    INTEGER,       ALLOCATABLE :: NROWSA( : )     ! no. rows
    INTEGER,       ALLOCATABLE :: NVARSA( : )     ! no. variables
    INTEGER,       ALLOCATABLE :: SDATEA( : )     ! start date
    INTEGER,       ALLOCATABLE :: STIMEA( : )     ! start time
    INTEGER,       ALLOCATABLE :: NLAYSA( : )     ! number of layers in the file
    INTEGER,       ALLOCATABLE :: NFILES( : )     ! number of files in each fileset
    INTEGER,       ALLOCATABLE :: ID1(:), ID2(:)     ! temperature bin indexes
    LOGICAL,       ALLOCATABLE :: USEFIRST(:)     ! true: use first time step of file
    LOGICAL,       ALLOCATABLE :: LVOUTA( :,: )     ! iff out var in input file
    REAL,          ALLOCATABLE :: ADJ_FACTOR( : )     ! adjustment factors
    REAL,          ALLOCATABLE :: TA( : )         ! ambient temperatures
    REAL,          ALLOCATABLE :: TABASE( : )     ! base ambient temperatures
    REAL,          ALLOCATABLE :: QV( : )         ! specific humidity
    REAL,          ALLOCATABLE :: RA( : )         ! aerodynamic resistance
    REAL,          ALLOCATABLE :: RC( : )         ! precipitatoin (cm)
    REAL,          ALLOCATABLE :: RABASE( : )     ! base aerodynamic resistance
    REAL,          ALLOCATABLE :: TEMPS( : )      ! temperature bins for layered emissions table
    REAL,          ALLOCATABLE :: METRWC( : )     ! met adjustment factor for RWC sector
    REAL,          ALLOCATABLE :: METNH3( : )     ! met adjustment factor for NH3 sector
    REAL,          ALLOCATABLE :: NOXGAS( : )     ! Humidity NOx adjustment factor for gasoline
    REAL,          ALLOCATABLE :: NOXDIS( : )     ! Humidity NOx adjustment factor for diesel
    CHARACTER(16), ALLOCATABLE :: VNAMEA( :,: )     ! variable names
    CHARACTER(16), ALLOCATABLE :: VUNITA( :,: )     ! variable units
    CHARACTER(80), ALLOCATABLE :: VDESCA( :,: )     ! var descrip
    CHARACTER(16), ALLOCATABLE :: IONAME( : )     ! IOAPI 16chr 2-d input file names
    CHARACTER(16), ALLOCATABLE :: METADJ( : )     ! Met adjustment Y|N
    CHARACTER(16), ALLOCATABLE :: NOXADJ( : )     ! NOx adjustment for mobile source Y|N
    CHARACTER(32), ALLOCATABLE :: FNAME ( : )     ! 2-d input file names
    CHARACTER(32), ALLOCATABLE :: ADJ_LFN( : )    ! Species name
    CHARACTER(16), ALLOCATABLE :: ADJ_SPC( : )    ! logicalFileName
    CHARACTER(49), ALLOCATABLE :: ADJ_LFNSPC( : )     ! concatenated {logicalFileName}_{Species}
    CHARACTER(32), ALLOCATABLE :: TAG_LFN( : )    ! Tagging species name
    CHARACTER(16), ALLOCATABLE :: TAG_SPC( :    )     ! Tagging logicalFileName
    CHARACTER(49), ALLOCATABLE :: TAG_LFNSPC( : )     ! Tagging concatenated {logicalFileName}_{Species}
    CHARACTER(16), ALLOCATABLE :: TAG_APPEND( : )      ! Tagging index
    CHARACTER(66), ALLOCATABLE :: TAG_LFNSPCTAG( : )      ! Tagging index(logicalfile+species+tag)
    CHARACTER(16)                 VNAMEP( MXVARS3 )     ! pt variable names
    CHARACTER(16)                 VUNITP( MXVARS3 )     ! pt variable units
    CHARACTER(80)                 VDESCP( MXVARS3 )     ! pt var descrip

    !.....   Intermediate output variable arrays
    INTEGER       INDXN ( MXVARS3 )     ! sorting index for OUTIDX
    INTEGER       OUTIDX( MXVARS3 )     ! index to master model species list

    CHARACTER(16) OUTNAM( MXVARS3 )     ! unsorted output variable names
    CHARACTER(16) VUNITU( MXVARS3 )     ! unsorted output variable units
    CHARACTER(80) VDESCU( MXVARS3 )     ! unsorted output variable descriptions

    LOGICAL       LVOUTP( MXVARS3 )     ! iff output var exists in point input

    !.....   Logical names and unit numbers

    INTEGER       ADEV            ! unit for logical names list for SEG
    INTEGER       IDEV            ! unit for logical names list for 2d files
    INTEGER       LDEV            ! unit for log file
    INTEGER       RDEV            ! unit for merge report file
    INTEGER       ODEV            ! unit for QA report file
    INTEGER       SDEV            ! unit for overall QA report file
    INTEGER       TDEV            ! unit for taggin input file
    INTEGER       GDEV            ! unit for taggin species QA file
    CHARACTER(16) ONAME           ! Merged output file name
    CHARACTER(16) PNAME           ! Point source input file name

    !.....   Other local variables
    INTEGER       C, DD, F, I, J, K, L, L1, L2, N, NL, V, S, T     ! pointers and counters

    INTEGER       ADJ                        ! tmp adjustment factor main index
    INTEGER       ADJ1                       ! tmp adjustment factor index 1
    INTEGER       ADJ2                       ! tmp adjustment factor index 2
    INTEGER       TAG                        ! tmp tagging species index
    INTEGER       DUMMY                      ! dummy value for use with I/O API functions
    INTEGER       EDATE                      ! ending julian date
    INTEGER       ETIME                      ! ending time HHMMSS
    INTEGER    :: G_SDATE = 0                ! start date from environment
    INTEGER    :: G_STIME = 0                ! start time from environment
    INTEGER    :: G_NSTEPS = 1               ! number of time steps from environment
    INTEGER    :: G_TSTEP = 0                ! time step from environment
    INTEGER       ICNTFIL                    ! tmp count of fileset file count
    INTEGER       IOS                        ! i/o status
    INTEGER       IREC                       ! line number count
    INTEGER       JDATE                      ! iterative julian date
    INTEGER       JTIME                      ! iterative time HHMMSS
    INTEGER       LB                         ! leading blanks counter
    INTEGER       LE                         ! location of end of string
    INTEGER       MXNF                       ! tmp no. of 2-d input files
    INTEGER       MXNFIL                     ! max no. of 2-d input files
    INTEGER       MXNFAC                     ! max no. of adjustment factors
    INTEGER       MXNTAG                     ! max no. of tagging species
    INTEGER    :: NADJ = 0                   ! no. of adjustment factors
    INTEGER    :: NTAG = 0                   ! no. of tagging species
    INTEGER       NTEMPS,IDX1,IDX2,LNOX      ! no. of 3-d layered emissions table with idxes
    INTEGER       NFILE                      ! no. of 2-d input files
    INTEGER       NSTEPS                     ! no. of output time steps
    INTEGER       NVOUT                      ! no. of output variables
    INTEGER       RDATE                      ! reference date
    INTEGER       SAVLAYS                    ! number of layers
    INTEGER       SDATE                      ! starting julian date
    INTEGER       SECS                       ! tmp seconds
    INTEGER       SECSMAX                    ! seconds maximum
    INTEGER       SECSMIN                    ! seconds minimum
    INTEGER       STIME                      ! starting time HHMMSS
    INTEGER       STEPS                      ! tmp number of steps
    INTEGER       TIMET                      ! tmp time from seconds
    INTEGER       TSTEP                      ! time step
    INTEGER       VLB                        ! VGLVS3D lower bound

    REAL       :: FACS = 1.0                 ! adjustment factor
    REAL          SUME2D                     ! SUM( E2D(1:NGRID) )
    REAL          RATIO,TEMPVAL,TMPVAL       ! ratio and temp val
    REAL          X, Y, CX, CY               ! temp values
    REAL          RWC_TEMP                   ! RWC related values
    REAL          BASE                       ! AGNH3 base temp ana wsp
    REAL          QVMOL, A, B                ! NOx humidity correctino eqs vars

    CHARACTER(16)  FDESC                     ! tmp file description
    CHARACTER(16)  MRGFDESC                  ! name for file description EV
    CHARACTER(128) METADESC                  ! output meta description from MRGFDESC
    CHARACTER(16)  METNAM, NOXNAM            ! tmp met and nox adj names
    CHARACTER(16)  IO_NAM                    ! tmp 16 chr logical file name
    CHARACTER(32)  NAM                       ! tmp logical file name
    CHARACTER(32)  LNAM                      ! tmp previous file name
    CHARACTER(16)  VNM                       ! tmp variable name
    CHARACTER(16)  TVNM                      ! tmp2 variable name
    CHARACTER(16)  TVARNAME, RANAME, QVNAME, RCNAME      ! temp variable name
    CHARACTER(49)  LFNSPC                    ! tmp spec and file name
    CHARACTER(66)  LFNSPCTAG                 ! tmp speC, file, and tag name
    CHARACTER(256) LINE                      ! input buffer
    CHARACTER(256) NAMBUF                    ! tmp buffer for logical file name
    CHARACTER(256) MESG                      ! message field
    CHARACTER(80)  NAME1                     ! tmp file name component
    CHARACTER(15)  RPTCOL                    ! single column in report line
    CHARACTER(10)  EFMT                      ! output emissions foamat
    CHARACTER(100) REPFMT                    ! output emissions foamat
    CHARACTER(300) REPFILE                   ! name of report file
    CHARACTER(300) RPTLINE                   ! line of report file
    CHARACTER(16)  SPCTMP                    ! tmp species name
    CHARACTER(32)  LFNTMP                    ! tmp file name

    LOGICAL    :: EFLAG   = .FALSE.       ! error flag
    LOGICAL    :: ETFLAG  = .FALSE.       ! true: 3-d emissions table
    LOGICAL    :: FIRST3D = .TRUE.    ! true: first 3-d file not yet input
    LOGICAL    :: LFLAG   = .FALSE.       ! true  if 3-d file input
    LOGICAL    :: TFLAG   = .FALSE.       ! true: grid didn't match
    LOGICAL    :: MRGDIFF = .FALSE.       ! true: merge files from different days
    LOGICAL    :: RWCFLAG = .FALSE.       ! true: adjust RWC`emissions with met
    LOGICAL    :: NH3FLAG = .FALSE.       ! true: adjust NH3 `emissions with met
    LOGICAL    :: METFLAG = .FALSE.       ! true: adjust emissions with met
    LOGICAL    :: MOBFLAG = .FALSE.       ! true: adjust mobile emissions with met
    LOGICAL    :: NOXFLAG = .FALSE.       ! true: adjust emissions with met
    LOGICAL    :: NOXADJEQS = .FALSE.     ! true: apply the latest MOVES3 NOx humidity correction equations

    !***********************************************************************
    !   begin body of program MRGGRID

    LDEV = INIT3()

    !.....  Write out copyright, version, web address, header info, and prompt
    !       to continue running the program.
    CALL INITEM( LDEV, CVSW, PROGNAME )

    !.....  Read names of input files and open files
    MESG = 'Enter logical name for 2-D AND 3-D GRIDDED INPUTS list'

    IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE.,    &
                        'FILELIST', PROGNAME   )

    !.....  Get environment variables
    MESG = 'Merge files from different days into single file'
    MRGDIFF = ENVYN( 'MRG_DIFF_DAYS', MESG, .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_DIFF_DAYS"', 2 )
    END IF

    IF( MRGDIFF ) THEN
    !.....  Get date and time settings from environment
        CALL GETM3EPI( -1, G_SDATE, G_STIME, G_TSTEP, G_NSTEPS )
    END IF

    !.....  Determine maximum number of input files in file
    MXNFIL = GETFLINE( IDEV, 'List of files to merge' )

    !.....  Write message out about MXVARS3
    WRITE( MESG,94010 ) 'Mrggrid compiled with I/O API MXVARS3 =',    &
                        MXVARS3
    CALL M3MSG2( MESG )

    !.....  Allocate memory for arrays that just depend on the maximum number
    !       of files
    ALLOCATE( NFILES( MXNFIL ),    &
              DURATA( MXNFIL ),    &
              NCOLSA( MXNFIL ),    &
              NROWSA( MXNFIL ),    &
              NLAYSA( MXNFIL ),    &
              NVARSA( MXNFIL ),    &
              SDATEA( MXNFIL ),    &
              STIMEA( MXNFIL ),    &
               FNAME( MXNFIL ),    &
              IONAME( MXNFIL ),    &
              METADJ( MXNFIL ),    &
              NOXADJ( MXNFIL ),    &
            USEFIRST( MXNFIL ),    &
              LVOUTA( MXVARS3,MXNFIL ),    &
              VNAMEA( MXVARS3,MXNFIL ),    &
              VUNITA( MXVARS3,MXNFIL ),    &
              VDESCA( MXVARS3,MXNFIL ),    &
              VGLVS( 0:MXLAYS3 ), STAT=IOS )
    CALL CHECKMEM( IOS, 'NFILES...VGLVS', PROGNAME )
    VGLVS = 0.

    !.....  Loop through input files and open them
    F = 0
    IREC = 0
    DO

        !.....  Read file names - exit if read is at end of file
        READ( IDEV, 93000, END=27, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF ( IOS .NE. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'I/O error', IOS, 'reading file list at line', IREC
            CALL M3MESG( MESG )
            CYCLE
        END IF

        !.....  Skip blank and comment lines
        IF ( BLKORCMT( LINE ) ) CYCLE

        SEGMENT = ''
        CALL PARSLINE( LINE, 3, SEGMENT )

        F = F + 1

        IF( F .LE. MXNFIL ) THEN

            LB = LBLANK(   SEGMENT( 1 ) )
            LE = LEN_TRIM( SEGMENT( 1 ) )
            FNAME( F ) = SEGMENT( 1 )( LB+1:LE )

            !.....  Re-set logical file name
            IF( LE > 16 ) THEN
                CALL GETENV( FNAME( F ), NAMBUF )

                WRITE( IO_NAM,93030 ) FNAME( F )( 1:13 ) // '_', F

                IF( .NOT. SETENVVAR( IO_NAM, NAMBUF ) ) THEN
                    MESG = 'Could not set logical name for ' // 'file ' // NAMBUF
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                IONAME( F ) = IO_NAM

            ELSE

                IO_NAM = TRIM( FNAME( F ) )
                IONAME( F ) = IO_NAM

            ENDIF

            METADJ( F ) = SEGMENT( 2 )
            NOXADJ( F ) = SEGMENT( 3 )
            L1 = LEN_TRIM( METADJ( F ) )
            L2 = LEN_TRIM( NOXADJ( F ) )

            IF( L1 > 0 ) METFLAG = .TRUE.
            IF( METADJ( F ) == 'RWC'    ) RWCFLAG = .TRUE.
            IF( METADJ( F ) == 'AGNH3'  ) NH3FLAG = .TRUE.
            IF( METADJ( F ) == 'METEMIS' ) MOBFLAG = .TRUE.
            IF( METFLAG .AND. L2 > 0 ) NOXFLAG = .TRUE.

            IF ( .NOT. OPENSET( IONAME(F), FSREAD3, PROGNAME ) ) THEN

                MESG = 'Could not open file "' // TRIM( FNAME(F) ) // '".'
                CALL M3MSG2( MESG )
                MESG = 'Ending program "MRGGRID".'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF      !  if open3() failed

            !.....  Store whether it's a fileset file or not
            I = INDEX1( FNAME(F), MXFILE3, RNAMES )
            NFILES( F ) = SIZE( FILE_INFO( I )%LNAMES )

        END IF

    END DO
27  CONTINUE

    NFILE = F

    IF( NFILE .GT. MXNFIL ) THEN
        WRITE( MESG,94010 )    &
          'INTERNAL ERROR: Dimension mismatch.  Input file count:',    &
          NFILE, 'program allows (MXNFIL):', MXNFIL
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    ELSEIF( NFILE .EQ. 0 ) THEN
        MESG = 'No input files in list    !'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    ENDIF

    !.....  Get met adjustment variable setting from FILELIST input file
    !.....  Adjust mobile emissiosn with meteorology
    IF( METFLAG ) THEN
        MESG = 'CRITICAL: Meteorological adjustment to gridded emissions ' //    &
               'has been enabled    !'
        CALL M3MESG( MESG )

        MESG = 'Specifies ambient temperature variable name for ' //    &
               'meteorological adjustment'
        CALL ENVSTR( 'TEMP_VAR', MESG, 'TEMP2',TVARNAME, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "OUTZONE"', 2 )
        END IF
        CALL UPCASE( TVARNAME )

        MESG = 'Specifies aerodynamic resistance variable name for ' //    &
               'meteorological adjustmen for NH3 emissions'
        CALL ENVSTR( 'RA_VAR', MESG, 'RADYNI', RANAME, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "TEMP_VAR"', 2 )
        END IF
        CALL UPCASE( RANAME )

        MESG = 'Specifies precipitation variable name for ' //    &
               'meteorological adjustmen for NH3 emissions'
        CALL ENVSTR( 'RC_VAR', MESG, 'RC', RCNAME, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "RC_VAR"', 2 )
        END IF
        CALL UPCASE( RCNAME )

        IF( NOXFLAG ) THEN
            MESG = 'CRITICAL: Humidity correction for NOx emissions ' //    &
               'from mobile sources has been enabled    !'
            CALL M3MESG( MESG )

            MESG = 'Specifies specific humidity variable name for ' //    &
                   'humidity adjustmen for mobile NOx emissions'
            CALL ENVSTR( 'QV_VAR', MESG, 'QV',QVNAME, IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "QV_VAR"', 2 )
            END IF
            CALL UPCASE( QVNAME )

            NOXADJEQS = ENVYN( 'USE_MOVES3_NOX_ADJ_EQS',          &
                       'Use the MOVES3 NOx humidity correction equations',    &
                       .FALSE., IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "USE_MOVES3_NOX_ADJ_EQS"', 2 )
            END IF
        END IF

        !.....  RWC source meteorology adjustment
        MESG = 'Enter RWC Ambient Threshold Temperature in unit of Fahrenheit (F)'
        RWC_TEMP = ENVREAL( 'RWC_TEMP_THRESHOLD', MESG, 45.0, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "RWC_TEMP_THRESHOLD"', 2 )
        END IF

    END IF

    !.....  Get environment variable settings for adjustment factor input file
    CALL ENVSTR( 'ADJ_FACS', MESG, ' ', NAME1 , IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "ADJ_FACS"', 2 )
    END IF

    !.....  Determine maximum number of input files in file
    IF( IOS < 0 ) THEN     !  failure to open
        ADEV = IOS
        MXNFAC = 1
        MESG = 'NOTE : No adjustment factors were applied because'//    &
               ' there is no ADJ_FACS environment variable defined'
        CALL M3MSG2( MESG )

    ELSE
        MESG = 'Enter logical name for a list of adjustment factors'
        ADEV = PROMPTFFILE( MESG,.TRUE.,.TRUE.,'ADJ_FACS',PROGNAME )
        MXNFAC = GETFLINE( ADEV, 'List of adjustment factos' )

        !....  Write summary of sector specific factor adjustment output
        ODEV = PROMPTFFILE(    &
           'Enter logical name for the MRGGRID QA REPORT file',    &
           .FALSE., .TRUE., 'REPMERGE_ADJ', PROGNAME )

        MXNF = 0
        DO         ! head of report file
            READ( ODEV, 93000, END=222 ) LINE
            MXNF = MXNF + 1
        ENDDO
222     CONTINUE

        !.....  Write header line to report
        IF( MXNF == 0 ) THEN
            WRITE( ODEV,93000 ) '#MRGGRID logical file QA Report'
            WRITE( ODEV,93000 ) '#COLUMN_TYPES=Int(4)|Varchar(32)|' //    &
                         'Varchar(16)|Real(8)|Real(8)|Real(8)|Real(8)'
            WRITE( ODEV,93000 ) 'DATE,FileName,Species,Factor,'//    &
                              'Before,After,Ratio'
        END IF

        !....  Write summary of overall factor adjustment output by species
        SDEV = PROMPTFFILE(    &
           'Enter logical name for the MRGGRID Overall REPORT file',    &
           .FALSE., .TRUE., 'REPMERGE_SUM', PROGNAME )

        MXNF = 0
        DO             ! head of report file
            READ( SDEV, 93000, END=333 ) LINE
            MXNF = MXNF + 1
        ENDDO
333     CONTINUE

        !.....  Write header line to report
        IF( MXNF == 0 ) THEN
            WRITE( SDEV,93000 ) '#MRGGRID Overall Summary by Species'
            WRITE( SDEV,93000 ) '#COLUMN_TYPES=Int(4)|Varchar(16)|' //    &
                              'Real(8)|Real(8)|Real(8)'
            WRITE( SDEV,93000 ) 'DATE,Species,Before,After,Ratio'
        END IF

    END IF

    !.....  Allocate memory for arrays that just depend on the maximum number
    !       of adjustment factors in ADJ_FACS input file.
    ALLOCATE( ADJ_LFN( MXNFAC ),    &
              ADJ_SPC( MXNFAC ),    &
           ADJ_LFNSPC( MXNFAC ),    &
           ADJ_FACTOR( MXNFAC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'ADJ_FACTOR', PROGNAME )

    ADJ_SPC = ' '
    ADJ_LFN = ' '
    ADJ_LFNSPC = ' '
    ADJ_FACTOR = 0.0

    !.....  Define a number of adjustment factors
    IF( ADEV < 0 ) THEN
        NADJ = 1
    ELSE
    !.....  Store a list of adjustment factors
        CALL READ_ADJ_FACS( NADJ )
    END IF

    !.....  Duplicate Check of ADJ_FACS file
    DO  F = 1, NADJ
        LFNSPC = ADJ_LFNSPC( F )
        DD = 0
        DO I = 1, NADJ
            IF( LFNSPC == ADJ_LFNSPC( I ) ) DD = DD + 1
        END DO

        IF( DD > 1 ) THEN
            MESG = 'ERROR: Duplicate entries of '// TRIM(ADJ_SPC(F))    &
                 // ' species from the ' // TRIM( ADJ_LFN(F) ) //    &
                ' file in the ADJ_FACS file.' // LFNSPC
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.
        END IF
    ENDDO

    !.....  Give error message and end program unsuccessfully
    IF( EFLAG ) THEN
        MESG = 'ERROR: Duplicate entries in the ADJ_FACS file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.....  Allocate arrays that will store sector-specific daily/gridded total emissinos
    ALLOCATE( BEFORE_ADJ( NADJ ),    &
               AFTER_ADJ( NADJ ), STAT=IOS )
    CALL CHECKMEM( IOS, 'AFTER_ADJ', PROGNAME )
    BEFORE_ADJ = 0.0
    AFTER_ADJ  = 0.0

    !.....  Get environment variable settings for tagging species input file
    CALL ENVSTR( 'TAG_SPECIES', MESG, ' ', NAME1 , IOS )
    !.....  Determine maximum number of input files in file
    IF( IOS < 0 ) THEN     !  failure to open
        TDEV = IOS
        MESG = 'NOTE : No tagging species were available because'//    &
               ' there is no TAG_SPECIES environment variable defined'
        CALL M3MSG2( MESG )

    ELSE

        MESG = 'NOTE : Tagging species based upon TAG_SPECIES file'
        CALL M3MSG2( MESG )

        !.....  Store a list of tagging species
        CALL READ_TAG_SPECIES( NTAG )

        !....  Write summary of sector specific factor adjustment output
        GDEV = PROMPTFFILE(    &
           'Enter logical name for the MRGGRID Tagging REPORT file',    &
           .FALSE., .TRUE., 'REPMERGE_TAG', PROGNAME )

        !.....  Write header line to report
        WRITE( GDEV,93000 ) '#MRGGRID Tagging species Report'
        WRITE( GDEV,93000 ) '#COLUMN_TYPES=Varchar(32)|' //    &
                            'Varchar(16)|Varchar(16)'
        WRITE( GDEV,93000 ) 'FileName,OriginalSpecies,TaggedSpecies'

    END IF

    !.....  Determine I/O API layer storage lower bound
    VLB = LBOUND( VGLVS3D,1 )

    !.....  Get file descriptions and store for all input files
    !.....  Loop through 2D input files
    NLAYS = 1
    DO F = 1, NFILE

        NAM    = FNAME( F )
        IO_NAM = IONAME( F )       ! retrieve 16 char ioapi local file name
        METNAM = METADJ( F )       ! 'METEMIS', 'AGNH3', or 'RWC'

        ICNTFIL = ALLFILES
        IF( NFILES( F ) .EQ. 1 ) ICNTFIL = 1       ! send ALLFILES if more than one file, send 1 otherwise
        IF ( .NOT. DESCSET( IO_NAM, ICNTFIL ) ) THEN
            MESG = 'Could not get description of file "'  //    &
                    NAM( 1:LEN_TRIM( NAM ) ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSE
        !.....  Reset NLAYS3D and VGLV3D when processing layered ETABLES from Movesmrg
            IF( METNAM == 'METEMIS' ) THEN           ! when RatePer (RP) modes layered mobile emission table are processed
                NTEMPS = NLAYS3D
                IF( .NOT. ALLOCATED( TEMPS ) ) THEN
                    ALLOCATE( TEMPS( NTEMPS ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'TEMPS', PROGNAME )
                END IF
                TEMPS( 0:NTEMPS ) = VGLVS3D( 0:NLAYS3D )
                NLAYS3D = 1
                VGLVS3D = 0.0
                ETFLAG = .TRUE.
            END IF

            NROWSA( F ) = NROWS3D
            NCOLSA( F ) = NCOLS3D
            NLAYSA( F ) = NLAYS3D
            NVARSA( F ) = NVARSET
            SDATEA( F ) = SDATE3D
            STIMEA( F ) = STIME3D
            DURATA( F ) = MXREC3D

            IF( F == 1 ) TSTEP = TSTEP3D

            DO V = 1, NVARSET
                VNAMEA( V,F ) = VNAMESET( V )
                VUNITA( V,F ) = VUNITSET( V )
                VDESCA( V,F ) = VDESCSET( V )

                IF( TDEV > 0 ) THEN
                !.....  Set tmp variables
                    VNM = VNAMESET( V )

                !....  Search tagged species for the current file
                    LFNSPC = TRIM( NAM ) // '~' // TRIM( VNM )

                    TAG = INDEX1( LFNSPC, NTAG, TAG_LFNSPC )
                    ADJ1= INDEX1( LFNSPC, NADJ, ADJ_LFNSPC )

                    !.....  Assign tagged species for the current species
                    IF( TAG > 0 ) THEN
                        TVNM = TRIM( VNM ) // '_' // TRIM( TAG_APPEND( TAG ) )

                        MESG = 'NOTE : Appending a tag (' //    &
                            TRIM( TAG_APPEND(TAG) ) // ') to the '    &
                            //'species ' //TRIM( VNM )// ' from the '    &
                            //TRIM( NAM )// ' file'
                        CALL M3MSG2( MESG )

                        !.....  Replace a species name with a tagged one
                        VNAMEA( V,F ) = TVNM

                        !.....  Write the changes to tagging summary report
                        WRITE( GDEV, 92000 ) NAM, VNM, TVNM

                        !.....  Error and Warning messages for the tagged species
                        !       before you apply the adjustment factors if necessary

                        LFNSPC = TRIM( NAM ) // '~' // TRIM( TVNM )

                        ADJ2 = INDEX1( LFNSPC, NADJ, ADJ_LFNSPC )

                        IF( ADJ1 > 0 .AND. ADJ2 < 1 ) THEN
                            MESG ='WARNING : Adjustment factor ' //    &
                                ' for the species ' // TRIM( VNM )    &
                                // ' from file ' // TRIM( NAM ) //    &
                               ' will be skipped due to ' //    &
                                'the change of species name to ' //    &
                                TRIM( TVNM )
                            CALL M3MSG2( MESG )

                        END IF

                    END IF

                END IF

            END DO

        END IF

        !.....  Search for tag species in the logical file
        IF( TDEV > 0 ) THEN
            DO N = 1, NTAG
                LFNTMP = TAG_LFN( N )      ! retriev logical file name from TAG_SPECIES

                IF( LFNTMP == NAM ) THEN

                    SPCTMP = TAG_SPC( N )     ! retrieve spcieces name from TAG_SPECIES

                    K = INDEX1( SPCTMP, NVARSET, VNAMESET )

                    IF( K <= 0 ) THEN
                        EFLAG = .TRUE.
                        MESG = 'ERROR: The species ' //    &
                            TRIM( TAG_SPC(N) )//' you want to tag '//    &
                            'is not available from file '//TRIM(NAM)
                        CALL M3MSG2( MESG )
                    END IF

                END IF

            END DO

        END IF

        !.....  Search for adj factor species in the logical file
        IF( ADEV > 0 ) THEN

            DO J = 1, NADJ

                LFNTMP = ADJ_LFN( J )         ! retriev logical file from ADJ_FACS

                IF( LFNTMP == NAM ) THEN
                    SPCTMP = ADJ_SPC( J )                     ! retrieve spcieces name from ADJ_FACS

                    IF( TDEV > 0 ) THEN

                        !....  Search tagged species for the current file
                        LFNSPC = TRIM( LFNTMP ) // '~' // TRIM( SPCTMP )
                        TAG = INDEX1( LFNSPC, NTAG, TAG_LFNSPCTAG )

                        !.....  Assign adjustment factor for the current species
                        IF( TAG > 0 ) SPCTMP = TAG_SPC( TAG )

                    END IF

                    K = INDEX1( SPCTMP, NVARSET, VNAMESET )
                    IF( K <= 0 ) THEN
                        EFLAG = .TRUE.
                        MESG = 'ERROR: The species '//TRIM(ADJ_SPC(J))//    &
                            ' you want to adjust is not available in the '//    &
                            TRIM(NAM) // ' file'
                        CALL M3MSG2( MESG )
                    END IF

                END IF

            END DO

        END IF

        !.....  Compare all other time steps back to first file.
        !.....  They must match exactly.
        IF( TSTEP3D /= TSTEP ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Time step', TSTEP3D,    &
                   'in file "' // TRIM( NAM ) //    &
                   '" is inconsistent with first file value of', TSTEP
            CALL M3MSG2( MESG )
        END IF

        !.....  Compare all other grids back to first grid.
        !.....  They must match exactly.
        WRITE( FDESC, '(A,I3.3)' ) 'FILE', F
        TFLAG = .FALSE.
        CALL CHKGRID( FDESC, 'GRID', 0, TFLAG )

        IF( TFLAG ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: File "' // TRIM( NAM ) //    &
              '" (NX,NY)  : (', NCOLSA( F ), ',', NROWSA( F ), ')'//    &
              CRLF() // BLANK10 // 'is inconsistent with first ' //    &
              'file (NX,NY) : (', NCOLS, ',', NROWS, ')'
            CALL M3MSG2( MESG )
        END IF

        !.....  Compare layer structures for 3-d files. The number of layers do
        !       not need to match, but the layer structures do need to match.
        NLAYS = MAX( NLAYS, NLAYSA( F ) )
        IF( NLAYSA( F ) .GT. 1 ) THEN
            LFLAG = .TRUE.

                !.....  For the first file that is 3-d, initialize output layer
                !       structure
            IF ( FIRST3D ) THEN

                NLAYS = NLAYSA( F )
                VGTYP = VGTYP3D
                VGTOP = VGTOP3D
                VGLVS( 0:NLAYS ) = VGLVS3D( 0+VLB:NLAYS+VLB )                   ! array
                FIRST3D = .FALSE.

                !.....  For additional 3-d files, compare the layer structures
            ELSE

                !.....  Check vertical type
                IF( VGTYP3D .NE. VGTYP ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 ) 'ERROR: Vertical ' //    &
                           'coordinate type', VGTYP3D,    &
                           'in file "'// TRIM( NAM ) //    &
                           '" is inconsistent with first 3-d'//    &
                           'file value of', VGTYP
                    CALL M3MSG2( MESG )
                END IF

                !.....  Check vertical top
                IF( VGTOP3D .NE. VGTOP ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 ) 'ERROR: Vertical top value', VGTOP3D,    &
                           'in file "'// TRIM( NAM ) //    &
                           '" is inconsistent with first 3-d'//    &
                           'file value of', VGTOP
                    CALL M3MSG2( MESG )
                END IF

                !.....  Loop through layers of current file F
                DO NL = 0, NLAYSA( F )

                    !.....  For layers that are common to this file and previous files
                    IF( NL .LE. NLAYS ) THEN

                        IF( VGLVS3D( NL+VLB ) .NE. VGLVS( NL )) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG, 94020 ) 'ERROR: Layer', NL,    &
                              'in file "'// TRIM( NAM ) //    &
                              '" with level value', VGLVS3D(NL+VLB),    &
                              CRLF()//BLANK10//'is inconsistent '//    &
                              'with first file value of', VGLVS(NL)
                            CALL M3MSG2( MESG )
                        END IF

                    !.....  Add additional layers from current file to output
                    !       layer structure
                    ELSE
                        VGLVS( NL ) = VGLVS3D( NL+VLB )
                    END IF

                END DO                ! End checking layers

                !.....  Reset the global number of layers as the maximum between
                !       the current file and all previous files
                NLAYS = MAX( NLAYS, NLAYSA( F ) )

            END IF        ! End first 3-d file or not
        END IF            ! End 3-d files

    END DO                ! End loop through files

    !.....  Give error message and end program unsuccessfully
    IF( EFLAG ) THEN
        MESG = 'Inconsistent time step, grid, species, or layers '//    &
                'among the files    !'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.....  Check that environment settings are consistent with files
    IF( MRGDIFF ) THEN
        IF( TSTEP /= G_TSTEP ) THEN
            WRITE( MESG,94010 ) 'ERROR: Value for G_TSTEP ',    &
                G_TSTEP, 'is inconsistent with the time step' //    &
                CRLF() // BLANK10 // 'of the input files', TSTEP
            CALL M3MSG2( MESG )

            MESG = 'Inconsistent environment settings'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
    END IF

    !.....  Deterimine output date, time, and number of time steps
    SDATE = G_SDATE
    STIME = G_STIME
    NSTEPS = G_NSTEPS
    CALL SETOUTDATE( SDATE, STIME, NSTEPS, NFILE, SDATEA,    &
                     STIMEA, DURATA, FNAME, MRGDIFF, USEFIRST )

    !.....  Build master output variables list
    NVOUT = 0

    !.....  Loop through input files and build an output variable list
    !.....  Loop through variables in the files
    DO F = 1, NFILE
    DO V = 1, NVARSA( F )

        VNM = VNAMEA( V,F )

        !.....  Look for variable name in output list
        K = INDEX1( VNM, NVOUT, OUTNAM  )              ! look in output list

        !.....  If its not in the output list, add it
        IF( K .LE. 0 ) THEN
            NVOUT = NVOUT + 1
            INDXN ( NVOUT ) = NVOUT
            OUTNAM( NVOUT ) = VNM
            VDESCU( NVOUT ) = VDESCA( V,F )
            VUNITU( NVOUT ) = VUNITA( V,F )

        !.....  If variable is in the output list, check the units
        ELSE
            IF ( VUNITA( V,F ) .NE. VUNITU( K ) ) THEN
                EFLAG = .TRUE.
                L  = LEN_TRIM( VNM )
                L1 = LEN_TRIM( VUNITA( V,F ) )
                L2 = LEN_TRIM( VUNITU( K )   )
                WRITE( MESG,94010 ) 'ERROR: Variable "' // TRIM( VNM ) //       &
                       '" in file', F, 'has units "'// TRIM( VUNITA( V,F ) ) // &
                       '"' // CRLF() // BLANK10 //    &
                       'that are inconsistent with a '//    &
                       'previous file that had units "' //    &
                       TRIM( VUNITU( K ) )// '" for this variable'
                CALL M3MSG2( MESG )

            END IF      ! End check of units

        END IF      ! End variable in output list already or not

    END DO          ! End loop through variables in this file
    END DO              ! End loop through files.

    !.....  Give error message and end program unsuccessfully
    IF( EFLAG ) THEN
        MESG = 'Inconsistent units for common variables among '//    &
               'the files    !'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.....  Sort output variables into alphabetical order
    CALL SORTIC( NVOUT, INDXN, OUTNAM )

    !.....  Set up for opening output file...
    !.....  Get grid information
    ICNTFIL = ALLFILES
    IF( NFILES( 1 ) .EQ. 1 ) ICNTFIL = 1       ! send ALLFILES if more than one file, send 1 otherwise

    IF( .NOT. DESCSET( IONAME( 1 ), ICNTFIL ) ) THEN
        MESG = 'Could not get description of file "'  //    &
                FNAME( 1 )( 1:LEN_TRIM( FNAME(1) ) ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ENDIF

    SDATE3D = SDATE
    STIME3D = STIME
    NVARS3D = NVOUT

    !.....  Set up layer structure for output file
    NLAYS3D = NLAYS
    VGTOP3D = VGTOP
    VGTYP3D = VGTYP
    VGLVS3D = 0.     ! initialize array
    DO NL = 0, NLAYS
        VGLVS3D( NL+VLB ) = VGLVS( NL )
    END DO

    !.....  Set the EV name for met data description
    MESG = 'Setting for the environment variable name for meta ' //    &
           'file description for output file'
    CALL ENVSTR( 'MRG_FILEDESC', MESG, ' ', MRGFDESC, IOS  )
    FDESC3D( 1 ) = 'Merged emissions output file from Mrggrid'

    IF( IOS >= 0 ) THEN
        MESG = 'Use this meta file description for output file'
        CALL ENVSTR( MRGFDESC, MESG, ' ', METADESC, IOS )
        FDESC3D( 1 ) = METADESC
    END IF

    !....  Update variable names in sorted order, and also
    !....  set up logical arrays for which files have which species
    DO V = 1, NVOUT
        VNM = OUTNAM( INDXN( V ) )

        VNAME3D( V ) = VNM                  ! store sorted output vars, etc.
        VDESC3D( V ) = VDESCU( INDXN( V ) )
        UNITS3D( V ) = VUNITU( INDXN( V ) )
        VTYPE3D( V ) = M3REAL

        DO F = 1, NFILE
            LVOUTA( V,F ) = .FALSE.
            J = INDEX1( VNM, NVARSA( F ), VNAMEA( 1,F ) )
            IF( J .GT. 0 ) LVOUTA( V,F ) = .TRUE.
        END DO
    END DO

    !.....  Allocate memory for the number of grid cells and layers
    NGRID = NROWS * NCOLS
    ALLOCATE( E2D( NGRID ), STAT=IOS )
    CALL CHECKMEM( IOS, 'E2D', PROGNAME )
    ALLOCATE( EOUT( NGRID, NLAYS ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EOUT', PROGNAME )

    IF( ETFLAG ) THEN    ! 3-d emissions table
        ALLOCATE( E3D( NGRID,NTEMPS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'E3D', PROGNAME )
        E3D = 0.0
    END IF

    !.....  Prompt for and open output file
    ONAME = PROMPTMFILE(    &
            'Enter logical name for MERGED GRIDDED OUTPUT file',    &
            FSUNKN3, 'OUTFILE', PROGNAME )

    !.....  Prompt for and open report file
    IF( MRGDIFF ) THEN
        RDEV = PROMPTFFILE(    &
               'Enter logical name for the MRGGRID REPORT file',    &
               .FALSE., .TRUE., 'REPMRGGRID', PROGNAME )

        !.....  Write header line to report
        WRITE( RPTLINE,93010 ) 'Output date'
        WRITE( RPTCOL,93010 ) 'Output time'
        RPTLINE = TRIM( RPTLINE ) // RPTCOL

        DO F = 1, NFILE
            NAM = FNAME( F )
            WRITE( RPTCOL,93010 ) TRIM( NAM ) // ' date'
            RPTLINE = TRIM( RPTLINE ) // RPTCOL
        END DO

        WRITE( RDEV,93000 ) TRIM( RPTLINE )
    END IF

    !.....  Warning missing logical file names from the ADJ_FACS list
    LNAM = ' '
    IF( ADEV > 0 ) THEN
        DO I = 1, NADJ
            NAM = ADJ_LFN( I )
            L = INDEX1( NAM, NFILE, FNAME )
            IF( L <= 0 ) THEN
                MESG = 'WARNING: The logical file '//TRIM(NAM) //    &
                    ' in the adjustment file (ADJ_FACS) is not '//    &
                    'found in the FILELIST on DATE : '//    &
                    MMDDYY(SDATE)
                IF( LNAM /= NAM ) THEN
                    CALL M3MSG2( MESG )
                    LNAM = NAM
                END IF
            END IF
        END DO
    END IF

    !.....  Warning missing logical file names from the TAG_SPECIES list
    LNAM = ' '
    IF( TDEV > 0 ) THEN
        DO I = 1, NTAG
            NAM = TAG_LFN( I )
            L = INDEX1( NAM, NFILE, FNAME )
            IF( L <= 0 ) THEN
                MESG = 'WARNING: The logical file '//TRIM(NAM) //    &
                    ' in the tagging file (TAG_SPECIES) is not '//    &
                    'found in the FILELIST on DATE : '//    &
                    MMDDYY(SDATE)
                IF( LNAM /= NAM ) THEN
                    CALL M3MSG2( MESG )
                    LNAM = NAM
                END IF
            END IF
        END DO
    END IF

    !.....  Allocate arrays that will store overall daily/gridded total emissinos by species
    ALLOCATE( BEFORE_SPC( NVOUT ),    &
               AFTER_SPC( NVOUT ), STAT=IOS )
    CALL CHECKMEM( IOS, 'AFTER_SPC', PROGNAME )
    BEFORE_SPC = 0.0
    AFTER_SPC  = 0.0

    !.....  Loop through hours
    JDATE = SDATE
    JTIME = STIME
    FACS  = 1.0

    IF( METFLAG ) THEN
        ALLOCATE( METRWC( NGRID ),    &
                  METNH3( NGRID ),    &
                  NOXGAS( NGRID ),    &
                  NOXDIS( NGRID ),    &
                     ID1( NGRID ),    &
                     ID2( NGRID ),    &
                      TA( NGRID ),    &
                      QV( NGRID ),    &
                  TABASE( NGRID ),    &
                      RA( NGRID ),    &
                  RABASE( NGRID ),    &
                      RC( NGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RC', PROGNAME )
        METRWC = 1.0
        METNH3 = 1.0
        NOXGAS = 1.0
        NOXDIS = 1.0
        ID1 = 0
        ID2 = 0
        TA = 0.0
        TABASE = 0.0
        RA = 0.0
        RABASE = 0.0
        RC = 0.0
        IF( .NOT. OPEN3( 'METCRO2D', FSREAD3, PROGNAME ) ) THEN
            MESG = 'Could not open METCRO2D meteorology file '
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
        END IF
        IF( NOXFLAG ) THEN
            QV = 0.0
            IF( .NOT. OPEN3( 'METCRO3D', FSREAD3, PROGNAME ) ) THEN
                MESG = 'Could not open METCRO3D meteorology file '
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            END IF
        END IF

    END IF

    DO T = 1, NSTEPS

    !.....  Apply the met adjustment for mobile, NH3 and RWC sources
        IF( METFLAG ) THEN
    !.....  Read current meteorology file
            IF( .NOT. READ3( 'METCRO2D', TVARNAME, 1,    &
                              JDATE, JTIME, TA ) ) THEN
                MESG = 'Could not read ' // TRIM( TVARNAME ) //    &
                       ' from METCRO2D input file'
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            END IF
            IF( T == 1 ) TABASE = TA

            IF( NH3FLAG ) THEN
                IF( .NOT. READ3( 'METCRO2D', RANAME, 1,    &
                                  JDATE, JTIME, RA ) ) THEN
                    MESG = 'Could not read aerodynamic resistance'//    &
                           ' variable from METCRO2D input file'
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF
                IF( T == 1 ) RABASE = RA

                IF( .NOT. READ3( 'METCRO2D', RCNAME, 1,    &
                                  JDATE, JTIME, RC ) ) THEN
                    MESG = 'Could not read precipitation (cm)'//    &
                           ' variable from METCRO2D input file'
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF
            END IF

            IF( MOBFLAG .AND. NOXFLAG ) THEN
                IF( .NOT. READ3( 'METCRO3D', QVNAME,1,    &
                                  JDATE, JTIME, QV ) ) THEN
                    MESG = 'Could not read QV specific humidity ' //    &
                           'variable from METCRO3D input file'
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF
            END IF

            METRWC = 1.0
            METNH3 = 1.0
            IF( RWCFLAG ) THEN

                DO S = 1, NGRID      ! loop over grid cells
                    TEMPVAL = 1.8 * TA( S ) - 459.67      ! K --> F
    !       zero-out the RWC emissions if ambient temp is higher than RWC threshold temp.
                    IF( TEMPVAL > RWC_TEMP ) METRWC( S ) = 0.2     !
                END DO
            ENDIF

            IF( NH3FLAG ) THEN      ! BASH NH3 calculation
                DO S = 1, NGRID      ! loop over grid cells
                    BASE =(161500.0/TABASE(S))*EXP(-10380.0/TABASE(S))
                    TEMPVAL=(161500.0/TA(S))*EXP(-10380.0/TA(S))
                    METNH3( S ) = 1.0 + 0.2 * ( (TEMPVAL-BASE) / BASE )
                    IF( RC( S ) > 0.1 ) METNH3(S) = METNH3(S) * 0.5      !  Reduce the NH3 emission by 50% if there is rain
                END DO
            ENDIF

            IF( MOBFLAG ) THEN
                NOXGAS = 1.0
                NOXDIS = 1.0

                DO S = 1, NGRID      ! loop over grid cells

                    TEMPVAL = ( TA( S ) - CTOK ) * CTOF + 32.      ! convert K to F degree

                    IDX1 = 0
                    IDX2 = 0
                    DO L1 = NTEMPS, 1, -1
                        IF( TEMPVAL < TEMPS( L1 ) ) CYCLE

                        IF( TEMPVAL == TEMPS( L1 ) ) THEN
                            IDX1 = L1
                            IDX2 = L1
                        ELSE
                            IDX1 = L1
                            IDX2 = L1 + 1
                        END IF
                        EXIT
                    END DO

                    IF( IDX1 == 0 ) THEN
                        IF( TEMPVAL >= (TEMPS( 1 ) - 50.)) THEN
                            IDX1 = 1
                            IDX2 = 1
                        END IF
                    END IF

                    IF( IDX2 > NTEMPS ) THEN
                        IF( TEMPVAL <= (TEMPS( NTEMPS ) + 50.) ) THEN
                            IDX1 = NTEMPS
                            IDX2 = NTEMPS
                        END IF
                    END IF

                    ID1( S ) = IDX1
                    ID2( S ) = IDX2

    !.....  NOx humidity correction
                    IF( NOXFLAG ) THEN
                        QV( S ) = QV( S ) * 1000.0     ! convert from kg water/kg dry air to g/kg
                        IF( NOXADJEQS ) THEN       ! Use newer MOVES3 NOx humidity correction eqs
                            A = MIN( QV( S ), 17.71 )
                            B = MAX( 3.0, A )
                            NOXGAS( S ) = 1.0 - 0.0329 * ( B - 10.71 )       ! NOx humidity correction for Gasoline fuel
                            QVMOL = QV( S ) * 0.001607524      ! convert from g of water/kg of dry air to moles of water/moles of dry air
                            A = MIN( QVMOL, 0.035 )
                            B = MAX( 0.002, A )
                            NOXDIS( S ) = 1.0 / ( 9.953 * B  + 0.832 )      ! Nox humidity correction for Diesel fuel
                        ELSE       ! Use older MOVES NOx humidity correction eqs
                            QV( S ) = QV( S ) * 7.0
                            A = MIN( QV( S ), 124.0)
                            B = MAX( 21.0, A ) - 75.0
                            NOXGAS( S ) = 1.0 - ( B * 0.0038 )     ! NOx humidity correction for Gasoline fuel
                            NOXDIS( S ) = 1.0 - ( B * 0.0026 )     ! NOx humidity correction for Diesel fuel
                        END IF
                    END IF

                END DO      !  end parallel loop on S

            END IF
        END IF

    !.....  Loop through species
        DO V = 1, NVOUT

            VNM = VNAME3D( V )

    !.....  Output array
            EOUT = 0.       ! array
            LNOX = 0

            DO F = 1, NFILE

    !.....  Set read date
                IF( MRGDIFF ) THEN
                    IF( USEFIRST( F ) ) THEN
                        DUMMY = 0
                        STEPS = SEC2TIME(    &
                                SECSDIFF(    &
                                  SDATE, DUMMY, JDATE, DUMMY ) )
                        RDATE = SDATEA( F )
                        CALL NEXTIME( RDATE, DUMMY, STEPS )
                    ELSE
                        RDATE = JDATE
                    END IF
                ELSE
                    RDATE = JDATE
                END IF

    !.....  Set tmp variables
                NL     = NLAYSA( F )       ! no of layers
                NAM    = FNAME ( F )       ! input file name
                IO_NAM = IONAME( F )       ! retrieve 16 char ioapi input file name
                METNAM = METADJ( F )       ! metadj variable name
                NOXNAM = NOXADJ( F )       ! NOx adustment eqs names
                LNOX   = LEN_TRIM( NOXNAM )

    !.....  Search adjustment factor for the current file
                LFNSPC = TRIM( NAM ) // '~' // TRIM( VNM )

    !.....  Assign adjustment factor for the current species
                ADJ = INDEX1( LFNSPC, NADJ, ADJ_LFNSPC )
                IF( ADJ > 0 ) THEN
                    FACS = ADJ_FACTOR( ADJ )
                    WRITE( MESG,93011 )'Apply adjustment factor' ,    &
                        FACS, ' to the '  // TRIM( VNM ) //    &
                        ' species from the '//TRIM( NAM )// ' file'
                    CALL M3MSG2( MESG )
                ELSE
                    FACS = 1.0
                END IF

    !.....  Search tagged species for the current file
                IF( TDEV > 0 ) THEN
                    TAG = INDEX1( LFNSPC, NTAG, TAG_LFNSPCTAG )
                    IF( TAG > 0 ) THEN
                        TVNM = TAG_SPC( TAG )
                    ELSE
                        TVNM = VNM
                    END IF
                ELSE
                    TVNM = VNM

                END IF

    !.....  If file has species, read (do this for all files)...
                IF( LVOUTA( V,F ) ) THEN

                    ICNTFIL = ALLFILES
                    IF( NFILES( F ) .EQ. 1 ) ICNTFIL = 1       ! send ALLFILES if more than one file, send 1 otherwise

    !.....  Read, and add
                    DO K = 1, NL
                        IF( .NOT.    &
                            READSET( IO_NAM, TVNM, K, ICNTFIL, RDATE, JTIME, E2D )) THEN
                            MESG = 'Could not read "' // TVNM // '" from file "' // TRIM( NAM ) // '".'
                            CALL M3EXIT( PROGNAME, RDATE, JTIME, MESG, 2 )
                        ENDIF

    !.....  Apply the met adjustment for mobile, NH3 and RWC sources
                        IF( LEN_TRIM( METNAM ) > 0 ) THEN

                            IF( METNAM == 'RWC' ) THEN
                                E2D = METRWC * E2D

                            ELSE IF(  METNAM == 'AGNH3' ) THEN
                                E2D = METNH3 * E2D

                            ELSE IF( METNAM == 'METEMIS' ) THEN

                                IF( .NOT. READSET( IO_NAM, TVNM, -1, ICNTFIL,    &
                                                   RDATE, JTIME, E3D )) THEN
                                    MESG = 'Could not read "' // TVNM // '" from file "' // TRIM( NAM )// '".'
                                    CALL M3EXIT( PROGNAME, RDATE, JTIME, MESG, 2 )
                                ENDIF

                                DO S = 1, NGRID      ! loop over grid cells
                                    TEMPVAL = ( TA( S ) - CTOK ) * CTOF + 32.      ! convert K to F degree
                                    IDX1 = ID1( S )
                                    IDX2 = ID2( S )
                                    RATIO = 1.0
                                    IF( IDX1 .NE. IDX2 ) THEN
                                        RATIO = ( TEMPVAL - TEMPS(IDX1) ) /               &    ! temp interpolation fraction
                                                ( TEMPS(IDX2) - TEMPS(IDX1) )
                                    END IF
                                    E2D(S) = E3D(S,IDX1) + RATIO * (E3D(S,IDX2)-E3D(S,IDX1))

    !.....  NOx humidity correction
                                    IF( LNOX > 0 ) THEN
                                        IF( TVNM=='NO' .OR. TVNM=='NO2' .OR. TVNM=='HONO' ) THEN
                                            RATIO = 1.0
                                            IF( NOXNAM == 'GASOLINE') RATIO =  NOXGAS( S )    ! Applying NOx humidity correction for gasoline
                                            IF( NOXNAM == 'DIESEL'  ) RATIO =  NOXDIS( S )    ! Applying NOx humidity correction for diesel
                                            E2D( S ) = RATIO * E2D( S )
                                        END IF
                                    END IF
                                END DO
                            END IF
                        END IF

    !.....  Logical file specific summary
                        IF( ADJ > 0 ) THEN
                            SUME2D = SUM( E2D(1:NGRID) )
                            BEFORE_ADJ( ADJ ) = BEFORE_ADJ( ADJ ) + SUME2D
                            AFTER_ADJ ( ADJ ) = AFTER_ADJ ( ADJ ) + SUME2D*FACS
                            BEFORE_SPC( V )   = BEFORE_SPC( V )   + SUME2D
                            AFTER_SPC ( V )   = AFTER_SPC ( V )   + SUME2D*FACS
                        END IF

                        EOUT( 1:NGRID,K ) = EOUT( 1:NGRID,K ) +    &
                                             E2D( 1:NGRID ) * FACS

                    END DO      !  end loop on K

                END IF      ! if pollutant is in this file

            END DO          ! parallel loop over input files

    !.....  Build report line if needed
            IF( MRGDIFF .AND. V == 1 ) THEN
                WRITE( RPTLINE,93020 ) JDATE
                WRITE( RPTCOL,93020 ) JTIME
                RPTLINE = TRIM( RPTLINE ) // RPTCOL
                DO F = 1, NFILE
                    IF( USEFIRST( F ) ) THEN
                        DUMMY = 0
                        STEPS = SEC2TIME( SECSDIFF( SDATE, DUMMY, JDATE, DUMMY ) )
                        RDATE = SDATEA( F )
                        CALL NEXTIME( RDATE, DUMMY, STEPS )
                    END IF
                    WRITE( RPTCOL,93020 ) RDATE
                    RPTLINE = TRIM( RPTLINE ) // RPTCOL
                END DO
            END IF

    !.....  Write species/hour to output file
            IF( .NOT. WRITE3( ONAME, TVNM, JDATE, JTIME, EOUT )) THEN
                MESG = 'Could not write "'// TVNM// '" to file "'//    &
                        ONAME( 1:LEN_TRIM( ONAME ) ) // '".'
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            END IF

        END DO       ! loop through variables

    !.....  Write this time step to report
        IF( MRGDIFF ) THEN
            WRITE( RDEV,93000 ) TRIM( RPTLINE )
        END IF

        CALL NEXTIME( JDATE, JTIME, TSTEP )

    END DO       ! loop through timesteps

    !....  Write summary of sector specific factor adjustment output
    !       Columns: Date, Sector, Species, value before, value after, ratio of before/after
    !       Later we can add the total amount of the adjusted species summed accross all of the input files
    !       and the total amount of the adjusted species in the output file.

    !.....  Write header line to report
    DO F = 1, NADJ

        VNM = ADJ_SPC( F )     ! species name
        NAM = ADJ_LFN( F )     ! logical file name
        FACS= ADJ_FACTOR( F )       ! adjustment factor

        IF( BEFORE_ADJ( F ) == 0.0 ) CYCLE

        RATIO = ( AFTER_ADJ( F ) / BEFORE_ADJ( F ) )

        REPFMT = "(I8,2(',',A),',',F10.6,',',"

    !.....  Define the format of real values
        CALL GET_FORMAT( VNM, BEFORE_ADJ( F ), EFMT )
        REPFMT = TRIM( REPFMT ) // TRIM( EFMT )

        CALL GET_FORMAT( VNM, AFTER_ADJ( F ), EFMT )
        REPFMT = TRIM( REPFMT ) // TRIM( EFMT )
        REPFMT = TRIM( REPFMT ) // "F10.6)"

        WRITE( RPTLINE,REPFMT ) SDATE, NAM, VNM, FACS,    &
                              BEFORE_ADJ( F ), AFTER_ADJ( F ), RATIO

        IF( ADEV > 0 .AND. RATIO /= 1.0 ) THEN
            WRITE( ODEV,93000 ) TRIM( RPTLINE )
        ENDIF
    END DO

    !.....  Write header line to overall summary report
    DO V = 1, NVOUT

        VNM   = VNAME3D( V )     ! species name

        IF( BEFORE_SPC( V ) == 0.0 ) CYCLE

        RATIO = ( AFTER_SPC( V ) / BEFORE_SPC( V ) )

        REPFMT = "( I8,',',A,',',"

    !.....  Define the format of real values
        CALL GET_FORMAT( VNM, BEFORE_SPC( V ), EFMT )
        REPFMT = TRIM( REPFMT ) // TRIM( EFMT )

        CALL GET_FORMAT( VNM, AFTER_SPC( V ), EFMT )
        REPFMT = TRIM( REPFMT ) // TRIM( EFMT )
        REPFMT = TRIM( REPFMT ) // "F10.6)"

        WRITE( RPTLINE,REPFMT ) SDATE, VNM, BEFORE_SPC(V),    &
                                AFTER_SPC(V),RATIO

        IF( SDEV > 0 .AND.  RATIO /= 1.0 ) THEN
            WRITE( SDEV,93000 ) TRIM( RPTLINE )
        END IF
    END DO

    !..... Normal Completion
    CALL M3EXIT( PROGNAME, 0, 0, ' ', 0)

    !******************  FORMAT  STATEMENTS   ******************************

    !.....   Formatted file I/O formats.... 93xxx

92000 FORMAT(  A,',',A,',',A  )

93000 FORMAT(  A )

93010 FORMAT( A15 )

93011 FORMAT(  A, F8.5, A )

93020 FORMAT( I15 )

93030 FORMAT( A,I2.2 )

    !.....   Internal buffering formats.... 94xxx

94010 FORMAT( 10( A, :, I7, :, 1X ) )

94020 FORMAT( A, :, I3, :, 1X, 10 ( A, :, F8.5, :, 1X ) )


    !*****************  INTERNAL SUBPROGRAMS  ******************************

CONTAINS

    !----------------------------------------------------------------------

    !.....  This internal subprogram determines the format to output
    !       emission values
    SUBROUTINE GET_FORMAT( VBUF, VAL, FMT )

        !.....  Subroutine arguments
        CHARACTER(*), INTENT (IN) :: VBUF
        REAL        , INTENT (IN) :: VAL
        CHARACTER(*), INTENT(OUT) :: FMT

        !----------------------------------------------------------------------

        !.....  Value is too large for
        IF( VAL .GT. 999999999. ) THEN
            FMT = "E10.3,',',"

            WRITE( MESG,95020 ) 'WARNING: "' //TRIM( VBUF )//    &
                '"Emissions value of', VAL, CRLF()// BLANK10//    &
                '" is too large for file format, so writing ' //    &
                'in scientific notation for source'

            CALL M3MESG( MESG )

        ELSE IF( VAL .GT. 99999999. ) THEN
            FMT = "F10.0,',',"

        ELSE IF( VAL .GT. 9999999. ) THEN
            FMT = "F10.1,',',"

        ELSE IF( VAL .GT. 999999. ) THEN
            FMT = "F10.2,',',"

        ELSE IF( VAL .GT. 0. .AND. VAL .LT. 1. ) THEN
            FMT = "E10.4,',',"

        ELSE
            FMT = "F10.3,',',"

        END IF

        RETURN

95020   FORMAT( 10( A, :, E12.5, :, 1X ) )

    END SUBROUTINE GET_FORMAT

    !*****************  INTERNAL SUBPROGRAMS  ******************************

    !----------------------------------------------------------------------
    !.....  This internal subprogram determines to store a list of
    !       adjustment factos  from ADJ_FACS input file.
    SUBROUTINE READ_ADJ_FACS( F )

        !.....  Subroutine arguments
        INTEGER,     INTENT( OUT ) :: F

        CHARACTER( 32 ) SEGMENT( 6 )              ! line parsing array

        !.....

        !.....  Loop through input files and open them
        IREC = 0
        F = 0
        DO

            !.....  Read file names - exit if read is at end of file
            READ( ADEV, 93000, END = 30, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                WRITE( MESG,94010 )    &
                    'I/O error', IOS,    &
                    'reading adustment factor file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                CYCLE
            END IF

            !....  Skip blank and comment lines
            IF ( BLKORCMT( LINE ) ) CYCLE

            !....  Get line
            CALL PARSLINE( LINE, 3, SEGMENT )

            CALL UPCASE( SEGMENT( 1 ) )               ! species name
            CALL UPCASE( SEGMENT( 2 ) )               ! logical file name

            !.....  Search adjustment factor for the current file
            NAM = TRIM( SEGMENT( 2 ) )

            L = INDEX1( NAM, NFILE, FNAME )

            !.....  Skip EMF-specific header line
            IF( L <= 0 .AND. .NOT. CHKREAL( SEGMENT( 3 ) ) ) CYCLE

            F = F + 1

            ADJ_SPC( F ) = TRIM( SEGMENT( 1 ) )
            ADJ_LFN( F ) = TRIM( SEGMENT( 2 ) )
            ADJ_LFNSPC( F ) = TRIM( SEGMENT( 2 ) ) // '~' // TRIM( SEGMENT( 1 ) )
            ADJ_FACTOR( F ) = STR2REAL( SEGMENT( 3 ) )

            IF( ADJ_FACTOR( F ) < 0 ) THEN
                MESG = 'WARNING: ' // TRIM( ADJ_SPC(F) ) //     &
                   ' emissions from the ' //TRIM(NAM)//         &
                   ' file will be substracted due to a negative' //    &
                   ' adjustment factor'
                CALL M3MSG2( MESG )

            ELSE IF( ADJ_FACTOR( F ) == 0 ) THEN
                MESG = 'WARNING: ' // TRIM( ADJ_SPC(F) ) //    &
                   ' emissions from the ' //TRIM(NAM)// ' file' //    &
                   ' will be zero due to a zero adjustment factor'
                CALL M3MSG2( MESG )

            END IF

        END DO

30      CONTINUE

        RETURN

93000   FORMAT(  A )

94010   FORMAT( 10( A, :, I7, :, 1X ) )

    END SUBROUTINE READ_ADJ_FACS

    !*****************  INTERNAL SUBPROGRAMS  ******************************

    !----------------------------------------------------------------------
    !.....  This internal subprogram determines to store a list of
    !       adjustment factos  from ADJ_FACS input file.
    SUBROUTINE READ_TAG_SPECIES( F )

        !.....  Subroutine arguments
        INTEGER,     INTENT( OUT ) :: F

        CHARACTER( 32 ) SEGMENT( 6 )              ! line parsing array

        !.....

        MESG = 'Enter logical name for a list of tagging species'
        TDEV = PROMPTFFILE( MESG,.TRUE.,.TRUE.,'TAG_SPECIES',PROGNAME )
        MXNTAG = GETFLINE( TDEV, 'List of tagging species' )

        !.....  Allocate memory for arrays that just depend on the maximum number
        !       of adjustment factors in ADJ_FACS input file.
        ALLOCATE( TAG_LFN( MXNTAG ),    &
                  TAG_SPC( MXNTAG ),    &
               TAG_LFNSPC( MXNTAG ),    &
               TAG_APPEND( MXNTAG ),    &
               TAG_LFNSPCTAG( MXNTAG ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TAG_LFNSPCTAG', PROGNAME )

        TAG_SPC = ' '
        TAG_LFN = ' '
        TAG_LFNSPC = ' '
        TAG_APPEND = ' '
        TAG_LFNSPCTAG = ' '

        !.....  Loop through input files and open them
        IREC = 0
        F = 0
        DO

            !.....  Read file names - exit if read is at end of file
            READ( TDEV, 93000, END = 40, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                WRITE( MESG,94010 )    &
                    'I/O error', IOS,  'reading APPEND_TAGS file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                CYCLE
            END IF

            !....  Skip blank and comment lines
            IF ( BLKORCMT( LINE ) ) CYCLE

            !....  Get line
            CALL PARSLINE( LINE, 3, SEGMENT )
            CALL UPCASE( SEGMENT( 1 ) )               ! logical file name
            CALL UPCASE( SEGMENT( 2 ) )               ! species name
            CALL UPCASE( SEGMENT( 3 ) )               ! tagging name

            !.....  Skip EMF-specific header line
            IF( TRIM( SEGMENT( 1 ) ) == 'SECTOR' ) THEN
                IF( TRIM( SEGMENT( 2 ) ) == 'SPECIES' ) THEN
                    IF( TRIM( SEGMENT( 3 ) ) == 'TAG' ) THEN
                        CYCLE
                    END IF
                END IF
            END IF

            !.....  store tagging species
            F = F + 1

            TAG_LFN( F ) = TRIM( SEGMENT( 1 ) )
            TAG_SPC( F ) = TRIM( SEGMENT( 2 ) )
            TAG_LFNSPC( F ) = TRIM( SEGMENT( 1 ) ) // '~' //    &
                              TRIM( SEGMENT( 2 ) )
            TAG_APPEND( F ) = TRIM( SEGMENT( 3 ) )
            TAG_LFNSPCTAG(F)= TRIM( SEGMENT( 1 ) ) // '~' //    &
                              TRIM( SEGMENT( 2 ) ) // '_' //    &
                              TRIM( SEGMENT( 3 ) )

            L1 = LEN_TRIM( SEGMENT( 2 ) )
            L2 = LEN_TRIM( SEGMENT( 3 ) )
            L  = L1 + L2 + 1

            IF( L > 16 ) THEN
                MESG = 'ERROR: Can not append a tag to the species '    &
                    // TRIM( TAG_SPC( F ) ) // ' from the ' //    &
                    TRIM( NAM ) // ' file due to exceeding ' //    &
                    'max size(16-char) of tagged species name.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        END DO

40      CONTINUE

        RETURN

93000   FORMAT(  A )

94010   FORMAT( 10( A, :, I7, :, 1X ) )

    END SUBROUTINE READ_TAG_SPECIES

END PROGRAM MRGGRID

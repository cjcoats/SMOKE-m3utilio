
PROGRAM SMK2EMIS

!***********************************************************************
!  program body starts at line 152
!
!  DESCRIPTION:
!       Use SMOKE NetCDF gridded input to produce UAM EMISSIONS input
!       files.  Outputs 2-d files only.
!
!  PRECONDITIONS REQUIRED:
!       M3IO single layer gridded data
!       Tested for lat-lon and UTM projections only
!       Tested for CAMx and UAM-V compatiblity
!       Not fully tested for UAM-IV compatiblity
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!       Models-3 I/O
!       GETYN, PROMPTFFILE, PROMPTMFILE, NEXTIME
!
!  REVISION  HISTORY:
!       Prototype  4/00 by JMV
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

!.........  MODULES for public variables
!.........  This module contains the information about the source category
    USE MODINFO, ONLY: CRL

    IMPLICIT NONE

!...........   INCLUDES:
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   EXTERNAL FUNCTIONS and their descriptions:
    INTEGER           , EXTERNAL :: GETFLINE

!...........   PARAMETERS and their descriptions:

    CHARACTER(16), PARAMETER :: PROGNAME = 'SMK2EMIS'   !  program name
    CHARACTER(50), PARAMETER :: CVSW =&
    &     '$Name SMOKEv5.0_Jun2023$' ! CVS release tag

!.........  Allocatable arrays
    CHARACTER(4), ALLOCATABLE :: SPCNAM ( :,: )  ! output species names
    REAL,         ALLOCATABLE :: EMIS   ( :,: )  ! input and output variable

    CHARACTER(NAMLEN3), ALLOCATABLE :: OUTNAMS( : ) ! output spec names

!.........  Static arrays
    CHARACTER(4)    CFTYPE( 10 )  ! UAM filename
    CHARACTER(4)    CFNOTE( 60 )  ! UAM note

    CHARACTER(NAMLEN3) SEGMENT( 2 )

!.........  Unit numbers and logical file names
    INTEGER :: LDEV          ! unit number for log file
    INTEGER :: MDEV          ! unit number for mapping file
    INTEGER :: ODEV          ! UAM file desigator for output file

    CHARACTER(16)  ENAME     !  NetCDF logical input file name

!.........  Other local variables

    INTEGER    I, J, K, L1, L2, HR   ! indices and counters

    INTEGER    BTIME         ! beginning time of output file (YYYYDDD)
    INTEGER    ETIME         ! ending time of output file (HHMMSS)
    INTEGER    IBD           ! start date for each output time period (YYDDD)
    INTEGER    IED           ! end date for output time period (YYDDD)
    INTEGER    IOS           ! i/o status
    INTEGER    IREC          ! line count
    INTEGER    IUTMZONE      ! UTM zone
    INTEGER :: IXY = 0       ! segment origin always 0
    INTEGER    JDATE         ! date variable (YYYYDDD)
    INTEGER    JTIME         ! time variable (HHMMSS)
    INTEGER :: NLAYS = 1     ! number of vertical layers in UAM file
    INTEGER    NLINE         ! no. lines input ascii file
    INTEGER :: NSEG = 1      ! number of segments always 1
    INTEGER    NSTEPS        ! number of time steps
    INTEGER    NSPECS        ! number of species
    INTEGER    NZLOWR, NZUPPR ! UAM-IV diffbreak layer variables
    INTEGER    NDAYS, NHRS
    INTEGER    SDATE         ! UAM file starting date
    INTEGER    SDATEP1       ! date variable plus one time step (YYYYDDD)
    INTEGER    STIME         ! UAM file starting hour
    INTEGER    STIMEP1       ! time variable plus one time step (HHMMSS)
    INTEGER    TSTEP         ! time step (HHMMSS)
    INTEGER    TZONE         ! time zone of output

    LOGICAL :: EFLAG   = .FALSE. ! true: error found
    LOGICAL :: MAPFLAG = .FALSE. ! true: use variable mapping input file

    REAL       BT            ! Beginning time for a time step period
    REAL       ET            ! Ending time for a time step period
    REAL       HTSUR, HTLOWR, HTUPPR  ! UAM-IV layer ht variables
    REAL       UTMX, UTMY    !
    REAL       XORIG, YORIG  ! Grid origin
    REAL       XCELL, YCELL  ! Grid resolution

    CHARACTER      TEMPCH  !  For integer conversion
    CHARACTER(10)  HDRKEY  !  UAM file name label for output file
    CHARACTER(44)  NOTEDEF !  UAM file note default
    CHARACTER(44)  UNOTE   !  UAM file note from env variable
    CHARACTER(60)  FNOTE   !  UAM file description
    CHARACTER(80)  LINE    !  line buffer
    CHARACTER(256) MESG    !  Temporary message array

!***********************************************************************
!   begin body of program SMK2EMIS

    LDEV = INIT3()

!.........  Write out copyright, version, web address, header info, and prompt
!           to continue running the program.

    CALL INITEM( LDEV, CVSW, PROGNAME )

!.........  Get environment variables that control this program
    MAPFLAG = ENVYN( 'SMK2EMIS_VMAP_YN', 'Use names mapping file',&
    &                 .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK2EMIS_VMAP_YN"', 2 )
    END IF

!.........  Set source category based on environment variable setting
    CALL GETCTGRY

!.........  Prompt for name of NetCDF input file

    ENAME = PROMPTMFILE(&
    &     'Enter logical name for SMOKE gridded input (NetCDF) file',&
    &      FSREAD3, CRL // 'GTS_L', PROGNAME )

!.........  Prompt for variable mapping file, if needed
    IF( MAPFLAG ) THEN

        MDEV = PROMPTFFILE(&
        &       'Enter logical name for VARIABLE NAME MAP file',&
        &       .TRUE., .TRUE., 'VNAMMAP', PROGNAME )

    END IF

!.........  Prompt for keyword name of output file
!.........  Only EMISSIONS supported at this point

    MESG = 'UAM file label for output file'
    CALL ENVSTR( 'FLABEL', MESG, 'EMISSIONS', HDRKEY, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "FLABEL"', 2 )
    END IF

!.........  Convert UAM keyword name to character format (as in UAM)
    CFTYPE = ' '
    DO J = 1, 10
        CFTYPE( J ) = HDRKEY( J:J )
    END DO

!.........  Prompt for name of output file

    ODEV = PROMPTFFILE(&
    &       'Enter logical name for ' // HDRKEY //&
    &       ' UAM output file ',&
    &       .FALSE., .FALSE., 'UAM_' // CRL // 'GTS', PROGNAME )

!.........  Get the time zone for output of the emissions
    TZONE = ENVINT( 'OUTZONE', 'Output time zone', 0, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "FLABEL"', 2 )
    END IF

!.........   Get setting from environment variables
!.........   for UTM and UAM4 variables

    IUTMZONE = ENVINT( 'UTM_ZONE', 'UTM zone ', -9 , IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UTM_ZONE"', 2 )
    END IF

    NZLOWR = ENVINT( 'UAM4_LAYBELOW', 'Layers below diffbreak', 0 , IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM4_LAYBELOW"', 2 )
    END IF

    NZUPPR = ENVINT( 'UAM4_LAYABOVE', 'Layers above diffbreak',&
    &                 0 , IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM4_LAYABOVE"', 2 )
    END IF

    MESG = 'Height of surface layer [meters]'
    HTSUR = ENVREAL( 'UAM4_HTSFC', MESG, 0., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM4_HTSFC"', 2 )
    END IF

    MESG = 'Min. ht. of cells between sfc and diffbreak [meters]'
    HTLOWR = ENVREAL( 'UAM4_HTLOWR', MESG, 0., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM4_HTLOWR"', 2 )
    END IF

    MESG = 'Min. ht. of cells between diffbreak and top [meters]'
    HTUPPR = ENVREAL( 'UAM4_HTUPPR', MESG, 0., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM4_HTUPPR"', 2 )
    END IF

!.........  Create FNOTE from grid name and UAM_NOTE env. variable
!.........  Convert FNOTE to character format (as in UAM)

    NOTEDEF = ' UAM gridded emissions from ' // PROGNAME

    MESG = 'UAM file note for output file'
    CALL ENVSTR( 'UAM_NOTE', MESG, NOTEDEF , UNOTE, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM_NOTE"', 2 )
    END IF

    FNOTE  = GDNAM3D // UNOTE ( 1 : 44 )

    CFNOTE = ' '
    DO J = 1, 60
        CFNOTE( J ) = FNOTE( J:J )
    END DO

!.........  Read NetCDF header information:
!.........  coordinate info, date, time, timestep , and variables

    IF ( .NOT. DESC3( ENAME ) ) THEN

        MESG = 'Could not get description of file "' //&
        &       TRIM( ENAME ) // '"'
        CALL M3EXIT( PROGNAME , 0, 0, MESG, 2 )

    END IF

    NSPECS = NVARS3D
    NSTEPS = MXREC3D
    TSTEP  = TSTEP3D
    JDATE  = SDATE3D
    JTIME  = STIME3D
    XORIG  = XORIG3D
    YORIG  = YORIG3D
    XCELL  = XCELL3D
    YCELL  = YCELL3D

!...  If time step other than hourly then exit
!...  Assumes UAM/CAMx wants hourly emissions for now

    IF ( TSTEP .NE. 10000 ) THEN

        WRITE( MESG,94010 ) 'ERROR: Gridded SMOKE emissions '//&
        &        'time-step not hourly :  TSTEP =  ', TSTEP
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

!.............  Prompt for starting date of UAM output file

    SDATE  = GETNUM( JDATE, 9999999, JDATE,&
    &                'Enter starting date (YYYYDDD)' )

    IF( SDATE .NE. JDATE ) JTIME = 0  ! Reset because now could be 0 -24

!.............  Prompt for starting time of UAM output file

    STIME  = GETNUM( JTIME, 235959, JTIME,&
    &               'Enter starting time (HHMMSS)' )
    BTIME  = INT( STIME / TSTEP )

!.............  Prompt for time number of time steps
!.............  Number of hours to output to UAM file

    I = SECSDIFF( JDATE, JTIME, SDATE, STIME ) / 3600
    NSTEPS = NSTEPS - I
    MESG   = 'Enter number of time steps'
    NSTEPS = GETNUM( 1, NSTEPS, NSTEPS, MESG )

!.........  Calculate ending date and time

    NDAYS = NSTEPS / 24
    NHRS  = MOD( NSTEPS , 24 )

    IF ( SDATE .LT. 2000000 ) THEN
        IBD = SDATE  - 1900000
    ELSE
        IBD = SDATE  - 2000000
    END IF

    IED = IBD + NDAYS
    ETIME = MOD( ( BTIME + NHRS - 1 ), 24 )

!.........  Allocate memory for and read remapping information
    ALLOCATE( OUTNAMS( NSPECS ), STAT=IOS )
    CALL CHECKMEM( IOS, 'OUTNAMS', PROGNAME )

!.........  Initialize names map as input names
    OUTNAMS( 1:NSPECS ) = VNAME3D( 1:NSPECS )

    IF ( MAPFLAG ) THEN

        NLINE = GETFLINE( MDEV, 'Variable mapping file' )

        CALL M3MSG2( 'Reading variable mapping file...' )

        IREC = 0
        DO I = 1, NLINE

            READ( MDEV, 93000, END=999, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF( LINE .EQ. ' ' ) CYCLE

            IF ( IOS .NE. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG, 94010 )&
                &    'I/O error', IOS,&
                &    'reading variable mapping file at line', IREC
                CALL M3MESG( MESG )
                CYCLE

            END IF

            CALL PARSLINE( LINE, 2, SEGMENT )

!.................  Find variable name in list of species
            K = INDEX1( SEGMENT( 1 ), NSPECS, VNAME3D )

!.................  If not found, write message
            IF( K .LE. 0 ) THEN

                MESG = 'WARNING: Ignoring "' //&
                &       TRIM( SEGMENT( 1 ) ) //&
                &       '" in variable mapping file'
                CALL M3MSG2( MESG )

!.................  If found, write message and reassign
            ELSE

                MESG = 'NOTE: Renaming variable "' //&
                &       TRIM( OUTNAMS( K ) ) // '" as "' //&
                &       TRIM( SEGMENT( 2 ) ) // '"'
                CALL M3MSG2( MESG )

                OUTNAMS( K ) = SEGMENT( 2 )

            END IF

        END DO

    END IF

    IF ( EFLAG ) THEN

        MESG = 'Problem reading mapping names file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

!.........  Allocate for species name in character format
!.........  Allocate memory for emissions
    ALLOCATE( SPCNAM( NSPECS, 10  ),&
    &            EMIS( NCOLS3D, NROWS3D ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SPCNAM', PROGNAME )

!.........  Convert output species names to character format
    DO J = 1, 10
        DO I = 1, NSPECS
            SPCNAM( I,J ) = OUTNAMS( I )( J:J )
        END DO
    END DO

!.........  Write UAM Binary file  .........................................

!....   Write UAM EMISSIONS header

    WRITE(ODEV,ERR=7010) CFTYPE, CFNOTE, NSEG, NSPECS, IBD,&
    &                     REAL( BTIME ), IED , REAL( ETIME )

    WRITE(ODEV,ERR=7010) UTMX, UTMY, IUTMZONE,&
    &                     XORIG, YORIG, XCELL,&
    &                     YCELL, NCOLS3D, NROWS3D, NLAYS,&
    &                     NZLOWR, NZUPPR, HTSUR, HTLOWR, HTUPPR


    WRITE(ODEV,ERR=7010) IXY, IXY, NCOLS3D, NROWS3D

!.........  Write out character form of output variable names

    WRITE(ODEV,ERR=7010) ((SPCNAM(J,I),I=1,10),J=1,NSPECS)

!.........  Initialize time and date for time step after STIME and SDATE

    SDATEP1 = SDATE
    STIMEP1 = STIME
    CALL NEXTIME (SDATEP1, STIMEP1, TSTEP)

!.........  Loop over time steps

    DO HR = 1, NSTEPS

!.........  Initialize output file special time parameters
!.........  Note 20th and 21st centuries only are supported

        BT  = REAL( STIME/10000 )
        ET  = BT + 1.
        IF ( SDATE .LT. 2000000 ) THEN
            IBD = SDATE   - 1900000
            IED = SDATEP1 - 1900000
        ELSE
            IBD = SDATE   - 2000000
            IED = SDATEP1 - 2000000
        END IF

        WRITE (*,93200) SDATE, STIME

!.............  Set UAM formatted times and dates

        WRITE(ODEV,ERR=7011) IBD, BT, IED, ET

!.............  Loop over output variables

        DO K = 1, NSPECS

!.....................  Read input file for time and species of interest

            IF ( .NOT. READ3( ENAME, VNAME3D( K ), 1,&
            &                  SDATE, STIME, EMIS )) THEN
                CALL M3EXIT( PROGNAME , 0, 0,&
                &           'Error reading ' // VNAME3D( K ) //&
                &           ' from file ' // ENAME , 2 )
            END IF

!.............  Write UAM formatted output data for all variables

            WRITE(ODEV,ERR=7011) NSEG, (SPCNAM( K,J ),J=1,10),&
            &           (( EMIS( I,J ),I = 1,NCOLS3D ), J=1,NROWS3D )

        END DO       ! Output variables loop

!.............  Update times and dates for next pass through loop

        SDATE = SDATEP1
        STIME = STIMEP1
        CALL NEXTIME( SDATEP1, STIMEP1, TSTEP)

    END DO ! Time step loop


!.........  End of program

    CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

!.........  Error message for reaching the end of file too soon
999 MESG = 'End of file reached unexpectedly. ' //&
    &       'Check format of variable ' // CRLF() // BLANK5 //&
    &       'mapping file.'
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

!******************    ERROR MESSAGES     ******************************

7010 WRITE (*,*) 'ERROR:  WRITING UAM HEADER OF OUTPUT'
    STOP
7011 WRITE (*,*) 'ERROR:  WRITING UAM RECORD OUTPUT'
    STOP

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT ( A )

93100 FORMAT ( A1 )

93200 FORMAT (/, 5X, 'Writing out for date ',i7,' and time ',i6,'.')

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10 ( A, :, I5, :, 2X ) )


END PROGRAM SMK2EMIS


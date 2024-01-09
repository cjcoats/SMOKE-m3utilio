
SUBROUTINE MXGRPEMIS( NINVGRP, TSTEP, SDATE, STIME, NSTEPS )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !     This routine determines the maximum daily emissions during the
    !     period being modeled.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Written 7/2001 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",
    !       fix bug in SORTR2 call, and related changes
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
    !***********************************************************************
    USE M3UTILIO

    !.......   MODULES for public variables
    !.......   This module is the source inventory arrays
    USE MODSOURC, ONLY: CIFIP, TZONES

    !.......  This module contains arrays for plume-in-grid and major sources
    USE MODELEV, ONLY: LELVRNK, LPNGRNK, GROUPID, NEVPEMV, MXEMIS,  &
                       MXEIDX, MXRANK, EVPEMIDX, GINDEX, SRCXL,     &
                       SRCYL

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, CRL, NIPOL, NSRC, EINAM

    !.......  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: NCOLS, NROWS
    USE MODGRDLIB

    !.......This module is required by the FileSetAPI
    USE MODFILESET

    IMPLICIT NONE

    !.......   INCLUDES:
    INCLUDE 'EMCNST3.h90'       ! emissions constant parameters
    INCLUDE 'SETDECL.h90'       ! FileSetAPI variables and functions

    INCLUDE 'CONST3.EXT'        ! physical and mathematical constants from I/O API

    !.......   ARGUMENTS and their descriptions:
    INTEGER     , INTENT (IN) :: NINVGRP      ! no. of inventory groups
    INTEGER     , INTENT (IN) :: TSTEP        ! time step
    INTEGER  , INTENT(IN OUT) :: SDATE        ! Julian start date limit
    INTEGER  , INTENT(IN OUT) :: STIME        ! start time limit
    INTEGER  , INTENT(IN OUT) :: NSTEPS       ! number of time steps limit

    !.......   EXTERNAL FUNCTIONS and their descriptions:
    INTEGER, EXTERNAL :: GETFLINE

    !.......   Local allocatable arrays
    INTEGER :: DAYBEGT( NSRC )      ! day beginning time by source
    INTEGER :: DAYENDT( NSRC )      ! day ending time by source
    LOGICAL :: LDAYSAV( NSRC )       ! true: daylight savings used by source
    REAL    :: RGROUP ( NSRC )       ! tmp emissions by source
    REAL    :: EMIS   ( NSRC )       ! tmp emissions by source
    REAL    :: SRCSUM ( NSRC,NEVPEMV )     ! emissions summing by source
    REAL    :: MXGRPEM( NINVGRP,NEVPEMV )     ! emissions max for groups
    REAL    :: GRPSUM ( NINVGRP,NEVPEMV )     ! emissions summing by group

    INTEGER, ALLOCATABLE :: MXRECLST( : )      ! maximum hours per file in list
    INTEGER, ALLOCATABLE :: SDATELST( : )      ! start dates for files in list
    INTEGER, ALLOCATABLE :: STIMELST( : )      ! start times for files in list
    INTEGER, ALLOCATABLE :: TIMIDX  ( : )      ! index for sorting dates/times

    CHARACTER(512), ALLOCATABLE :: PTMPLIST( : )     ! list of PTMP files

    !.......   Logical file names and unit numbers
    INTEGER         TDEV        !  list file unit number
    CHARACTER(16)   TNAME       !  logical name for hourly inventory input file

    !.......   Local fixed arrays
    !.......   OTHER LOCAL VARIABLES and their descriptions:
    INTEGER     G, J, K, K2, L, L2, N, S, V, T        ! indices and counters

    INTEGER         COL               ! tmp grid column
    INTEGER         DAY               ! tmp day of week
    INTEGER         ED                ! tmp end date
    INTEGER         EDATE             ! ending date YYYYDDD
    INTEGER         ET                ! tmp end time
    INTEGER         ETIME             ! ending time HHMMSS
    INTEGER         FPTR              ! file index
    INTEGER         HDIF              ! date/time difference in hours
    INTEGER         IOS               ! i/o status
    INTEGER         ISECS             ! tmp number seconds
    INTEGER         JDATE             ! Julian date
    INTEGER         JTIME             ! time HHMMSS
    INTEGER         LDATE             ! date from previous iteration
    INTEGER         NFILES            ! number of list files
    INTEGER         NPTLINES          ! number of lines in list file TDEV
    INTEGER         NS                ! tmp number time steps
    INTEGER         PG                ! previous G
    INTEGER         ROW               ! tmp grid row
    INTEGER         SD                ! start date
    INTEGER         ST                ! start time

    LOGICAL :: EFLAG    = .FALSE.     ! true: error detected
    LOGICAL :: FILEOPEN = .FALSE.     ! true: an I/O API file is open

    CHARACTER(16)   INLGCNAM          ! input logical file name
    CHARACTER(512)  PTMPFILE          ! tmp physical file name
    CHARACTER(512)  PREVFILE          ! previous physical file name
    CHARACTER(256)  MESG
    CHARACTER(NAMLEN3) CBUF       ! tmp pollutant name

    CHARACTER(16), PARAMETER :: PROGNAME = 'MXGRPEMIS'       !  program name

    !***********************************************************************
    !   begin body of subroutine MXGRPEMIS

    !.......  Create note about why temporal file is being read in
    MESG = 'NOTE: Hourly emissions file is required to determine '//    &
           'actual daily emissions ' // CRLF() // BLANK10 //            &
           'for evaluating emissions-based criteria for elevated '//    &
           'and/or PinG ' // CRLF() // BLANK10 // 'source selection'
    CALL M3MSG2( MESG )

    !.......  Open hourly emissions file
    MESG = 'Enter logical name for POINT'//    &
           ' HOURLY EMISSIONS FILES LIST (or "NONE" for PTMP)'

    INLGCNAM = CRL//'TMPLIST'
    TDEV     = PROMPTFFILE( MESG, .TRUE., .TRUE., INLGCNAM, PROGNAME )
    TNAME    = 'PTMP'

    !.......  If we don't have a list, get physical file for PTMP to use
    IF( TDEV .LE. 0 ) THEN
        CALL ENVSTR( TNAME, 'Hourly point source',  ' ', PTMPFILE, IOS )
        IF( IOS .NE. 0 ) THEN
            MESG = 'No setting for hourly point emissions file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        NPTLINES = 1      ! fake setting as if 1 file listed in list file

    !.......  Otherwise, get the number of lines from the list file
    ELSE
        NPTLINES = GETFLINE( TDEV, INLGCNAM // ' file' )

    END IF

    !.......  Allocate arrays for maximum emission daily totals per source
    !           in local time zone for all pollutants used as a selection criteria
    !.......  Allocate global arrays for maximum emis by stack group in
    !           local time zone for all pollutants used as a selection criteria
    !.......   Allocate and initialize arrays for storing hourly input file info
    ALLOCATE( MXEMIS( NSRC,NEVPEMV ),   &
              MXEIDX( NSRC,NEVPEMV ),   &
              MXRANK( NSRC,NEVPEMV ),   &
            PTMPLIST( NPTLINES ),       &
            SDATELST( NPTLINES ),       &
            STIMELST( NPTLINES ),       &
            MXRECLST( NPTLINES ),       &
              TIMIDX( NPTLINES ), STAT=IOS )
    CALL CHECKMEM( IOS, 'MXEMIS...TIMIDX', PROGNAME )

    RGROUP  = FLOAT( GROUPID )
    MXEMIS  = 0.            ! array
    MXEIDX  = 0             ! array
    SRCSUM  = 0.            ! array
    DAYBEGT = 0             ! array
    DAYENDT = 0             ! array
    MXGRPEM = 0.            ! array
    GRPSUM  = 0.            ! array
    LDAYSAV = .FALSE.       ! array

    !.......  Initialize sorting indices
    DO S = 1, NSRC
        MXEIDX( S,1:NEVPEMV ) = S       ! array
    END DO

    PTMPLIST = ' '       ! array
    SDATELST = 0         ! array
    STIMELST = 0         ! array
    MXRECLST = 0         ! array
    TIMIDX   = 0         ! array

    !.......   Read lines from hourly emissions list file
    IF( TDEV .GT. 0 ) THEN
        CALL RDLINES( TDEV, INLGCNAM//' file', NPTLINES, PTMPLIST )
        CLOSE( TDEV )
    ELSE
        PTMPLIST( 1 ) = PTMPFILE
    END IF

    MESG = 'Checking hourly emissions files...'
    CALL M3MSG2( MESG )

    K = 0

    !.......  Loop through met files and check time period
    DO N = 1, NPTLINES

        !.......  Close previous file if needed
        IF( FILEOPEN ) THEN
            IF( .NOT. CLOSESET( TNAME ) ) THEN
                MESG = 'Could not close hourly emissions file ' // TRIM( PTMPFILE )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ELSE
                FILEOPEN = .FALSE.
            END IF
        END IF

        !.......  Get physical file name for current iteration
        PTMPFILE = PTMPLIST( N )

        !.......  Skip any blank lines
        IF( PTMPFILE .EQ. ' ' ) CYCLE

        !.......  Set logical file name
        IF( .NOT. SETENVVAR( TNAME, PTMPFILE ) ) THEN
            EFLAG = .TRUE.
            MESG = 'INTERNAL ERROR: Could not set logical name for file ' // TRIM( PTMPFILE )
            CALL M3MSG2( MESG )
            CYCLE
        END IF

        !.......  Try to open file
        IF( .NOT. OPENSET( TNAME, FSREAD3, PROGNAME ) ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: Could not open hourly emissions file ' //TRIM( PTMPFILE )
            CALL M3MSG2( MESG )
            CYCLE
        ELSE
            FILEOPEN = .TRUE.
        END IF

        !.......  Read description of file
        IF( .NOT. DESCSET( TNAME,ALLFILES ) ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: Could not get description of file ' // TRIM( PTMPFILE )
            CALL M3MSG2( MESG )
            CYCLE
        END IF

        !.......  Store start date/time and end date/time
        K = K + 1
        TIMIDX  ( K ) = K
        PTMPLIST( K ) = PTMPLIST( N )          ! remove blanks from list
        SDATELST( K ) = SDATE3D
        STIMELST( K ) = STIME3D
        MXRECLST( K ) = MXREC3D

        !.......  Check the number of sources
        CALL CHKSRCNO( CATEGORY, TNAME, NROWS3D, NSRC, EFLAG )

        !.......  Check the variable number and names
        IF( NVARSET .NE. NIPOL ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: The number of ' //          &
                  'pollutants in the PTMP file (', NVARSET,         &
                  CRLF() // BLANK10 // 'is inconsitent with '//     &
                  'the number in the PNTS file.'
            CALL M3MSG2( MESG )

        !.......  If number is okay, check names
        ELSE

            DO V = 1, NIPOL
                IF( EINAM( V ) .NE. VNAMESET( V ) ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Pollutant "'//          &
                          TRIM( EINAM( V ) )// '" is in a different '// &
                          'order in PNTS than in PTMP.'
                    CALL M3MSG2( MESG )

                END IF
            END DO          ! Loop on pollutants

        END IF              ! If number of pollutants the same

    END DO         ! End loop over list of files

    NFILES = K

    !.......  Sort files by date and time
    CALL SORTI2( NFILES, TIMIDX, SDATELST, STIMELST )

    !.......  Determine which files will be used for each day and check that
    !           there is coverage for all days in period
    DO N = 1, NFILES - 1

        K  = TIMIDX( N )            ! pull out files in sorted order
        K2 = TIMIDX( N+1 )

        !....... Compute end time for file on current iteration
        EDATE = SDATELST( K )
        ETIME = STIMELST( K )
        CALL NEXTIME( EDATE, ETIME, ( MXRECLST(K)-1 ) * TSTEP3D )

        !.......  Compare end date/time to next file start date/time...
        HDIF = SECSDIFF( EDATE, ETIME, SDATELST(K2), STIMELST(K2) )
        HDIF = HDIF / 3600          ! integer math

        !.......  If end date/time >= start date/time of next file, trucate this
        !               file info to not duplicate hours on next day.
        IF( HDIF .LE. 0 ) THEN
            MXRECLST( K ) = MXRECLST( K ) - HDIF - 1

        !.......  If end date/time < start date/time of next file, give error
        !               and continue.
        ELSE IF ( HDIF .GT. 1 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: List of hourly emissions files have a', HDIF,   &
                   CRLF() // BLANK10 // 'hour gap between files:' //                    &
                   CRLF() // BLANK10 // TRIM( PTMPLIST( K ) )// ' and' //               &
                    CRLF()// BLANK10 // TRIM( PTMPLIST( K2 ) )
            CALL M3MSG2( MESG )

        END IF

    END DO

    !.......  Stop if errors found
    IF ( EFLAG ) THEN

        MESG = 'Problem with hourly emissions input files.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    !.......  If no errors, then set full duration of input files and compare
    !           to settings of the episode duration
    ELSE

        !.......  Initialize end time based on subroutine arguments
        EDATE = SDATE
        ETIME = STIME
        CALL NEXTIME( EDATE, ETIME, NSTEPS*TSTEP )

        !.......  Initialize input file start time by first file after sorting
        SDATE3D = SDATELST( TIMIDX( 1 ) )
        STIME3D = STIMELST( TIMIDX( 1 ) )
        MXREC3D = SUM( MXRECLST )               ! array

        !.......  Ensure time information is consistent with arguments
        ISECS = SECSDIFF( SDATE, STIME, SDATE3D, STIME3D )

        IF( ISECS .GT. 0 ) THEN          ! SDATE3D/STIME3D are later
            SDATE = SDATE3D
            STIME = STIME3D
        END IF

        ED = SDATE3D
        ET = STIME3D
        CALL NEXTIME( ED, ET, ( MXREC3D-1 ) * TSTEP3D )

        ISECS = SECSDIFF( EDATE, ETIME, ED, ET )

        IF( ISECS .LT. 0 ) THEN          ! ED/ET are earlier
            EDATE = ED
            ETIME = ET
        END IF

        NS = SECSDIFF( SDATE, STIME, EDATE, ETIME )/ 3600

        IF( NS .LE. 0 ) THEN
            MESG = 'Because of file ' // TNAME //    &
                   ', dates and times do not overlap at all'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSE IF ( NS .NE. NSTEPS ) THEN
            WRITE( MESG,94010 ) 'WARNING: The contents of the ' //      &
                   'PTMP file overlap with the requested '//            &
                   CRLF()// BLANK10 // 'date and time for only ', NS, 'hours.'
            CALL M3MSG2( MESG )

        END IF

        NSTEPS = NS

    END IF                  ! If header was read okays

    !.......  If errors found so far, then abort
    IF ( EFLAG ) THEN
        MESG = 'Problem with hourly emissions files.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    !.......  Otherwise, output status message
    ELSE
        MESG = 'Computing source and group maximum daily emissions totals...'
        CALL M3MSG2( MESG )
    END IF

    !.......  Loop through hours
    JDATE = SDATE
    JTIME = STIME
    LDATE = -9
    PREVFILE = ' '
    DO T = 1, NSTEPS

        !.......  When new day...
        IF ( JDATE .NE. LDATE ) THEN

        !.......  Write message for day of week and date
            DAY = WKDAY( JDATE )
            MESG = 'Processing ' // DAYS( DAY ) // ' ' // MMDDYY( JDATE )
            CALL M3MSG2( MESG )

        !.......  Create array of which sources are affected by daylight
        !                  savings
            CALL GETDYSAV( NSRC, CIFIP, LDAYSAV )

        !.......  Set start and end hours of day for all sources
            CALL SETSRCDY( NSRC, JDATE, TZONES, LDAYSAV, .FALSE., DAYBEGT, DAYENDT )
        END IF

        !.......  Determine input file for this hour
        DO N = 1, NFILES
            K = TIMIDX( N )
            SD = SDATELST( K )
            ST = STIMELST( K )
            ED = SD
            ET = ST
            CALL NEXTIME( ED, ET, ( MXRECLST(K)-1 ) * TSTEP3D )

            IF( JDATE .GE. SD .AND. JDATE .LE. ED .AND.    &
                JTIME .GE. ST .AND. JTIME .LE. ET       ) THEN
                FPTR = K
                EXIT           ! exit loop to set file to use
            END IF

        END DO

        !.......  Get physical file name for current iteration
        PTMPFILE = PTMPLIST( FPTR )

        !.......  Close previous file if needed
        IF( PTMPFILE .NE. PREVFILE ) THEN
            IF( FILEOPEN ) THEN
                IF( .NOT. CLOSESET( TNAME ) ) THEN
                    MESG = 'Could not close hourly emissions file ' // PREVFILE
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                ELSE
                    FILEOPEN = .FALSE.
                END IF
            END IF

            PREVFILE = PTMPFILE

        END IF

        !.......  Set logical file name
        IF( .NOT. SETENVVAR( TNAME, PTMPFILE ) ) THEN
            EFLAG = .TRUE.
            MESG = 'INTERNAL ERROR: Could not set logical name for file ' // PTMPFILE
            CALL M3MSG2( MESG )
            CYCLE
        END IF

        !.......  Try to open file
        IF( .NOT. OPENSET( TNAME, FSREAD3, PROGNAME ) ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: Could not open hourly emissions file ' // PTMPFILE
            CALL M3MSG2( MESG )
            CYCLE
        ELSE
            FILEOPEN = .TRUE.
        END IF

        !.......  Loop through pollutants that are used as selection criteria
        DO K = 1, NEVPEMV

            !.......  Set global emissions variable index
            V    = EVPEMIDX( K )
            CBUF = EINAM( V )

            !.......  Read emissions value
            IF ( .NOT. READSET( TNAME, CBUF, 1, ALLFILES, JDATE, JTIME, EMIS ) ) THEN
                MESG = 'Could not read ' // TRIM( CBUF ) // ' from ' // TNAME
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

            END IF

            !.......  Loop through sources and add up daily total based on
            !                   local day
            PG   = -9
            DO J = 1, NSRC

                S = GINDEX ( J )
                G = GROUPID( S )

                !.......  Initialize group total and source total if we are
                !                       at the first hour of the day for this source
                IF ( DAYBEGT( S ) .EQ. JTIME ) THEN

                    !...........  Only initialize group total if the group is new
                    IF ( G .GT. 0 .AND. G .NE. PG ) THEN
                        GRPSUM( G,: ) = 0.
                    END IF

                    SRCSUM( S,: ) = 0.

                END IF

                !.......  If source is in a group, add emissions
                !.......  Also store  group number for next iteration
                IF ( G .GT. 0 ) THEN
                    GRPSUM (G,K) = GRPSUM( G,K ) + EMIS( S )
                    MXGRPEM(G,K) = MAX( MXGRPEM(G,K), GRPSUM(G,K) )
                    PG = G

                !.......  Otherwise, add emissions for this time step for source
                ELSE
                    SRCSUM( S,K ) = SRCSUM( S,K ) + EMIS( S )
                    MXEMIS( S,K ) = MAX( MXEMIS(S,K), SRCSUM(S,K) )

                END IF


            END DO       ! End loop on sources

        END DO          ! End loop on pollutants used as selection criteria

        LDATE = JDATE
        CALL NEXTIME( JDATE, JTIME, TSTEP )

    END DO           ! End loop on time steps

    !.......  Replace source emissions with group emissions and exclude
    !         sources that are outside of the grid...
    !.......  Loop through pollutants that are used as selection criteria
    DO K = 1, NEVPEMV

        DO S = 1, NSRC

            IF( .NOT. INGRID( SRCXL( S ), SRCYL( S ), NCOLS, NROWS, COL, ROW ) ) THEN
                MXEMIS( S,K ) = AMISS3
                CYCLE
            END IF

            G = GROUPID( S )
            IF ( G .LE. 0 ) CYCLE      ! Skip sources without group

            MXEMIS( S,K ) = MXGRPEM( G,K )

        END DO      ! end loop on sources

    END DO          ! End loop on pollutants used as selection criteria

    !.......  Loop through pollutants that are used as selection criteria
    !           and sort max daily emissions in ascending order.
    IF ( LELVRNK .OR. LPNGRNK ) THEN

        !.......  Write status message
        MESG = 'Ranking sources and groups by emissions values...'
        CALL M3MSG2( MESG )

        DO K = 1, NEVPEMV

            !.......  Note that all sources in a group have the same group-total
            !                   emissions value at this point.
            CALL SORTR2( NSRC, MXEIDX(1,K), MXEMIS(1,K), RGROUP )

            !.......  Reset indices if a source is in a group, invert the
            !         order, and change the index to a ranking.
            N = 0
            PG = -9
            DO J = NSRC, 1, -1

                S = MXEIDX ( J,K )
                G = GROUPID( S )

                !.......  If source is not a group or group has changed
                IF ( G .LE. 0 .OR. G .NE. PG ) N = N + 1

                !.......  Reset ranked number for source such that all
                !                       members of a group will have the same rank
                MXRANK( S,K ) = N

                PG = G

            END DO

        END DO

    END IF

    !.......  Deallocate local memory
    DEALLOCATE( MXRECLST, SDATELST, STIMELST, TIMIDX, PTMPLIST, MXEIDX )

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10 ( A, :, I8, :, 2X  ) )

END SUBROUTINE MXGRPEMIS

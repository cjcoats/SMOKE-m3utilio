
SUBROUTINE OPENTMPIN( UFLAG, PFLAG, ENAME, ANAME, DNAME, HNAME,    &
                      SDEV, CDEV, HDEV, KDEV, PYEAR )

    !***********************************************************************
    !  subroutine body starts at line 123
    !
    !  DESCRIPTION:
    !      This subroutine opens the file or files for output from the tmppoint
    !      program. It also populates the MODINFO module with the source-category
    !      specific information. It also determines the name of the temperature
    !      variable.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !      Created 7/99 by M. Houyoux
    !
    !     Revised ??/???? by ??
    !
    !     Version 07/2014 by C.Coats for  new GENTPRO CSV profiles and cross-references
    !     These are purely local to SUBROUTINE PROCTPRO()
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !**************************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !        System
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
    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, CRL, NSRC, NIACT, INVPIDX, EANAM,    &
                       EINAM, NIPOL, NIPPA, NPPOL, NPACT, ACTVTY

    IMPLICIT NONE

    !.......   INCLUDES

    INCLUDE 'EMCNST3.h90'       ! emissions constat parameters
    INCLUDE 'SETDECL.h90'       !  FileSetAPI variables and functions

    !.......   SUBROUTINE ARGUMENTS

    LOGICAL     , INTENT    (IN) :: UFLAG     ! use uniform temporal profile
    LOGICAL     , INTENT(IN OUT) :: PFLAG     ! use episode time periods
    CHARACTER(*), INTENT(IN OUT) :: ENAME     ! name for I/O API inven input
    CHARACTER(*), INTENT(IN OUT) :: ANAME     ! name for ASCII inven input
    CHARACTER(*), INTENT   (OUT) :: DNAME     ! day-spec file
    CHARACTER(*), INTENT   (OUT) :: HNAME     ! hour-spec file
    INTEGER     , INTENT   (OUT) :: SDEV      ! unit no.: ASCII inven file
    INTEGER     , INTENT   (OUT) :: CDEV      ! unit no.: region codes file
    INTEGER     , INTENT   (OUT) :: HDEV      ! unit no.: holidays file
    INTEGER     , INTENT   (OUT) :: KDEV      ! unit no.: time periods file
    INTEGER     , INTENT   (OUT) :: PYEAR     ! projected year

    !.......   EXTERNAL FUNCTIONS and their descriptions:

    CHARACTER(NAMLEN3), EXTERNAL :: GETCFDSC
    INTEGER           , EXTERNAL :: GETIFDSC
    LOGICAL           , EXTERNAL :: USEEXPGEO

    !.......   Other local variables

    INTEGER         IDEV            ! tmp unit number if ENAME is map file
    INTEGER         IOS             ! status from environment variables
    INTEGER         I,J               ! index
    INTEGER         L               ! string length

    LOGICAL      :: DFLAG = .FALSE.          ! true: day-specific  file available
    LOGICAL      :: EFLAG = .FALSE.          ! true: error found
    LOGICAL      :: HFLAG = .FALSE.          ! true: hour-specific file available
    LOGICAL         OFLAG           ! true: average day emissions needed
    LOGICAL         XFLAG           ! true: use daylight time exemptions file

    CHARACTER(16)   INAME           ! tmp name for inven file of unknown fmt
    CHARACTER(256)  MESG            ! message buffer

    CHARACTER(NAMLEN3)  NAMBUF     ! file name buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'OPENTMPIN'     ! program name

    !***********************************************************************
    !   begin body of subroutine OPENTMPIN

    !.......  Get environment variables that control program behavior
    DFLAG = ENVYN( 'DAY_SPECIFIC_YN', 'Use day-specific data', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "DAY_SPECIFIC_YN"', 2 )
    END IF

    HFLAG = ENVYN( 'HOUR_SPECIFIC_YN', 'Use hour-specific data', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "HOUR_SPECIFIC_YN"', 2 )
    END IF


    !.......  Waring message for imcompatiblity of FILL_ANNUAL setting from Smkinven

    !.......  Set average day emissions flag (INVPIDX)
    OFLAG = ENVYN( 'SMK_AVEDAY_YN', MESG, .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK_AVEDAY_YN"', 2 )
    END IF
    IF( OFLAG ) THEN
        MESG = 'CRITICAL: If FILL_ANNUAL was set to Y in SMKINVEN run, ' // &
               CRLF() // BLANK10 // 'SMK_AVEDAY_YN must be set to N in TEMPORAL run.'
        CALL M3MSG2( MESG )
    END IF

    !.......  Prompt for and open inventory file
    INAME = ENAME
    MESG = 'Enter logical name for the MAP INVENTORY file'
    IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INAME, PROGNAME )

    !.......  Open and read map file
    CALL RDINVMAP( INAME, IDEV, ENAME, ANAME, SDEV )

    IF( DFLAG ) THEN
        NAMBUF = PROMPTMFILE(    &
                 'Enter logical name for DAY-SPECIFIC file',    &
                 FSREAD3, CRL // 'DAY', PROGNAME )
        DNAME = NAMBUF
    ELSE
        DNAME = 'NONE'
    END IF

    IF( HFLAG ) THEN
        NAMBUF = PROMPTMFILE(    &
                 'Enter logical name for HOUR-SPECIFIC file',    &
                 FSREAD3, CRL // 'HOUR', PROGNAME )
        HNAME = NAMBUF
    ELSE
        HNAME = 'NONE'
    END IF

    !.......  Open the time periods that Temporal should process

    MESG = 'Enter logical name for Episode Time Periods file inputs list (or "NONE")'
    KDEV = PROMPTFFILE( MESG, .TRUE., .TRUE.,'PROCDATES', PROGNAME)
    IF( KDEV == -2 ) THEN
        PFLAG = .FALSE.
    ELSE
        PFLAG = .TRUE.
    END IF

    !.......  Store source-category-specific header information,
    !        including the inventory pollutants in the file (if any).  Note that
    !        the I/O API head info is passed by include file and the
    !        results are stored in module MODINFO.
    !.......  Set average day emissions flag (INVPIDX)
    IF( OFLAG ) INVPIDX = 1
    CALL GETSINFO( ENAME )

    !.......  Reset activity to pollutant to create hourly VMT without running EMISFAC
    !        To support MOVES-SMOKE integratoin tool approach.
    !        Only no-SPEED activity data (i.e. VMT) will be treated like pollutant for later Temporal program
    IF( NIACT > 0 ) THEN
        IF( ALLOCATED( EINAM ) ) DEALLOCATE( EINAM )
        ALLOCATE( EINAM( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EINAM', PROGNAME )

        NIPOL = NIACT
        NPPOL = 2
        NIACT = 0
        NPACT = 0
        EINAM = EANAM
        ACTVTY = ' '
    END IF

    PYEAR = GETIFDSC( FDESC3D, '/PROJECTED YEAR/', .FALSE. )

    !.......  Store non-category-specific header information
    NSRC = NROWS3D

    !.......  Open region codes file for determining daylight savings time status
    IF( .NOT. USEEXPGEO() ) THEN
        CDEV = PROMPTFFILE(    &
               'Enter logical name for COUNTRY, STATE, AND ' //    &
               'COUNTY file', .TRUE., .TRUE., 'COSTCY', PROGNAME )
    END IF

    !.......  Open holidays file for determining holidays by region
    HDEV = PROMPTFFILE(    &
               'Enter logical name for HOLIDAYS file',    &
               .TRUE., .TRUE., 'HOLIDAYS', PROGNAME )

    !.......  Abort if error was found
    IF ( EFLAG ) THEN

        MESG = 'Problem with input files'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94000 FORMAT( I2.2 )

94010 FORMAT( 10( A, :, I8, :, 1X ) )

    !******************  INTERNAL SUBPROGRAMS   ******************************

CONTAINS

    !.......  This internal subprogram tries to retrieve the I/O API header
    !        and aborts if it was not successful
    SUBROUTINE RETRIEVE_IOAPI_HEADER( FILNAM )

        !.......  Subprogram arguments
        CHARACTER(*) FILNAM

        !----------------------------------------------------------------------

        IF ( .NOT. DESC3( FILNAM ) ) THEN

            MESG = 'Could not get description of file "' // TRIM( FILNAM ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

    END SUBROUTINE RETRIEVE_IOAPI_HEADER

END SUBROUTINE OPENTMPIN


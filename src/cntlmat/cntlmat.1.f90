
PROGRAM CNTLMAT

!***********************************************************************
!  program body starts at line 136
!
!  DESCRIPTION:
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Copied from CTLPMAT.F version 4.4 by M. Houyoux 3/99
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
!**********************************************************************

    USE M3UTILIO
!.........  MODULES for public variables
!.........  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, CRL, NSRC, BYEAR

!.........This module is required by the FileSetAPI
    USE MODFILESET

    IMPLICIT NONE

!...........   INCLUDES:
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
    INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

!...........   EXTERNAL FUNCTIONS and their descriptions:
    INTEGER, EXTERNAL :: GETIFDSC

!...........  LOCAL PARAMETERS and their descriptions:

    CHARACTER(50), PARAMETER ::&
    &CVSW = '$Name SMOKEv5.0_Jun2023$' ! CVS release tag

!...........   LOCAL VARIABLES and their descriptions:

!.........  Array that contains the names of the inventory variables needed for
!           this program
    CHARACTER(NAMLEN3) IVARNAMS( MXINVARR )

!...........   Arrays for each of the packets
    INTEGER      PKTCNT( NPACKET )  ! count of records per packet
    INTEGER      PKTBEG( NPACKET )  ! first line of packet in input file
    INTEGER      XRFCNT( NPACKET )  ! cross-reference info memory needs

!...........   Allocatable local arrays

    LOGICAL, ALLOCATABLE :: TFLAG( : )  !  flags:  track these sources

!...........   Logical names and unit numbers
    INTEGER         CDEV         !  control file unit no.
    INTEGER         CTMPDEV      !  file unit no. for tmp CTL file
    INTEGER         GTMPDEV      !  file unit no. for tmp CTG file
    INTEGER         IDEV         !  tmp unit number if inven is map-formatted
    INTEGER         LDEV         !  log file unit no.
    INTEGER         LTMPDEV      !  file unit no. for tmp ALW file
    INTEGER         MTMPDEV      !  file unit no. for tmp MACT file
    INTEGER         PTMPDEV      !  file unit no. for tmp PROJ file
    INTEGER         SDEV         !  ASCII part of inventory unit no.
    INTEGER      :: WDEV = 0     !  warnings/error unit no.

    CHARACTER(16)   ANAME   ! logical name for ASCII inventory input file
    CHARACTER(16)   ENAME   ! logical name for i/o api inventory input file
    CHARACTER(16)   INAME   ! tmp logical name for inven file of unknown fmt
    CHARACTER(16)   MNAME   ! logical name for multiplicative control matrix
    CHARACTER(16)   PNAME   ! logical name for projection matrix
    CHARACTER(16)   SNAME   ! logical name for ascii inventory input file

!...........   Other local variables
    INTEGER      :: CPYEAR = -1  !  control packet year to project to
    INTEGER         IOS          !  I/O status
    INTEGER         ENLEN        !  length of the emissions inven name
    INTEGER         NCPE         !  no control packet entries
    INTEGER         NINVARR      !  number inventory variables to input
    INTEGER      :: PYEAR   = 0  !  projected year of inventory
    INTEGER         SYEAR        !  year for projecting from

    LOGICAL      :: CFLAG   = .FALSE.  ! true: control cntls in use
    LOGICAL      :: EFLAG   = .FALSE.  ! true: error has occurred
    LOGICAL      :: GFLAG   = .FALSE.  ! true: ctg cntls in use
    LOGICAL      :: JFLAG   = .FALSE.  ! true: projections in use
    LOGICAL      :: KFLAG   = .FALSE.  ! true: tracking file in use
    LOGICAL      :: LFLAG   = .FALSE.  ! true: allowable cntls in use
    LOGICAL      :: LCTMP   = .FALSE.  ! true: tmp control file written
    LOGICAL      :: LPTMP   = .FALSE.  ! true: tmp projection file written
    LOGICAL      :: MFLAG   = .FALSE.  ! true: mact cntls is use
    LOGICAL      :: RFLAG   = .FALSE.  ! true: reactivty cntls in use
    LOGICAL      :: YFLAG   = .FALSE.  ! true: projection entries have years

    CHARACTER(7)    ACTION             ! buffer for PKTLOOP action
    CHARACTER(300)  MESG               ! message buffer

    CHARACTER(16) :: PROGNAME = 'CNTLMAT'   !  program name

!***********************************************************************
!   begin body of program CNTLMAT

    LDEV = INIT3()

!.........  Write out copyright, version, web address, header info, and prompt
!           to continue running the program.
    CALL INITEM( LDEV, CVSW, PROGNAME )

!.........  Get environment variable values...
    KFLAG = ENVYN( 'CONTROL_TRACKING',&
    &               'Use a special file to track specific sources',&
    &               .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "CONTROL_TRACKING"', 2 )
    END IF

!.........  Warning if tracking file use is attempted
    IF( KFLAG ) THEN
        MESG = 'WARNING: Specific source tracking has not been ' //&
        &       'implemented. ' // CRLF() // BLANK10 //&
        &       'CONTROL_TRACKING variable will have no effect.'
        CALL M3MSG2( MESG )
    END IF

!.........  Set source category based on environment variable setting
    CALL GETCTGRY

!.........  Get inventory file names given source category
    CALL GETINAME( CATEGORY, ENAME, ANAME )

!.........  Prompt for and open inventory file
    INAME = ENAME
    MESG = 'Enter logical name for the MAP INVENTORY file'
    IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INAME, PROGNAME )

!.........  Open and read map file
    CALL RDINVMAP( INAME, IDEV, ENAME, ANAME, SDEV )

    CDEV = PROMPTFFILE(&
    &         'Enter logical name for ASCII CONTROL PACKETS file',&
    &         .TRUE., .TRUE., 'GCNTL', PROGNAME )

    WDEV = PROMPTFFILE(&
    &         'Enter logical name for output WARNINGS/ERRORS file',&
    &         .FALSE., .TRUE., CRL//'CTLWARN', PROGNAME )

!.........  Store source-category-specific header information,
!           including the inventory pollutants in the file (if any).  Note that
!           the I/O API head info is passed by include file and the
!           results are stored in module MODINFO.
    CALL GETSINFO( ENAME )

!.........  Check for future-year inventory file
    PYEAR = GETIFDSC( FDESC3D, '/PROJECTED YEAR/', .FALSE. )

!.........  Set starting year on which to base possible projections created
!           in this program
    SYEAR = BYEAR
    IF( PYEAR .GT. 0 ) SYEAR = PYEAR

!.........  Set inventory variables to read for all source categories
    IVARNAMS( 1 ) = 'INVYR'
    IVARNAMS( 2 ) = 'CSCC'
    IVARNAMS( 3 ) = 'CSOURC'

!.........  Set inventory variables to read for specific source categories
    IF( CATEGORY .EQ. 'AREA' ) THEN
        NINVARR = 7
        IVARNAMS( 4 ) = 'CISIC'
        IVARNAMS( 5 ) = 'CMACT'
        IVARNAMS( 6 ) = 'CSRCTYP'
        IVARNAMS( 7 ) = 'CIFIP'

    ELSE IF( CATEGORY .EQ. 'MOBILE' ) THEN
        NINVARR = 5
        IVARNAMS( 4 ) = 'CVTYPE'
        IVARNAMS( 5 ) = 'CIFIP'

    ELSE IF( CATEGORY .EQ. 'POINT' ) THEN
        NINVARR = 7
        IVARNAMS( 4 ) = 'CISIC'
        IVARNAMS( 5 ) = 'CMACT'
        IVARNAMS( 6 ) = 'CSRCTYP'
        IVARNAMS( 7 ) = 'CIFIP'

    END IF

!.........  Allocate memory for and read required inventory characteristics
    CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

!.........  Build unique lists of SCCs per SIC from the inventory arrays
    CALL GENUSLST

!.........  Initialize arrays
    PKTCNT = 0  ! array
    PKTBEG = 0  ! array
    XRFCNT = 0  ! array

!.........  Allocate memory for control packet information in input file.
    CALL ALOCPKTS( CDEV, WDEV, SYEAR, CPYEAR, PKTCNT,&
    &               PKTBEG, XRFCNT, LPTMP, LCTMP )

!.........  Set the flags that indicate which packets are valid
    GFLAG = ( PKTCNT( 1 ) .GT. 0 )
    CFLAG = ( PKTCNT( 2 ) .GT. 0 )
    LFLAG = ( PKTCNT( 3 ) .GT. 0 )
    RFLAG = ( PKTCNT( 5 ) .GT. 0 )
    JFLAG = ( PKTCNT( 6 ) .GT. 0 )
    MFLAG = ( PKTCNT( 7 ) .GT. 0 )

!.........  Cannot have projection packet and invalid projection year
    IF( JFLAG .AND. CPYEAR .LT. 1900 ) THEN
        MESG = 'Misformatted projection packet header is returning '&
        &       // 'year < 1900'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Process packets: this means read packet, sort it, group it into
!           grouped, x-ref structures in MODXREF, and assign to sources, and
!           write outputs for non-pollutant specific packets.
!.........  For control matrices that depend on pollutants, temporary files
!           will be written if there is more than one pollutant group, and
!           these will be used to store the control data index information for
!           each packet type while determining the pollutants to use in opening
!           the final output files.

    ACTION = 'PROCESS'
    CALL PKTLOOP( CDEV, PTMPDEV, CTMPDEV, GTMPDEV, LTMPDEV, MTMPDEV,&
    &              WDEV, CPYEAR, ACTION, ENAME, PKTCNT, PKTBEG,&
    &              XRFCNT, LPTMP, LCTMP )

!.........  Process projection matrix that depends on pol/acts...
    IF( JFLAG .AND. LPTMP ) THEN

        CALL GENPROJ( PTMPDEV, CPYEAR, ENAME )

    ELSE IF ( JFLAG ) THEN
        MESG = 'WARNING: No records from the projection packet '//&
        &       'matched the inventory.'
        CALL M3MSG2( MESG )

    END IF

!.........  Process control matrices that depend on pollutants...

!.........  Multiplicative matrix
    IF( LCTMP .AND.&
    &  ( CFLAG .OR. GFLAG .OR. LFLAG .OR. MFLAG ) ) THEN

!.............  Write-out control matrix
        NCPE = MAX( PKTCNT( 2 ), PKTCNT( 7 ) )
        CALL GENMULTC( CTMPDEV, GTMPDEV, LTMPDEV, MTMPDEV,&
        &               NCPE, PYEAR, ENAME, MNAME, CFLAG, GFLAG,&
        &               LFLAG, MFLAG )
    ELSE IF ( CFLAG .OR. GFLAG .OR. LFLAG .OR. MFLAG ) THEN
        MESG = 'WARNING: No records from any control packet '//&
        &       'matched the inventory.'
        CALL M3MSG2( MESG )

    END IF

!.........  Do a final check that at least 1 tmp file had records written.
!           This is supposed to be checked by ALOCPKTS, but add here in case it
!           fails for some reason (like bug/code change error).
    IF( .NOT. LPTMP .AND. .NOT. LCTMP ) THEN

        MESG = 'ERROR: No records from any packet matched the '//&
        &       'inventory sources.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

!.........  Post-process temporary files to create final report file
    CALL WCNTLREP( CTMPDEV, GTMPDEV, LTMPDEV, MTMPDEV )

!.........  Successful completion
    CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END PROGRAM CNTLMAT

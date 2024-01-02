
SUBROUTINE RDINVMAP( INAME, IDEV, ENAME, ANAME, SDEV )

    !***********************************************************************
    !  program body starts at line
    !
    !  DESCRIPTION:
    !     Reads in inventory map file contents and checks them to ensure
    !     that they are consistent.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version ??/???? by ???
    !       Version 11/2023 by CJC:  USE M3UTILIO, ".f90" source format, and
    !       related changes
    !****************************************************************************/
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
    !.......  This module contains the information about the source category
    USE MODINFO,    ONLY: NMAP, MAPNAM, MAPFIL
    USE MODFILESET, ONLY: ALLFILES

    IMPLICIT NONE

    !.......   INCLUDES:

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters
    INCLUDE 'SETDECL.h90'       !  FileSetAPI variables and functions

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: INAME            ! map inventory name
    INTEGER     , INTENT (IN) :: IDEV             ! map inventory unit number
    CHARACTER(*), INTENT (IN) :: ENAME            ! I/O API inven logical name
    CHARACTER(*), INTENT( IN) :: ANAME            ! ASCII inven logical name
    INTEGER     , INTENT(OUT) :: SDEV             ! ASCII inven unit number

    !.......   EXTERNAL FUNCTIONS and their descriptions:

    INTEGER, EXTERNAL :: GETFLINE
    INTEGER, EXTERNAL :: GETIFDSC

    !.......   Other local variables
    INTEGER          L, L1, N, S        ! counters and indices
    INTEGER          J1, J2, J3, J4, J5

    INTEGER          IOS            ! i/o status
    INTEGER          IREC           ! record counter
    INTEGER          NHEAD          ! number of lines in header
    INTEGER          NLINE          ! number of lines in file
    INTEGER          NPOLCNT        ! data counter
    INTEGER          NSRCASCI       ! no. sources in ASCII file
    INTEGER          NSRC_TMP       ! no. sources in data files

    LOGICAL       :: EFLAG = .FALSE.      ! true: error found
    LOGICAL       :: RPFLAG = .FALSE.     ! true: time to read data lines

    CHARACTER(256)     MESG       ! message buffer
    CHARACTER(1024)    LINE       ! line buffer
    CHARACTER(NAMLEN3) PNAME      ! logical file name for data files
    CHARACTER(PHYLEN3) APHYS      ! ASCII physical file name
    CHARACTER(PHYLEN3) EPHYS      ! I/O API physical file name
    CHARACTER(PHYLEN3) PATH       ! path name

    CHARACTER(16), PARAMETER :: PROGNAME = 'RDINVMAP'       !  program name

    !***********************************************************************
    !   begin body of subroutine RDINVMAP

    !......  Determine path of map file, which will set relative path
    !          for files listed in map
    MESG = 'Inventory map file'
    CALL ENVSTR( INAME, MESG, ' ', APHYS, IOS )
    IF( IOS .NE. 0 ) THEN
        MESG = 'Unable to evaluate environment variable "' // TRIM( INAME ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    L = LEN_TRIM( APHYS )
    DO N = L, 1, -1

        IF( APHYS( N:N ) .EQ. '/' .OR.      &
            APHYS( N:N ) .EQ. '\'      ) THEN
            PATH = APHYS( 1:N )
            EXIT
        END IF

    END DO

    !......  Read map file and check that all of the components are
    !          in place
    IREC = 0
    NPOLCNT = 0
    DO

        READ( IDEV, 93000, END=111, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF ( IOS .NE. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'I/O error', IOS,   &
                'inventory map file at line', IREC
            CALL M3MESG( MESG )
            CYCLE
        END IF

        !......  Skip blank lines (not that there should be any)
        IF( LINE .EQ. ' ' ) CYCLE

        L1 = LEN_TRIM( LINE )

        !...... Check for fields that we expect
        J1 = INDEX( LINE, '/TEXT/' )
        J2 = INDEX( LINE, '/IOAPI/' )
        J3 = INDEX( LINE, '/NDAT/' )
        J4 = INDEX( LINE, '/DATMAP/' )
        J5 = INDEX( LINE, '/END/' )

        !......  Proceed as needed for each packet...
        !......  For /TEXT/ packet, retrieve physical file name
        IF( J1 .GT. 0 ) THEN

            IF( L1 .LE. J1+7 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: map file format ' //    &
                       'corrupted for /TEXT/ packet at line', IREC
                CALL M3MSG2( MESG )
            ELSE
                APHYS = TRIM( PATH ) // LINE( J1+7:L1 )
            END IF

        !......  For /IOAPI/ packet, retrieve physical file name
        ELSE IF( J2 .GT. 0 ) THEN

            IF( L1 .LE. J2+8 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: map file format ' //    &
                       'corrupted for /IOAPI/ packet at line', IREC
                CALL M3MSG2( MESG )
            ELSE
                EPHYS = TRIM( PATH ) // LINE( J2+8:L1 )
            END IF

        !......  For /NDAT/ packet, retrieve number of pollutants
        ELSE IF( J3 .GT. 0 ) THEN

            IF( L1 .LE. J3+7 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: map file format ' //    &
                       'corrupted for /NPOL/ packet at line', IREC
                CALL M3MSG2( MESG )
            ELSE
                NMAP = STR2INT( LINE( J3+7:L1 ) )

            END IF

        !......  For /DATMAP/ packet, set special flag, allocate
        !              memory for pollutant names and files, and initialize
        ELSE IF( J4 .GT. 0 ) THEN
            RPFLAG = .TRUE.

            IF( ALLOCATED( MAPNAM ) ) DEALLOCATE( MAPNAM, MAPFIL )
            ALLOCATE( MAPNAM( NMAP ),       &
                      MAPFIL( NMAP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MAPFIL', PROGNAME )
            MAPNAM = ' '
            MAPFIL = ' '
            CYCLE

        !......  For /END/ packet, ensure count of pollutants in file
        !              matched the number listed in the /NPOL/ packet
        ELSE IF( J5 .GT. 0 ) THEN
            RPFLAG = .FALSE.

            IF( NPOLCNT .NE. NMAP ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Inventory map file '//      &
                       'is corrupted. Number of data values '//         &
                       CRLF()// BLANK10// 'expected by /NDAT/ was',     &
                       NMAP, ', but ', NPOLCNT,' records found in file.'
                CALL M3MSG2( MESG )
            END IF

        END IF

        !......  If pollutant read flag is true, then read pollutant
        !              name and store it's file.
        IF( RPFLAG ) THEN

            NPOLCNT = NPOLCNT + 1
            MAPNAM( NPOLCNT ) = LINE( 1:NAMLEN3 )
            MAPFIL( NPOLCNT ) = TRIM( PATH ) // LINE( NAMLEN3+2:L1 )

        END IF

    END DO
111 CONTINUE     ! exit from read loop at the end

    !......  If map inventory file empty then abort
    IF( EFLAG ) THEN
        MESG = 'ERROR: Map-formatted inventory file is empty.' //   &
               CRLF()// BLANK10// 'Please manually delete it '//    &
               'and the *.ncf and *_dat/*.ncf files that'//         &
               CRLF()// BLANK10// 'should be listed in it.'
        CALL M3MSG2( MESG )
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
    END IF

    !......  If error reading file, then abort
    IF( EFLAG ) THEN
        MESG = 'Corrupted map-formatted inventory file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !......  Make sure all of the lines were there
    IF( APHYS .EQ. ' ' .OR.         &
        EPHYS .EQ. ' ' .OR.         &
        NMAP  .EQ. 0   .OR.         &
        SIZE(MAPNAM).EQ.0 .OR.      &
        SIZE(MAPFIL).EQ.0 .OR.      &
        RPFLAG                  ) THEN
        MESG = 'Corrupted map formatted file - missing lines'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF (MAPNAM(1) .EQ. ' ' .OR.     &
        MAPFIL(1) .EQ. ' ') THEN
        MESG = 'Corrupted map formatted file - missing lines'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !......  Open ASCII file
    IF( .NOT. SETENVVAR( ANAME, APHYS ) ) THEN
        EFLAG = .TRUE.
        MESG = 'ERROR: Could not set logical file name for file:' //    &
               CRLF() // BLANK10 // TRIM( EPHYS )
        CALL M3MSG2( MESG )

    ELSE

        SDEV = GETEFILE( ANAME, .TRUE., .TRUE., PROGNAME )

        IF ( SDEV .LT. 0 ) THEN         !  failure to open

            EFLAG = .TRUE.
            MESG = 'ERROR: Could not open ASCII inventory file ' // &
               'listed in map ' // CRLF() // BLANK10 //             &
               'inventory for file name: ' //                       &
               CRLF() // BLANK10 // TRIM( APHYS )
            CALL M3MSG2( MESG )

        END IF          !  if getefile() failed

    END IF

    !......  Set environment variable for
    IF( .NOT. SETENVVAR( ENAME, EPHYS ) ) THEN
        EFLAG = .TRUE.
        MESG = 'ERROR: Could not set logical file name for file:' //    &
               CRLF() // BLANK10 // TRIM( EPHYS )
        CALL M3MSG2( MESG )

    !......  Open I/O API file
    ELSE IF( .NOT. OPENSET( ENAME, FSREAD3, PROGNAME ) ) THEN
        EFLAG = .TRUE.
        MESG = 'ERROR: Could not open I/O API inventory file ' //   &
               'listed in map ' // CRLF() // BLANK10 //             &
               'inventory for file name: ' //                       &
               CRLF() // BLANK10 // TRIM( EPHYS )
        CALL M3MSG2( MESG )

    END IF

    !......  Abort if errors opening input files
    IF( EFLAG ) THEN
        MESG = 'Problem attempting to open inventory files.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !......  Get number of sources from ASCII file
    NLINE = GETFLINE( SDEV, 'ASCII inventory file' )
    READ( SDEV, * ) NHEAD
    NSRCASCI = NLINE - NHEAD - 1

    !......  Rewind ASCII file
    REWIND( SDEV )

    !......  Get number of sources from I/O API file
    IF( .NOT. DESCSET( ENAME, ALLFILES ) ) THEN
        MESG = 'Could not get description of file "' // TRIM( ENAME ) // '".'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !......  Make sure that the number of sources is consistent between
    !          the I/O API and ASCII parts of the inventory file.
    IF( NROWS3D .NE. NSRCASCI ) THEN
        WRITE( MESG,94010 ) 'ERROR: Number of sources in I/O API '//    &
               'and ASCII inventory files ' // CRLF() // BLANK10 //     &
               'are inconsistent.  I/O API file has', NROWS3D,          &
               'but ACSII file has', NSRCASCI, '.'
        CALL M3MSG2( MESG )
        CALL M3EXIT( PROGNAME, 0, 0, 'Corrupted input files', 2 )
    END IF

    !......  Loop through all pollutant files...
    DO N = 1, NMAP

        !......  Set environment variable for input file
        PNAME = 'TMP_POL_FILE'
        IF( .NOT. SETENVVAR( PNAME, MAPFIL( N ) ) ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: Could not set logical file name ' //     &
                       CRLF() // BLANK10 // 'for file ' // TRIM( EPHYS )
            CALL M3MSG2( MESG )

        !......  Open pollutant file with physical file name
        ELSE IF( .NOT. OPENSET( PNAME, FSREAD3, PROGNAME ) ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: Could not open "'// TRIM( MAPNAM(N) ) // &
               '" file listed in ' // CRLF() // BLANK10 //          &
               'map-formatted inventory for file name: ' //         &
               CRLF() // BLANK10 // TRIM( MAPFIL( N ) )
            CALL M3MSG2( MESG )

        END IF

        !......  Get number of sources from pollutant file and compare
        IF( .NOT. DESCSET( PNAME, ALLFILES ) ) THEN
            MESG = 'Could not get description of file "' // TRIM( PNAME ) // '".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        !......  Check that the number of sources is consistent with
        !              the other files
        NSRC_TMP = GETIFDSC( FDESC3D, '/NSRC/', .TRUE. )
        IF( NSRC_TMP .NE. NSRCASCI ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Number of sources in main'//    &
               ' inventory files and "'// TRIM( MAPNAM(N) )// '"'//     &
               CRLF()// BLANK10 // 'file are inconsistent. ' //         &
               'Inventory files have', NSRCASCI, ' sources but "' //    &
               TRIM( MAPNAM(N) )// '" file has', NSRC_TMP, '.'
            CALL M3MSG2( MESG )
        END IF

        IF( .NOT. CLOSESET( PNAME ) ) THEN
            MESG = 'Could not close file:'//CRLF()//BLANK10// MAPFIL( N )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

    END DO       ! end of loop through data files

    IF( EFLAG ) THEN
        MESG = 'Inconsistent number of sources.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

    !.......   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I10, :, 1X ) )

94100 FORMAT( 1X,  A, 1X,  A, I7 )

END SUBROUTINE RDINVMAP



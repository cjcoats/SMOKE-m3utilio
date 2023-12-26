
PROGRAM INVSPLIT

!***********************************************************************
!  program body starts at line
!
!  DESCRIPTION:
!       This program splits an IDA and toxics inventory file into
!       multiple output files
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!       Models-3 I/O
!       PROMPTFFILE, ENVINT, STR2INT, FIND1, INDEX1
!
!  REVISION  HISTORY:
!       Created  12/02 by M Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!***********************************************************************
!
! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
!                System
! File: $Id$
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

!...........   MODULES for public variables
!.........  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, CRL

    IMPLICIT NONE

!...........   INCLUDES:

    INCLUDE 'EMCNST3.EXT'     ! emissions constant parameters

!...........   PARAMETERS and their descriptions:

    CHARACTER(50), PARAMETER :: SCCSW = '%W%'

!.........  EXTERNAL FUNCTIONS and their descriptions:
    INTEGER, EXTERNAL :: GETFLINE
    INTEGER, EXTERNAL :: GETFORMT

!...........   LOCAL PARAMETERS
    INTEGER,       PARAMETER :: NSEG     = 80
    CHARACTER(16), PARAMETER :: PROGNAME = 'INVSPLIT'   !  program name

!.........  Allocatable arrays...

!.........  Splits file arrays
    INTEGER, ALLOCATABLE :: SINDEX  ( : )      ! sorting index
    INTEGER, ALLOCATABLE :: STATLIST( : )      ! state codes
    INTEGER, ALLOCATABLE :: FNUMLIST( : )      ! file numbers
    INTEGER, ALLOCATABLE :: OUTNUM  ( : )      ! output file codes

!.........  Logical file names and unit numbers
    INTEGER, ALLOCATABLE :: ODEV( : )  ! output unit numbers
    INTEGER   IDEV         ! input inventory file
    INTEGER   LDEV         ! log file unit number
    INTEGER   SDEV         ! log file unit number

!.........  Local arrays
    CHARACTER(40) SEGMENT( NSEG )

!.........  Local variables

    INTEGER   I, CNT, F, J, PF         ! indices and counters

    INTEGER   IFMT         ! inventory format code
    INTEGER   IOS          ! i/o status
    INTEGER   IREC         ! record counter
    INTEGER   NLIST        ! number of entries in split file
    INTEGER   NOUT         ! number of output files
    INTEGER   STA          ! tmp state ID

    LOGICAL :: EFLAG = .FALSE.   ! true: error found

    CHARACTER(256)  MESG    ! temporary message array
    CHARACTER(2560) LINE    ! line buffer

    CHARACTER(NAMLEN3) INNAME   ! input file name
    CHARACTER(NAMLEN3) OUTNAME  ! output file names

!***********************************************************************
!   begin body of program INVSPLIT

    LDEV = INIT3()

!.........  Write out copyright, version, web address, header info, and prompt
!           to continue running the program.

    CALL INITEM( LDEV, SCCSW, PROGNAME )

!.........  Set source category based on environment variable setting
    CALL GETCTGRY

!.........  Find name name of raw inventory file
    J = INDEX1( CATEGORY, NCAT, CATLIST )
    IF( J .LE. 0 ) THEN
        MESG = 'INTERNAL ERROR: Do not know about category ' //&
        &       TRIM( CATEGORY ) // ' in program ' // PROGNAME
        CALL M3MSG2( MESG )
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ELSE
        INNAME = ANAMLIST( J )

    END IF

!.........  Prompt for name of input inventory file
    MESG = 'Enter logical name of the RAW ' //&
    &       TRIM( CATEGORY ) // ' AVERAGE INVENTORY ' // 'file'

    IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INNAME, PROGNAME )

!.........  Prompt for name of input splits definitions file
    MESG = 'Enter logical name of the SPLITS DEFINITIONS file'
    SDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'SPLITS', PROGNAME )

!.........   Get format of input file
    IFMT = GETFORMT( IDEV, -1 )

!.........   Get size for splits file
    NLIST = GETFLINE( SDEV, 'Splits definitions' )

!.........  Allocate memory for arrays
    ALLOCATE( STATLIST( NLIST ),&
    &          FNUMLIST( NLIST ),&
    &            SINDEX( NLIST ), STAT=IOS )
    CALL CHECKMEM( IOS, 'STATLIST,,,SINDEX', PROGNAME )

!.........  Read split definitions and create list of output file
!           numbers
    IREC = 0
    CNT = 0
    DO I = 1, NLIST

        READ( SDEV, 93000, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF ( IOS .GT. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010)&
            &    'I/O error', IOS, 'reading splits definitions '//&
            &    'file at line', IREC
            CALL M3MESG( MESG )
            CYCLE
        END IF

        IF( LINE( 1:1 ) .EQ. CINVHDR ) CYCLE

        CALL PARSLINE( LINE, NSEG, SEGMENT )

        CNT = CNT + 1
        IF( CNT .LE. NLIST ) THEN
            SINDEX  ( CNT ) = CNT
            STATLIST( CNT ) = STR2INT( SEGMENT( 1 ) )
            FNUMLIST( CNT ) = STR2INT( SEGMENT( 2 ) )
        END IF

    END DO
    NLIST = CNT

!.........  Exit if problem with splits input file
    IF( EFLAG ) THEN
        MESG = 'Problem reading splits definitions file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Create list of output file numbers and count the number
!           of outputs...
!.........  Sort the number of output files
    CALL SORTI1( NLIST, SINDEX, FNUMLIST )

!.........  Count number of outputs
    PF = IMISS3
    NOUT = 0
    DO I = 1, NLIST
        J = SINDEX( I )
        F = FNUMLIST( J )
        IF( F .NE. PF ) NOUT = NOUT + 1
        PF = F
    END DO

!.........  Allocate memory for output numbers
    ALLOCATE( OUTNUM( NOUT ),&
    &            ODEV( NOUT ), STAT=IOS )
    CALL CHECKMEM( IOS, 'OUTNUM...ODEV', PROGNAME )

!.........  Store the output numbers
    PF = IMISS3
    NOUT = 0
    DO I = 1, NLIST
        J = SINDEX( I )
        F = FNUMLIST( J )
        IF( F .NE. PF ) THEN
            NOUT = NOUT + 1
            OUTNUM( NOUT ) = F
        END IF
        PF = F
    END DO

!.........  Prompt for output file names
    DO I = 1, NOUT
        WRITE( OUTNAME, '(A,I2.2)' ) 'OUTFILE', OUTNUM( I )
        WRITE( MESG, '(A,1X,I2.2)' ) 'Enter logical name ' //&
        &       'for output file', OUTNUM( I )
        ODEV( I ) = PROMPTFFILE( MESG, .FALSE., .TRUE., OUTNAME,&
        &                         PROGNAME )
    END DO

!.........  Loop through input file and write output files
    DO

        READ( IDEV, 93000, END=199, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF ( IOS .GT. 0 ) THEN
            WRITE( MESG, 94010 )&
            &    'I/O error', IOS,&
            &    'reading inventory file at line', IREC
            CALL M3MESG( MESG )
            EXIT
        END IF

!.............  Skip header lines and write them to the output file
        IF( LINE( 1:1 ) .EQ. CINVHDR ) THEN
            DO I = 1, NOUT
                WRITE( ODEV( I ), '(A)' ) TRIM( LINE )
            END DO
            CYCLE
        END IF

!.............  Convert state code to integer
        SELECT CASE ( IFMT )

          CASE ( ORLFMT, ORLNPFMT )
            CALL PARSLINE( LINE, NSEG, SEGMENT )
            STA = INT( STR2INT( SEGMENT(1) ) / 1000 )

          CASE ( FF10FMT )
            CALL PARSLINE( LINE, NSEG, SEGMENT )
            STA = INT( STR2INT( SEGMENT(2) ) / 1000 )

          CASE DEFAULT
            WRITE( MESG,94010 ) 'Cannot split file with ' //&
            &                    TRIM( FMTNAMES( IFMT ) ) // 'format'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END SELECT

!.............  Find state in list
        J = FIND1( STA, NLIST, STATLIST )

!.............  Error if state is not found
        IF( J .LE. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: State', STA, 'is not ' //&
            &       'found in splits file.' // CRLF() // BLANK10 //&
            &       'Ensure that the state is in the list and  '//&
            &       'the states are ' // CRLF() // BLANK10 //&
            &       'listed in sorted order.'
            CALL M3MSG2( MESG )
        END IF

!.............  No point in writing output files if there is an error
        IF( EFLAG ) CYCLE

!.............  Check if output file number
        F = FNUMLIST( J )

!.............  Write line to correct output file
        WRITE( ODEV( F ), '(A)' ) TRIM( LINE )

    END DO

199 CONTINUE   ! exit from read loop

!.........  Check if error found
    IF( IOS .GT. 0 ) THEN
        MESG = 'Problem reading input file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Normal completion of program
    CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )


!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10 ( A, :, I5, :, 2X ) )


END PROGRAM INVSPLIT



SUBROUTINE RDEPROC( FDEV )

!***********************************************************************
!
!  DESCRIPTION:
!     This subroutine reads the emission processes file, which contains columns
!     for the activity, associated process, and associated pollutants.  If
!     there is more than one process per activity, then these are listed on
!     separate lines in the file. If there is more than one pollutant per
!     activity and process, then these are listed in additional columns.
!     During the read, the column number is set dynamically.  Only processes
!     for activities that are in the inventory are read.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 10/99 by M. Houyoux
!       Modified  3/14 by B.H. Baek
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
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

!.........  MODULES for public variables
!.........  This module contains emission factor tables and related
    USE MODEMFAC, ONLY: EMTNAM, EMTIDX, NETYPE, MXETYPE, NEPOL, EMTPOL

!.........  This module contains the information about the source category
    USE MODINFO, ONLY: NIACT, ACTVTY

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   EXTERNAL FUNCTIONS and their descriptions:
    INTEGER, EXTERNAL :: GETFLINE
    INTEGER, EXTERNAL :: GETNLIST
    LOGICAL, EXTERNAL :: BLKORCMT

!...........   SUBROUTINE ARGUMENTS
    INTEGER, INTENT (IN) :: FDEV   ! cross-reference file unit no.

!...........   Local allocatable arrays
    INTEGER           , ALLOCATABLE :: INDX   ( : )  ! POLA sorting indx
    CHARACTER(NAMLEN3), ALLOCATABLE :: SEGMENT( : )  ! line segments
    CHARACTER(NAMLEN3), ALLOCATABLE :: POLNAM ( : )  ! pol names
    CHARACTER(NAMLEN3), ALLOCATABLE :: POLA   ( : )  ! all unsorted pols

!...........   Local parametes
    CHARACTER    , PARAMETER :: CONTCHAR = '\'  ! line continuation character
    CHARACTER(16), PARAMETER :: PROGNAME = 'RDEPROC' ! program name

!...........   Other local variables
    INTEGER         I, J, K, L, M, V    !  counters and indices

    INTEGER         IOS     !  i/o status
    INTEGER         L1, L2  !  tmp string lengths
    INTEGER         LJ      !  string length for emis type joiner
    INTEGER         MXCOLS  !  maximum number of columns in the file
    INTEGER         MXPOL   !  maximum number of pollutants per process
    INTEGER         NCOLS   !  no. columns in a row
    INTEGER         NLINES  !  number of lines in file
    INTEGER         NPOL    !  no. pollutants in a row
    INTEGER         NPUNSRT !  no. pols
    INTEGER         NTLINES !  no. of lines taking into account continuation lines

    LOGICAL      :: FNDPOL   = .FALSE.   ! true: found pollutant on master list
    LOGICAL      :: NEWLINE  = .TRUE.    ! true: current line is new (not continued)
    LOGICAL      :: EFLAG    = .FALSE.   ! true: error found

    CHARACTER(3000)    LINE     !  line buffer
    CHARACTER(300)     MESG     !  message buffer
    CHARACTER(NAMLEN3) ACT      !  tmp activity name
    CHARACTER(NAMLEN3) CPOL     !  tmp pollutant name
    CHARACTER(NAMLEN3) LPOL     !  tmp pollutant name from previous iter

!***********************************************************************
!   begin body of subroutine RDEPROC

!.........  Get the number of lines for the file and allocate array so that
!           the type of the line can be stored
    NLINES = GETFLINE( FDEV, 'Emission processes file' )
    NTLINES = NLINES

!.........  Read through file to determine the maximum no. of columns and
!           check packet information
    MXCOLS   = 0

    DO I = 1, NLINES

        READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE

        IF ( IOS .NE. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 )&
            &    'I/O error', IOS,&
            &    'reading emission processes file at line', I
            CALL M3MESG( MESG )
            CYCLE
        END IF

!.............  Skip any blank and comment lines
        IF( BLKORCMT( LINE ) ) THEN
            NTLINES = NTLINES - 1
            CYCLE
        END IF

        LINE = ADJUSTL( LINE )
        L1 = LEN_TRIM( LINE )

!.............  If previous line was continued, add this lines columns to number
!               from previous line; otherwise, count columns in just this line
        IF( .NOT. NEWLINE ) THEN
            NCOLS = NCOLS + GETNLIST( L1, LINE )
        ELSE
            NCOLS = GETNLIST( L1, LINE )
        END IF

!.............  Check for continuation character in current line
        IF( LINE( L1:L1 ) == CONTCHAR ) THEN
            NEWLINE = .FALSE.
            NTLINES = NTLINES - 1
            NCOLS = NCOLS - 1
        ELSE
            NEWLINE = .TRUE.
        END IF

        IF( NCOLS > MXCOLS ) THEN
            MXCOLS = NCOLS
        END IF

    END DO

!.........  Check for errors so far
    IF( EFLAG ) THEN
        MESG = 'Problem reading emission processes file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    MXPOL = MXCOLS - 1

!.........  Allocate memory for parsing line segements and storing pollutants
!.........  Allocate memory for emission types and count for each per activity
    ALLOCATE( SEGMENT( MXCOLS ),&
    &           POLNAM( MXPOL ),&
    &           EMTNAM( NTLINES*MXPOL, NIACT ),&
    &           EMTIDX( NTLINES*MXPOL, NIACT ),&
    &           NETYPE( NIACT ), STAT=IOS )
    CALL CHECKMEM( IOS, 'NETYPE', PROGNAME )

    EMTNAM = ' '  ! array
    EMTIDX = 0    ! array
    NETYPE = 0    ! array

!.........  Rewind file
    REWIND( FDEV )

    NEWLINE  = .TRUE.

!.........  Store contents of emissions processes file in output order
    J = 0
    L = 0
    M = 0
    DO I = 1, NLINES

        READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE

!.............  Skip blank and comment lines
        IF( BLKORCMT( LINE ) ) CYCLE

        LINE = ADJUSTL( LINE )
        L1 = LEN_TRIM( LINE )

!.............  Separate line into segments
        NCOLS = GETNLIST( L1, LINE )
        CALL PARSLINE( LINE, NCOLS, SEGMENT )

        IF( NEWLINE ) THEN
            ACT = SEGMENT( 1 )

            IF( SEGMENT( NCOLS ) == CONTCHAR ) THEN
                NPOL = NCOLS - 2
                POLNAM( 1:NPOL ) = SEGMENT( 2:NCOLS-1 )
            ELSE
                NPOL = NCOLS - 1
                POLNAM( 1:NPOL ) = SEGMENT( 2:NCOLS )
            END IF
        ELSE
            IF( SEGMENT( NCOLS ) == CONTCHAR ) THEN
                NPOL = NCOLS - 1
                POLNAM( 1:NPOL ) = SEGMENT( 1:NCOLS-1 )
            ELSE
                NPOL = NCOLS
                POLNAM( 1:NPOL ) = SEGMENT( 1:NCOLS )
            END IF
        END IF

!.............  Make sure activity is in the inventory
        K = INDEX1( ACT, NIACT, ACTVTY )

!.............  Store emission processes and associated pollutants
        IF( K .GT. 0 ) THEN

            DO V = 1, NPOL
                J = J + 1
                EMTNAM( J,K ) = POLNAM( V )
            END DO

            NETYPE( K ) = NETYPE( K ) + NPOL

        END IF

!.............  Check if current line is continued
        IF( SEGMENT( NCOLS ) == CONTCHAR ) THEN
            NEWLINE = .FALSE.
        ELSE
            NEWLINE = .TRUE.
        END IF

    END DO      ! End loop over file

!.........  Set the maximum number of emission types
    MXETYPE = MAXVAL( NETYPE )

!.........  Create a list of pollutants associated with the emission types...

!.........  Allocate memory for unsorted pollutants list
    ALLOCATE( POLA( MXETYPE * NIACT ),&
    &          INDX( MXETYPE * NIACT ), STAT=IOS )
    CALL CHECKMEM( IOS, 'INDX', PROGNAME )

!.........  Create unsorted pollutants list
    J = 0
    DO I = 1, NIACT

        DO K = 1, NETYPE( I )
            J = J + 1
            POLA( J ) = EMTNAM( K,I )
            INDX( J ) = J
        END DO

    END DO
    NPUNSRT = J

!.........  Sort pollutants
    CALL SORTIC( NPUNSRT, INDX, POLA )

!.........  Determine number of actual pollutants
    LPOL = '-9'
    K = 0
    DO I = 1, NPUNSRT

        J = INDX( I )
        CPOL = POLA( J )

        IF( CPOL .NE. LPOL ) THEN
            K = K + 1
            LPOL = CPOL
        END IF

    END DO

    NEPOL = K

!.........  Allocate memory for sorted pollutants
    ALLOCATE( EMTPOL( NEPOL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EMTPOL', PROGNAME )

!.........  Store sorted pollutants in a unique list
    LPOL = '-9'
    K = 0
    DO I = 1, NPUNSRT

        J = INDX( I )
        CPOL = POLA( J )

        IF( CPOL .NE. LPOL ) THEN
            K = K + 1
            EMTPOL( K ) = CPOL
            LPOL = CPOL
        END IF

    END DO

!.........  Create an index that references between the emission types and the
!           list of unique pollutants
    DO I = 1, NIACT

        DO K = 1, NETYPE( I )

            CPOL = EMTNAM( K,I )
            J = INDEX1( CPOL, NEPOL, EMTPOL )
            EMTIDX( K,I ) = J

        END DO

    END DO

!.........  Rewind file
    REWIND( FDEV )

!.........  Deallocate memory for local arrays
    DEALLOCATE( SEGMENT, POLNAM, POLA, INDX )

    RETURN

!.........  Error message for reaching the end of file too soon
999 MESG = 'End of file reached unexpectedly. ' //&
    &       'Check format of emission processes file.'
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDEPROC

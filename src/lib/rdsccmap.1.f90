
SUBROUTINE RDSCCMAP( ADEV )

!***********************************************************************
!  subroutine body starts at line 99
!
!  DESCRIPTION:
!       Reads SCC_MAP file that maps aggregated SCC to full SCCs.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:  none
!
!  REVISION  HISTORY:
!       Created  5/2014: by B.H. Baek
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
!***********************************************************************
    USE M3UTILIO

!.........  MODULES for public variables
!.........  This module is used for MOBILE6 setup information
    USE MODMOBIL, ONLY: NSCCMAP, SCCMAPLIST

    IMPLICIT NONE

!...........   INCLUDES

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   SUBROUTINE ARGUMENTS
    INTEGER, INTENT (IN) :: ADEV     ! COUNTY_FUELMONTH file unit no.

!...........   EXTERNAL FUNCTIONS and their descriptions:
    LOGICAL, EXTERNAL :: BLKORCMT
    INTEGER, EXTERNAL :: GETFLINE

!...........   Local allocatable arrays
    INTEGER, ALLOCATABLE :: IDX ( : )          ! index of SCCs
    CHARACTER( SCCLEN3 ), ALLOCATABLE :: SCCMAPRAW ( :,: )  ! raw mapped SCCs

!...........   Local arrays
    CHARACTER( SCCLEN3 )  SEGMENT( 3 )  ! parsed input line

!...........   Other local variables
    INTEGER I, J, N, NS               ! counters and indices

    INTEGER    IOS                    ! I/O status
    INTEGER :: IREC = 0               ! record counter
    INTEGER :: NLINES = 0             ! number of lines
    INTEGER :: NSCC, NFSCC            ! no of SCCs

    LOGICAL      :: EFLAG   = .FALSE.   ! true: error found

    CHARACTER(100)       LINE     !  line buffer
    CHARACTER(300)       MESG     !  message buffer
    CHARACTER( 8 )       CNFSCC   !  tmp buffer
    CHARACTER( SCCLEN3 ) CURSCC, PRVSCC, FULLSCC   ! current, previous, full SCCs

    CHARACTER(16) :: PROGNAME = 'RDSCCMAP'   ! program name

!***********************************************************************
!   begin body of subroutine RDSCCMAP

!.........  Get the number of lines in the file
    MESG = 'Aggregated SCC to full SCC map input'
    NLINES = GETFLINE( ADEV,MESG )

!.........  Allocate memory to store settings information
    ALLOCATE(       IDX( NLINES ),&
    &          SCCMAPRAW( NLINES,2 ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SCCMAPRAW', PROGNAME )

!.........  Initialize arrays
    IDX = 0
    NSCC = 0
    SCCMAPRAW = ' '
    DO I = 1, NLINES

        IDX( I ) = I

!.........  Read line
        READ( ADEV, 93000, END=999, IOSTAT=IOS ) LINE

        IREC = IREC + 1

        IF ( IOS /= 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010 )&
            &    'I/O error', IOS, 'reading reference county ' //&
            &    'fuel month settings file at line', IREC
            CALL M3MESG( MESG )
            CYCLE
        END IF

!.............  Skip blank/comment lines
        IF( BLKORCMT( LINE ) ) CYCLE

!.............  Parse the line into segments
        CALL PARSLINE( LINE, 2, SEGMENT )

!.............  Store values in unsorted array
        NSCC = NSCC + 1
        SCCMAPRAW( I,1 ) = SEGMENT( 2 )    ! referenced SCCs
        SCCMAPRAW( I,2 ) = SEGMENT( 1 )    ! full SCCs

        CALL PADZERO( SCCMAPRAW( I,1 ) )
        CALL PADZERO( SCCMAPRAW( I,2 ) )

    END DO  ! done reading SCC map input file

    CLOSE( ADEV )

!.........  Abort if error found while reading settings file
    IF( EFLAG ) THEN
        MESG = 'Problem reading SCC mapping input file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  If the no. of lines in the file is less than the no. of ref.
!              counties, then something is wrong, but we'll go through
!              the function to print out error messages
    ALLOCATE( SCCMAPLIST ( NSCC,3 ), STAT=IOS )
    CALL CHECKMEM( IOS, 'FMREFLIST', PROGNAME )
    SCCMAPLIST = ''

    CALL SORTIC( NLINES, IDX, SCCMAPRAW( :,1 ) )    ! sort ref SCCs

!.........  Count no of aggregated SCCs
    PRVSCC = ''
    NSCC   = 0
    DO I = 1, NLINES

        J = IDX( I )
        CURSCC  = SCCMAPRAW( J,1 )
        FULLSCC = SCCMAPRAW( J,2 )

        IF( CURSCC == ' ' ) CYCLE

        NSCC = NSCC + 1

        IF( PRVSCC /= CURSCC ) THEN
            NFSCC = 0
            NS = NSCC
        ELSE
            NFSCC = NFSCC + 1
        END IF

        WRITE( CNFSCC,'(I8)' ) NFSCC
        SCCMAPLIST(    NSCC,1 ) = CURSCC
        SCCMAPLIST(    NSCC,2 ) = FULLSCC
        SCCMAPLIST( NS:NSCC,3 ) = CNFSCC

        PRVSCC = CURSCC

    END DO

    NSCCMAP = NSCC

!.........  Deallocate local memory
    DEALLOCATE( IDX, SCCMAPRAW )

    RETURN

999 MESG = 'End of file reached unexpectedly. ' //&
    &       'Check format of COUNTY_FUELMONTH' // CRLF() // BLANK5 //&
    &       'input file.'
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx
93000 FORMAT( A )

!...........   Internal buffering formats............ 94xxx
94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDSCCMAP

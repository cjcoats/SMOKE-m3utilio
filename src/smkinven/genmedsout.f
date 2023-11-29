
        SUBROUTINE GENMEDSOUT( FDEV, FNAME, TZONE, TYPNAM )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine loops through a list of day-specific or hour-specific files
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutine
C
C  REVISION  HISTORY:
C       Created 12/2013 by B.H. Baek
C       Version 11/2023 by CJC:  USE M3UTILIO and related changes
C
C*************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$
C
C***************************************************************************
        USE M3UTILIO

C.........  MODULES for public variables
C.........  This module contains the information about the source category
        USE MODINFO, ONLY: NIPPA, NCHARS, NSRC

C.........  This module contains data for day- and hour-specific data
        USE MODDAYHR, ONLY: LPDSRC, IDXSRC, PDEMOUT

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: FDEV      ! hour-specific file unit no.
        CHARACTER(*), INTENT (IN) :: FNAME     ! logical file name
        INTEGER     , INTENT (IN) :: TZONE     ! output time zone
        CHARACTER(*), INTENT (IN) :: TYPNAM    ! 'day' or 'hour'

C.........  EXTERNAL FUNCTIONS
        INTEGER, EXTERNAL :: GETFLINE

C...........   Local file formats
        INTEGER, ALLOCATABLE, SAVE :: FILFMT( : )  ! file format code

C...........   Character strings of day- or hr-specific list file
        CHARACTER(300), ALLOCATABLE, SAVE :: LSTSTR( : )

C.........  Local arrays

C...........   Other local variables
        INTEGER          I, J, K, L, N, S, T

        INTEGER          IOS, IREC            ! i/o status
        INTEGER          NFILE, IFIL          ! number of MEDS files
        INTEGER          IDEV                 ! input file unit no.


        LOGICAL       :: LASTFLAG = .FALSE.  ! true: process last inv file
        LOGICAL       :: DFLAG    = .FALSE.  ! true: day-specific
        LOGICAL       :: EFLAG    = .FALSE.  ! true: error found

        CHARACTER(300) :: MESG = ' '         ! message buffer
        CHARACTER(200)   NAMTMP              ! file name buffer

        CHARACTER(16) :: PROGNAME = 'GENMEDSOUT' !  program name

C***********************************************************************
C   begin body of program GENMEDSOUT

C.........  Allocate memory for daily or hourly output arrays
        ALLOCATE( LPDSRC( NSRC ),
     &            PDEMOUT( NSRC,NIPPA ),
     &            IDXSRC( NSRC,1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IDXSRC', PROGNAME )
        LPDSRC  = .FALSE.
        PDEMOUT = 0.0
        DO S = 1, NSRC
            IDXSRC( S,1 ) = S
        END DO

C.........  Get number of lines of inventory files in list format
        MESG = 'Processing ' // TYPNAM // '-specific data...'
        NFILE = GETFLINE( FDEV, MESG )

C.........  Determine format of day- or hour-specific inputs and store file names
C.........  Allocate memory for storing file formats
        IF( ALLOCATED( FILFMT ) ) DEALLOCATE( FILFMT )
        IF( ALLOCATED( LSTSTR ) ) DEALLOCATE( LSTSTR )
        ALLOCATE( FILFMT( NFILE ),
     &            LSTSTR( NFILE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LSTSTR', PROGNAME )

        FILFMT = 0  ! array

C.........  Store lines of day-specific list file
        CALL RDLINES( FDEV, MESG, NFILE, LSTSTR )

C.........  Get file format for files listed in day-specific list file
        CALL CHKLSTFL( NFILE, FNAME, LSTSTR, FILFMT )

        IREC = 0
        DO IFIL = 1, NFILE

            NAMTMP = LSTSTR( IFIL )

C.............  If line is LIST header, skip to next line
            IF( INDEX( NAMTMP, 'LIST' ) > 0 ) CYCLE

C.............  Open files, and report status
            IDEV = JUNIT()
            OPEN( IDEV, ERR=1003, FILE=NAMTMP, STATUS='OLD' )

            WRITE( MESG,94010 ) 'Successful open ' //
     &             'for emissions file:' // CRLF() // BLANK5 //
     &             NAMTMP( 1:LEN_TRIM( NAMTMP ) )
            CALL M3MSG2( MESG )

C.............  Read day-orhour-specific MEDS files
            IF( IFIL == NFILE ) LASTFLAG = .TRUE.
            CALL RDMEDSPD( IDEV, TZONE, TYPNAM, LASTFLAG )

            CLOSE( IDEV )

        END DO     !  End of loop over time steps

C.............  Abort if error found
        IF ( EFLAG ) THEN
            MESG = 'Problem with input file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

1003    MESG = 'ERROR: Could not open file ' // CRLF() // BLANK10//
     &         NAMTMP
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx
94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE GENMEDSOUT


SUBROUTINE RDCEMSUM

!***************************************************************************
!  subroutine body starts at line 86
!
!  DESCRIPTION:
!      This subroutine reads the CEM summary file produced by the utility
!      CEMScan.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!      Subroutines: I/O API subroutine
!
!  REVISION  HISTORY:
!      Created 06/05 by C. Seppanen
!
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!***************************************************************************
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
!.........  This module contains data for day- and hour-specific data
    USE MODDAYHR, ONLY: NOBRLIST, OBRLIST, ANNNOX, ANNSO2, ANNGLOAD,&
    &                    ANNSLOAD, ANNHEAT

    IMPLICIT NONE

!.........  INCLUDES
    INCLUDE 'EMCNST3.EXT'   ! emissions constant parameters

!.........  SUBROUTINE ARGUMENTS

!.........  EXTERNAL FUNCTIONS
    LOGICAL, EXTERNAL :: BLKORCMT

!.........  File names and unit numbers
    INTEGER IDEV    ! input file unit

!.........  Other local variables
    INTEGER            I, N                   ! counters
    INTEGER            IOS                    ! I/O status

    REAL               NOXVAL                 ! tmp NOx value
    REAL               SO2VAL                 ! tmp SO2 value
    REAL               OPTIME                 ! tmp operating time (unused)
    REAL               GLOAD                  ! tmp gross load
    REAL               SLOAD                  ! tmp steam load
    REAL               HTINPUT                ! tmp heat input

    LOGICAL            EFLAG                  ! true: an error occurred

    CHARACTER(ORSLEN3) CORS                   ! tmp ORIS ID
    CHARACTER(BLRLEN3) BLID                   ! tmp boiler ID
    CHARACTER(256)     MESG                   ! message buffer

    CHARACTER(16) ::   PROGNAME = 'RDCEMSUM'  ! program name

!***********************************************************************
!   begin body of program RDCEMSUM

    EFLAG = .FALSE.

!.........  Prompt for input file
    MESG = 'Enter logical name of the CEM SUMMARY file'
    IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'CEMSUM', PROGNAME )

!.........  Loop through file and count number of valid combinations
    NOBRLIST = 0
    DO
        READ( IDEV, *, IOSTAT=IOS ) MESG

!.............  Check for I/O errors
        IF( IOS > 0 ) THEN
            MESG = 'Problem reading CEM summary file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!.............  Check for end of file
        IF( IOS < 0 ) EXIT

!.............  Skip blank lines and comments
        IF( BLKORCMT( MESG ) ) THEN
            CYCLE
        ELSE
            BACKSPACE( IDEV )
        END IF

        READ( IDEV, 93010 ) CORS, BLID, NOXVAL, SO2VAL,&
        &    OPTIME, GLOAD, SLOAD, HTINPUT

        NOBRLIST = NOBRLIST + 1
    END DO

    REWIND( IDEV )

!.........  Exit if error
    IF( EFLAG ) THEN
        MESG = 'Problem reading CEM summary file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Allocate memory to store data
    ALLOCATE( OBRLIST( NOBRLIST ),&
    &           ANNNOX( NOBRLIST ),&
    &           ANNSO2( NOBRLIST ),&
    &         ANNGLOAD( NOBRLIST ),&
    &         ANNSLOAD( NOBRLIST ),&
    &          ANNHEAT( NOBRLIST ), STAT=IOS )
    CALL CHECKMEM( IOS, 'OBRLIST...ANNHEAT', PROGNAME )

!.........  Loop through file and store data
    N = 0
    DO
        READ( IDEV, *, IOSTAT=IOS ) MESG

!.............  Check for end of file
        IF( IOS < 0 ) EXIT

        IF( BLKORCMT( MESG ) ) THEN
            CYCLE
        ELSE
            BACKSPACE( IDEV )
        END IF

        READ( IDEV, 93010 ) CORS, BLID, NOXVAL, SO2VAL,&
        &    OPTIME, GLOAD, SLOAD, HTINPUT

        N = N + 1

        OBRLIST ( N ) = ADJUSTR( CORS ) // ADJUSTR( BLID )
        ANNNOX  ( N ) = NOXVAL
        ANNSO2  ( N ) = SO2VAL
        ANNGLOAD( N ) = GLOAD
        ANNSLOAD( N ) = SLOAD
        ANNHEAT ( N ) = HTINPUT

    END DO

    NOBRLIST = N

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

93010 FORMAT( A6, 1X, A6, 6( 1X, E17.10 ) )

END SUBROUTINE RDCEMSUM

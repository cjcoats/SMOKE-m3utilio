
SUBROUTINE RDSICDSC( FDEV )

!**************************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine allocates memory for and reads the SIC descriptions.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 4/2004 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!**************************************************************************
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
!**************************************************************************
    USE M3UTILIO

!...........   Modules for public variables
!...........   This module contains the lists of unique source characteristics
    USE MODLISTS, ONLY: SICDESC, NINVSIC, INVSIC

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters

!...........   Subroutine arguments
    INTEGER, INTENT (IN) :: FDEV          ! file unit number

!...........   EXTERNAL FUNCTIONS and their descriptions:

    INTEGER, EXTERNAL :: GETFLINE
    LOGICAL, EXTERNAL :: BLKORCMT

!.........  Local parameters
    INTEGER, PARAMETER :: MXSEG = 5       ! max # of potential line segments

!.........  Other arrays
    CHARACTER( SDSLEN3 ) SEGMENT( MXSEG ) ! Segment of parsed lines

!...........   Local variables
    INTEGER         J, L, N               ! indices and counters

    INTEGER         ENDLEN                ! end length for reading descriptn
    INTEGER         IOS                   ! i/o status
    INTEGER      :: IREC = 0              ! record number
    INTEGER      :: NLINES = 0            ! number of lines in input file

    LOGICAL      :: DFLAG = .FALSE.       ! true: processing delimited format
    LOGICAL      :: FFLAG = .FALSE.       ! true: processing fixed format
    LOGICAL      :: EFLAG = .FALSE.       ! true: error found

    CHARACTER(256)  LINE                  ! Read buffer for a line
    CHARACTER(300)  MESG                  ! Message buffer

    CHARACTER(SICLEN3) CSIC               ! tmp SIC

    CHARACTER(16) :: PROGNAME = 'RDSICDSC'    !  program name

!***********************************************************************
!   Begin body of subroutine RDSICDSC

    REWIND( FDEV )  ! In case of multiple calls

!.........  Get the number of lines in the file
    NLINES = GETFLINE( FDEV, 'SIC Descriptions' )

!.........  Allocate memory for the SIC descriptions and initialize
    ALLOCATE( SICDESC( NINVSIC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SICDESC', PROGNAME )
    SICDESC = 'Description unavailable'          ! array

!.........  Read the SIC descriptions, and store with SIC
    ENDLEN = SICLEN3 + SDSLEN3
    DO N = 1, NLINES

        READ ( FDEV, 93000, END=998, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF ( IOS .GT. 0 ) THEN
            WRITE( MESG, 94010)&
            &      'I/O error', IOS, 'reading SIC '//&
            &      'description file at line', IREC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!.............  Check header to define the format (comma-delimited or fixed)
!               If it is not header, Skip blank and comment lines
        IF( INDEX( LINE, '#FIXED' ) > 0 ) THEN
            FFLAG = .TRUE.    ! processing fixed format
            CYCLE
        ELSE IF( INDEX( LINE, '#DELIMITED' ) > 0 ) THEN
            DFLAG = .TRUE.    ! processing delimited format
            CYCLE
        ELSE IF( BLKORCMT( LINE ) ) THEN
            CYCLE
        END IF

!.............  Left adjust line
        LINE = ADJUSTL( LINE )

!.............  Get SIC line
        IF( DFLAG ) THEN
            CALL PARSLINE( LINE, 2, SEGMENT )
            CSIC = SEGMENT( 1 )( 1:SICLEN3 )
            CALL PADZERO( CSIC )

!.................  Find SIC in inventory list, and if it's in the
!                   inventory, store the description.
            J = FINDC( CSIC, NINVSIC, INVSIC )

            IF ( J .GT. 0 ) THEN
                SICDESC( J ) = ADJUSTL( SEGMENT( 2 ) )
            END IF

        ELSE IF( FFLAG ) THEN
            CSIC = LINE( 1:SICLEN3-SICEXPLEN3 )
            CALL PADZERO( CSIC )

!.................  Find SIC in inventory list, and if it's in the
!                   inventory, store the description.
            J = FINDC( CSIC, NINVSIC, INVSIC )

            IF ( J .GT. 0 ) THEN
                SICDESC( J ) = ADJUSTL( LINE( SICLEN3-SICEXPLEN3+1:ENDLEN ) )
            END IF

        ELSE
            MESG = 'ERROR: Missing file format header in the SIC '//&
            &    'description file. Refer to Chapter 8 of the '//&
            &    'SMOKE manual for information on the SICDESC format'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        END IF

    END DO

!.........  Successful completion
    RETURN

!.........  Unexpected end of file
998 MESG = 'INTERNAL ERROR: Unexpected end of SIC description file'
    CALL M3MSG2( MESG )

    CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )


!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDSICDSC


        SUBROUTINE RDMACTDSC( FDEV )

C**************************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine allocates memory for and reads the MACT descriptions.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 7/2005 by B. Baek
C       Version 11/2023 by CJC:  USE M3UTILIO and related changes
C**************************************************************************
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
C**************************************************************************
        USE M3UTILIO

C...........   Modules for public variables
C...........   This module contains the lists of unique source characteristics
        USE MODLISTS, ONLY: MACTDESC, NINVMACT, INVMACT

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   Subroutine arguments
        INTEGER, INTENT (IN) :: FDEV          ! file unit number

C...........   EXTERNAL FUNCTIONS and their descriptions:

        INTEGER, EXTERNAL :: GETFLINE
        LOGICAL, EXTERNAL :: BLKORCMT

C.........  Local parameters
        CHARACTER(16), PARAMETER :: PROGNAME = 'RDMACTDSC'    !  program name
        INTEGER      , PARAMETER :: MXSEG = 5       ! max # of potential line segments

C.........  Other arrays
        CHARACTER( SDSLEN3 ) SEGMENT( MXSEG ) ! Segment of parsed lines

C...........   Local variables
        INTEGER         J, L, N               ! indices and counters

        INTEGER         ENDLEN                ! end length for reading descriptn
        INTEGER         IOS                   ! i/o status
        INTEGER      :: IREC = 0              ! record number
        INTEGER      :: NLINES = 0            ! number of lines in input file

        LOGICAL      :: EFLAG = .FALSE.       ! true: error found

        CHARACTER(256)  LINE                  ! Read buffer for a line
        CHARACTER(300)  MESG                  ! Message buffer

        CHARACTER( MACLEN3 ) TMACT            ! tmp MACT

C***********************************************************************
C   Begin body of subroutine RDMACTDSC

        REWIND( FDEV )  ! In case of multiple calls

C.........  Get the number of lines in the file
        NLINES = GETFLINE( FDEV, 'MACT Descriptions' )

C.........  Allocate memory for the MACT descriptions and initialize
        ALLOCATE( MACTDESC( NINVMACT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MACTDESC', PROGNAME )
        MACTDESC = 'Description unavailable'          ! array

C.........  Read the MACT descriptions, and store with MACT
        ENDLEN = MACLEN3 + SDSLEN3

        DO N = 1, NLINES

            READ ( FDEV, '( A )', END=998, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010)
     &                'I/O error', IOS, 'reading MACT '//
     &                'description file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

C.............  Skip blank and comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

C.............  Left adjust line
            LINE = ADJUSTL( LINE )

C.............  Get MACT line
            CALL PARSLINE( LINE, 2, SEGMENT )
            TMACT = SEGMENT( 1 )( 1:MACLEN3 )

            CALL PADZERO( TMACT )

C.............  Find MACT in inventory list, and if it's in the inventory,
C               store the description.
            J = FINDC( TMACT, NINVMACT, INVMACT )

            IF ( J .GT. 0 ) THEN
                MACTDESC( J ) = ADJUSTL( SEGMENT( 2 ) )
            END IF

        END DO

C.........  Successful completion
        RETURN

C.........  Unexpected end of file
998     MESG = 'INTERNAL ERROR: Unexpected end of MACT description file'
        CALL M3MSG2( MESG )

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE RDMACTDSC

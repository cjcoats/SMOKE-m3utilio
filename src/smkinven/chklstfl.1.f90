
SUBROUTINE CHKLSTFL( NLINE, FNAME, NLSTSTR, FILFMT )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine checks the format of the list-formatted inventory
!      file and returns the code for the type of files it contains.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 10/98 by M. Houyoux
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
!***************************************************************************
    USE M3UTILIO

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN) :: NLINE            ! number of lines in file
    CHARACTER(*), INTENT (IN) :: FNAME            ! logical name of file
    CHARACTER(*), INTENT (IN) :: NLSTSTR( NLINE ) ! contents of file by line
    INTEGER     , INTENT(OUT) :: FILFMT ( NLINE ) ! file format code

!...........   EXTERNAL FUNCTIONS and their descriptions:

    INTEGER, EXTERNAL :: GETFORMT
    INTEGER, EXTERNAL :: GETINVYR

!...........   File units and logical/physical names
    INTEGER      TDEV        !  emissions file in list format file

!...........   Other local variables
    INTEGER      I, J
    INTEGER      EXTFORMAT           !  used when a #LIST entry has format info
    INTEGER      IOS                 !  I/O status

    LOGICAL   :: EMSFLAG   = .FALSE. !  true: at least one file is EMS format
    LOGICAL   :: ORLFLAG  = .FALSE.  !  true: at least one file is ORL format

    CHARACTER(300)  INFILE      !  input file line buffer
    CHARACTER(500)  MESG        !  message buffer

    CHARACTER(16) :: PROGNAME =  'CHKLSTFL' ! program name

!***********************************************************************
!   begin body of subroutine CHKLSTFL

    EMSFLAG  = .FALSE.   ! Need to reset for each each subroutine call
    ORLFLAG  = .FALSE.
    EXTFORMAT = -1

!.........  Loop through lines of list-formatted file to check the formats
    DO J = 1, NLINE

!.............  Skip blank lines
        IF( NLSTSTR( J ) == ' ' ) CYCLE

!.............  Store the current line's file name
        INFILE = NLSTSTR( J )

!.............  Skip INVYEAR packet
        I = GETINVYR( INFILE )
        IF( I .GT. 0 ) THEN
            FILFMT( J ) = -1
            CYCLE
        END IF

!.............  Skip the date range packet
        I = INDEX( INFILE, 'DATERANGE' )
        IF( I .GT. 0 ) THEN
            FILFMT( J ) = -1
            CYCLE
        END IF

!.............  Check for #LIST entry
        I = INDEX( INFILE, '#LIST' )
        IF( I .GT. 0 ) THEN
            FILFMT( J ) = -1

            IF( INDEX( INFILE, 'EMS-95' ) > 0 ) THEN
                EXTFORMAT = EMSFMT

            ELSE IF( INDEX( INFILE, 'CEM' ) > 0 ) THEN
                EXTFORMAT = CEMFMT

            ELSE IF( INDEX( INFILE, 'FF10' ) > 0 ) THEN
                EXTFORMAT = FF10FMT

            ELSE IF( INDEX( INFILE, 'MEDS' ) > 0 ) THEN
                EXTFORMAT = MEDSFMT

            ELSE IF( INDEX( INFILE, 'GRID' ) > 0 ) THEN
                EXTFORMAT = NCDFMT

            ELSE IF( INDEX( INFILE, 'ORL' ) > 0 ) THEN

                IF( INDEX( INFILE, 'NONPOINT' ) > 0 ) THEN
                    EXTFORMAT = ORLNPFMT

                ELSE IF( INDEX( INFILE, 'FIRE' ) > 0 ) THEN
                    EXTFORMAT = ORLFIREFMT

                ELSE IF( INDEX( INFILE, 'FIREEMIS' ) > 0 ) THEN
                    EXTFORMAT = ORLDYFRFMT

                ELSE
                    EXTFORMAT = ORLFMT
                    WRITE( MESG,94010 ) 'WARNING: No format '//&
                    &  'included with #LIST ORL entry in '//&
                    &  'inventory input file.'// CRLF() // BLANK5 //&
                    &  'Assuming ORL point, nonroad, or on-road.'

                    CALL M3MSG2( MESG )

                END IF

            END IF

            CYCLE
        END IF

!.............  Open INFILE
        TDEV = JUNIT()

        OPEN( TDEV, FILE=INFILE, STATUS='OLD', IOSTAT=IOS )

!.............  Check for problems opening raw input file
        IF( IOS /= 0 ) THEN
            WRITE( MESG,94010 ) 'Problem at line ', J, 'of ' //&
            &   TRIM( FNAME ) // '.' // ' Could not open file:' //&
            &   CRLF() // BLANK5 // TRIM( INFILE )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!.............  Determine format of INFILE
        FILFMT( J ) = GETFORMT( TDEV, EXTFORMAT )
        CLOSE( TDEV )

!.............  Set flag based on format
        IF( FILFMT( J ) == EMSFMT ) EMSFLAG = .TRUE.
        IF( FILFMT( J ) == FF10FMT    .OR.&
        &    FILFMT( J ) == MEDSFMT    .OR.&
        &    FILFMT( J ) == ORLFMT     .OR.&
        &    FILFMT( J ) == ORLNPFMT   .OR.&
        &    FILFMT( J ) == ORLFIREFMT .OR.&
        &    FILFMT( J ) == ORLDYFRFMT     ) ORLFLAG = .TRUE.

!.............  Check that file formats are consistent
        IF( EMSFLAG .AND. FILFMT( J ) /= EMSFMT ) THEN
            WRITE( MESG,94010 )&
            &       'ERROR: In SMOKE list-formatted inventory file, '&
            &       // TRIM( FNAME ) // ', at least one file is ' //&
            &       CRLF() // BLANK10 // 'EMS-95 format ' //&
            &       'while another is not. When using an EMS-95 ' //&
            &       'formatted' // CRLF() // BLANK10 //&
            &       'inventory, all other inventory files must ' //&
            &       'also be in EMS-95 format.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF( ORLFLAG                   .AND.&
        &    FILFMT( J ) /= FF10FMT    .AND.&
        &    FILFMT( J ) /= MEDSFMT    .AND.&
        &    FILFMT( J ) /= ORLFMT     .AND.&
        &    FILFMT( J ) /= ORLNPFMT   .AND.&
        &    FILFMT( J ) /= ORLFIREFMT .AND.&
        &    FILFMT( J ) /= ORLDYFRFMT      ) THEN
            WRITE( MESG,94010 )&
            &       'ERROR: In SMOKE list-formatted inventory file, '&
            &       // TRIM( FNAME ) // ', at least one file is ' //&
            &       CRLF() // BLANK10 // 'ORL ' //&
            &       'format while another is not ' //&
            &       'ORL format.' // CRLF() // BLANK10 //&
            &       'When using ORL formatted ' //&
            &       'inventories, all other inventories ' //&
            &       CRLF() // BLANK10 // 'must also be ORL '&
            &       // 'format.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        EXTFORMAT = FILFMT( J )

    END DO     ! End of loop through list-formatted file

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE CHKLSTFL


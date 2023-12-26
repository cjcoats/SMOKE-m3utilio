
INTEGER FUNCTION GETFORMT( FDEV, EXTFORMAT )

    !***********************************************************************
    !  function body starts at line 109
    !
    !  DESCRIPTION:
    !      This function returns the format of the inventory file, assuming that
    !      the source category is already known.
    !
    !  PRECONDITIONS REQUIRED:
    !      Inventory file opened on unit FDEV
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !      Subroutines: I/O API subroutine
    !      Functions: I/O API functions, CHKINT, CHKREAL
    !
    !  REVISION  HISTORY:
    !       Created 11/98 by M. Houyoux
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

    !.......   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   EXTERNAL FUNCTIONS:

    !.......   SUBROUTINE ARGUMENTS
    INTEGER, INTENT(IN) :: FDEV           ! unit number of input point sources file
    INTEGER, INTENT(IN) :: EXTFORMAT      ! external format (used for list files)

    !.......   Other local variables
    INTEGER         L, L1, L2, L3       ! indices
    INTEGER         IREC        ! current record number
    INTEGER         IOS         ! I/O status

    CHARACTER(256)   MESG            ! message buffer
    CHARACTER(100)   LINE            ! buffer to read line into

    CHARACTER(16), PARAMETER :: PROGNAME = 'GETFORMT'      ! program name

    !***********************************************************************
    !   begin body of function GETFORMT

    GETFORMT = IMISS3

    !.....  Loop through lines of file
    IREC = 0
    DO

        READ( FDEV, 93000, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        !.....  Check for I/O errors
        IF( IOS > 0 ) THEN
            WRITE( MESG, 94010 )&
               'I/O error', IOS, 'reading inventory file at line', IREC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        !.....  Check if we've reached the end of the file
        IF( IOS < 0 ) EXIT

        !.....  Skip blank lines
        IF( LINE == ' ' ) CYCLE

        !.....  Check for header lines
        L1 = INDEX( LINE, CINVHDR )
        IF( L1 .EQ. 1 ) THEN
            !.....  Check if format is provided as a header entry
            L = INDEX( LINE, 'LIST' )
            IF( L .GT. 0 ) THEN
                GETFORMT = LSTFMT
                EXIT             ! To end of read loop
            END IF

            L = INDEX( LINE, 'EMS-95' )
            IF( L .GT. 0 ) THEN
                GETFORMT = EMSFMT
                EXIT             ! To end read loop
            END IF

            L = INDEX( LINE, 'CEM' )
            IF( L .GT. 0 ) THEN
                GETFORMT = CEMFMT
                EXIT             ! To end read loop
            END IF

            L = INDEX( LINE, 'FF10' )
            IF( L .GT. 0 ) THEN
                GETFORMT = FF10FMT
                EXIT             ! To end read loop
            END IF

            L = INDEX( LINE, 'MEDS' )
            IF( L .GT. 0 ) THEN
                GETFORMT = MEDSFMT
                EXIT             ! To end read loop
            END IF

            L = INDEX( LINE, 'GRID' )
            IF( L .GT. 0 ) THEN
                GETFORMT = NCDFMT
                EXIT             ! To end read loop
            END IF

            L = INDEX( LINE, 'ORL' )
            IF( L .GT. 0 ) THEN
                L = INDEX( LINE, 'NONPOINT' )
                IF( L .GT. 0 ) THEN
                    GETFORMT = ORLNPFMT
                ELSE

                    L2 = INDEX( LINE, 'FIRE' )
                    L3 = INDEX( LINE, 'FIREEMIS' )
                    IF( L2 .GT. 0 .AND. L3 .LE. 0 ) THEN
                        GETFORMT = ORLFIREFMT
                    ELSE IF( L3 .GT. 0 .AND. L2 .GT. 0 ) THEN
                        GETFORMT = ORLDYFRFMT
                    ELSE
                        GETFORMT = ORLFMT
                    END IF

                END IF
                EXIT               ! To end read loop

            END IF

            !.....  Otherwise, this is not a blank line, but also not a
            !               header line, so we're into the main body of the inventory
        ELSE
            EXIT

        END IF

    END DO         ! To head of read loop

    !.....  Rewind inventory file
    REWIND( FDEV )

    !.....  Check if an external format was passed in (used for list files)
    IF( GETFORMT == IMISS3 ) THEN
        IF( EXTFORMAT /= -1 ) THEN
            GETFORMT = EXTFORMAT
        END IF
    END IF

    !.....  If format has not been set, print error about missing header
    IF ( GETFORMT .EQ. IMISS3 ) THEN

        MESG = 'ERROR: Could not determine inventory file ' //  &
               'format due to missing or bad header ' //        &
               'information. ' // CRLF() // BLANK16 //          &
               'Valid headers are: ' // CRLF() // BLANK16 //    &
               '#LIST, #EMS-95, #FF10, #MEDS, #GRID, #ORL, or #CEM'
        CALL M3MSG2( MESG )
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Formatted file I/O formats.... 93xxx

93000 FORMAT( A )

    !.......   Internal buffering formats.... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END FUNCTION GETFORMT



SUBROUTINE GETNTISZ( FDEV, CATEGORY, NLINES )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine returns an exact number of records for a raw toxics
    !      inventory input file opened on unit FDEV.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created by M. Houyoux 1/99
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

    !.......   MODULES for public variables
    !.....  This module contains the lists of unique inventory information
    USE MODLISTS, ONLY: NUNIQCAS, UNIQCAS, UCASNKEP

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constat parameters

    !.......   SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN ) :: FDEV         !  unit number of input file
    CHARACTER(*), INTENT (IN ) :: CATEGORY     !  description of source category
    INTEGER     , INTENT (OUT) :: NLINES       !  total number of inventory records

    !.......   Local allocatable arrays
    CHARACTER(50)   SEGMENT( 10 )     ! pieces of input line

    !.......   Other local variables
    INTEGER         I                ! counter
    INTEGER         IOS              ! i/o status
    INTEGER      :: IREC    = 0      ! input line counter
    INTEGER         CASPOS           ! position of CAS number in line
    INTEGER         NSEG             ! number of segments in line

    CHARACTER(300)     LINE             ! input file line buffer
    CHARACTER(CASLEN3) CASNUM           ! CAS number from input line
    CHARACTER(300)     MESG             ! message buffer

    CHARACTER(16), PARAMETER ::   PROGNAME = 'GETNTISZ'     ! program name

    !***********************************************************************
    !   begin body of subroutine GETNTISZ

    !.....  Write message to screen
    CALL M3MSG2( 'Determining size of toxics inventory...' )

    !.....  Initialize counters
    NLINES = 0
    IREC = 0

    !.....  Allocate memory for line segments
    SELECT CASE( CATEGORY )
      CASE( 'AREA' )
        NSEG = 9
        CASPOS = 4

      CASE( 'MOBILE' )
        NSEG = 6
        CASPOS = 4

    END SELECT
    SEGMENT = ' '      ! array

    !.....  Loop through lines in file
    DO

        !.....  Read line of file
        READ( FDEV, 93000, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        !.....  Check I/O error status
        IF( IOS > 0 ) THEN
            WRITE( MESG, 94010 ) 'Error', IOS,  'reading ' //    &
                   'inventory file as character strings at line', IREC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSEIF( IOS < 0 ) THEN          ! reached end of file
            EXIT

        END IF

        !.....  Skip blank lines
        IF( LINE == ' ' ) CYCLE

        !.....  Skip header lines
        IF( LINE( 1:1 ) == CINVHDR ) CYCLE

        !.....  Extract CAS number from line
        CALL PARSLINE( LINE, NSEG, SEGMENT )
        CASNUM = ADJUSTL( SEGMENT( CASPOS ) )

        !.....  Find CAS number in unique sorted list
        I = FINDC( CASNUM, NUNIQCAS, UNIQCAS )
        IF( I < 1 ) CYCLE          ! skip entries with invalid CAS numbers

        !.....  Add number of records for this CAS
        NLINES = NLINES + UCASNKEP( I )

    END DO

    IF( NLINES .EQ. 0 ) THEN
        MESG = 'Inventory file has no valid lines of inventory data.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    REWIND( FDEV )

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Formatted file I/O formats.... 93xxx

93000 FORMAT( A )

    !.......   Internal buffering formats.... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE GETNTISZ

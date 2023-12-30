
SUBROUTINE ERRPKTS( PKTTYP, JT, JX, SKIPPOL, JTMAX, JXMAX, EFLAG )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      Check dimensions with maximum dimensions for packets, and reset
    !      maximum to actual in case any records were skipped.  Also, produce
    !      messages if pollutant was skipped or if the number of usable packets
    !      was zero.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 3/99 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
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

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS:

    CHARACTER(*), INTENT (IN) :: PKTTYP     ! packet type
    INTEGER     , INTENT (IN) :: JT         ! idx to control data tables
    INTEGER     , INTENT (IN) :: JX         ! idx to ungrouped cntl x-ref tables
    LOGICAL     , INTENT (IN) :: SKIPPOL    ! true: skipped rec is pol-spcfic
    INTEGER  , INTENT(IN OUT) :: JTMAX      ! max allowed JT
    INTEGER  , INTENT(IN OUT) :: JXMAX      ! max allowed JX
    LOGICAL     , INTENT(OUT) :: EFLAG      ! true: error occurred

    !.......   Other local variables
    INTEGER         L            !  counters and indices

    CHARACTER(256)  MESG         ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'ERRPKTS'     ! program name

    !***********************************************************************
    !   begin body of subroutine ERRPKTS

    !.......  Get length of packet name for messages
    L = LEN_TRIM( PKTTYP )

    !.......  Write warning message for pollutants in file that are
    !           not in master list
    IF( SKIPPOL ) THEN
        MESG = 'Pollutant-specific entries in the ' // PKTTYP // ' packet' // &
               CRLF()//BLANK10 // 'have been skipped.'
        CALL M3WARN( PROGNAME, 0, 0, MESG )
    END IF

    !.......  Error for overflow of cross-reference information
    IF( JX .GT. JXMAX ) THEN

        EFLAG = .TRUE.
        WRITE( MESG,94010 ) 'INTERNAL ERROR: dimension for ' // &
               'storing cross-reference for the ' //            &
               PKTTYP // CRLF() // BLANK10 //                   &
               'packet was', JXMAX, 'but actually needed', JX
        CALL M3MSG2( MESG )

    !.......  Otherwise, store final count
    ELSE
        JXMAX = JX

    END IF

    !.......  Error for overflow of control table information
    IF( JT .EQ. 0 ) THEN
        EFLAG = .TRUE.
        MESG = 'No usable ' // TRIM( PKTTYP ) // ' control packet entries'
        CALL M3MSG2( MESG )

    ELSE IF( JT .GT. JTMAX ) THEN

        EFLAG = .TRUE.
        WRITE( MESG,94010 ) 'INTERNAL ERROR: dimension for ' // &
               'storing control data for the ' //               &
               PKTTYP // CRLF() // BLANK10 //                   &
               'packet was', JTMAX, 'but actually needed', JT
        CALL M3MSG2( MESG )

    !.......  Otherwise, store final count
    ELSE
        JTMAX = JT

    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE ERRPKTS

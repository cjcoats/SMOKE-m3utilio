
SUBROUTINE WCNTLREP( CDEV, GDEV, LDEV, MDEV )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine write the report file for emission controls by source
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version ??/???? by ???
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***************************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !                System
    ! File: %W%
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
    ! Pathname: %P%
    ! Last updated: %G% %U%
    !
    !***************************************************************************
    USE M3UTILIO

    !.......  MODULES for public variables
    !.......  This module contains the inventory arrays
    USE MODSOURC, ONLY: CSOURC

    !.......  This module contains the control packet data and control matrices
    USE MODCNTRL, ONLY: NVCMULT, PNAMMULT, RPTDEV, PCTLFLAG

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, CRL, NSRC, NCHARS

    IMPLICIT NONE

    !.......   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS

    INTEGER     , INTENT (IN) :: CDEV       ! file unit no. for tmp CTL file
    INTEGER     , INTENT (IN) :: GDEV       ! file unit no. for tmp CTG file
    INTEGER     , INTENT (IN) :: LDEV       ! file unit no. for tmp ALW file
    INTEGER     , INTENT (IN) :: MDEV       ! file unit no. for tmp MACT file

    !.......  Local arrays
    INTEGER             OUTTYPES( NVCMULT,6 )     ! var type:int/real

    CHARACTER(NAMLEN3)  OUTNAMES( NVCMULT,6 )     ! var names
    CHARACTER(IOULEN3)  OUTUNITS( NVCMULT,6 )     ! var units
    CHARACTER(MXDLEN3)  OUTDESCS( NVCMULT,6 )     ! var descriptions

    !.......   Other local variables
    INTEGER          S, V      ! counters and indices

    INTEGER          CIDX       ! control plant index
    INTEGER          ODEV       ! file unit no. for output report

    REAL             E_IN       ! emissions before controls
    REAL             E_OUT      ! emissions after controls
    REAL             FAC        ! control factor

    CHARACTER(256)     MESG           ! message buffer
    CHARACTER(NAMLEN3) PNAM           ! tmp pollutant name

    CHARACTER(16), PARAMETER :: PROGNAME = 'WCNTLREP'     ! program name

    !***********************************************************************
    !   begin body of subroutine WCNTLREP

    !.......  Rewind temporary files
    IF( CDEV .GT. 0 ) REWIND( CDEV )
    IF( GDEV .GT. 0 ) REWIND( GDEV )
    IF( LDEV .GT. 0 ) REWIND( LDEV )
    IF( MDEV .GT. 0 ) REWIND( MDEV )

    !.......  Open reports file
    IF( MAX( CDEV, GDEV, LDEV, MDEV ) .GT. 0 ) THEN
        RPTDEV( 2 ) = PROMPTFFILE( 'Enter logical name for SUMMARY CONTROLS REPORT',&
                                   .FALSE., .TRUE., CRL // 'CSUMREP', PROGNAME )
        ODEV = RPTDEV( 2 )
    END IF

    !.......  For each pollutant that receives controls, obtain variable
    !             names for control efficiency, rule effectiveness, and, in the
    !             case of AREA sources, rule penetration. These variable names
    !             will be used in reading the inventory file.

    !.......  Check that NVCMULT does not equal 0, otherwise some systems will get confused
    IF( NVCMULT == 0 ) RETURN

    CALL BLDENAMS( CATEGORY, NVCMULT, 6, PNAMMULT, OUTNAMES, OUTUNITS, OUTTYPES, OUTDESCS )

    !.......  Read in indices from temporary files. No error checking is
    !             performed because it is assumed that the program has already
    !             successfully written the temporary files.

    !.......  Loop through pollutants
    DO V = 1, NVCMULT

    !.......  Loop through sources and output
        DO S = 1, NSRC

    !.......  If MACT packet applies for this pollutant
            IF( PCTLFLAG( V, 4 ) ) THEN

                READ( MDEV,* ) CIDX, PNAM, E_IN, E_OUT, FAC

                CALL CONTROL_MESG( ODEV, 'MACT', S, CIDX, PNAM, E_IN, E_OUT, FAC )

            END IF

    !.......  If CONTROL packet applies for this pollutant
            IF( PCTLFLAG( V, 1 ) ) THEN

                READ( CDEV,* ) CIDX, PNAM, E_IN, E_OUT, FAC

                CALL CONTROL_MESG( ODEV, 'CONTROL', S, CIDX, PNAM, E_IN, E_OUT, FAC )

            END IF

    !.......  If CTG packet applies for this pollutant
            IF( PCTLFLAG( V, 2 ) ) THEN
                READ( GDEV,* ) CIDX, PNAM, E_IN, E_OUT, FAC

                CALL CONTROL_MESG( ODEV, 'CTG', S, CIDX, PNAM, E_IN, E_OUT, FAC )

            END IF

    !.......  If ALLOWABLE packet applies for this pollutant
            IF( PCTLFLAG( V, 3 ) ) THEN

                READ( LDEV,* ) CIDX, PNAM, E_IN, E_OUT, FAC

    !.......  Interpret control information and write message
                CALL CONTROL_MESG( ODEV, 'ALLOWABLE', S, CIDX, PNAM, E_IN, E_OUT, FAC )

            END IF

    !.......  If REACTIVITY packet applies for this pollutant
    ! note: must be added

        END DO
    END DO

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

    !******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

    !.......  This internal subprogram
    SUBROUTINE CONTROL_MESG( FDEV, CTYPE, SMKID, CIDX, PNAM,&
                             EMIS_IN, EMIS_OUT, FACTOR )

        !.......  Subprogram arguments
        INTEGER     , INTENT (IN) :: FDEV             ! output unit number
        CHARACTER(*), INTENT (IN) :: CTYPE            ! control type
        INTEGER     , INTENT (IN) :: SMKID            ! smoke ID
        INTEGER     , INTENT (IN) :: CIDX             ! control index
        CHARACTER(*), INTENT (IN) :: PNAM             ! pollutant name
        REAL        , INTENT (IN) :: EMIS_IN          ! emissions before control
        REAL        , INTENT (IN) :: EMIS_OUT         ! emissions after control
        REAL        , INTENT (IN) :: FACTOR           ! EMIS_IN*FACTOR = EMIS_OUT

        !.......  Local variables
        INTEGER           L, L2

        INTEGER, SAVE ::  PSMKID = 0           ! previous smoke ID

        CHARACTER(300)    BUFFER                         ! message buffer
        CHARACTER(300)    MESG                           ! message buffer

        CHARACTER(NAMLEN3), SAVE :: LNAM          ! previous pollutant name

        !----------------------------------------------------------------------

        !.......  Skip records that are not controlled
        IF( CIDX .EQ. 0 ) RETURN

        !.......  Write message for pollutant for each new pollutant
        IF( PNAM .NE. LNAM ) THEN

            MESG = 'Source controls for pollutant ' // PNAM
            WRITE( FDEV,93000 ) ' '
            WRITE( FDEV, 93000 ) REPEAT( '-', 80 )
            WRITE( FDEV,93000 ) TRIM( MESG )
            WRITE( FDEV,93000 ) ' '

        END IF

        !.......  Only write out header line if source is different from previous
        !               or if the message for a new pollutant was written
        IF( SMKID .NE. PSMKID .OR. PNAM .NE. LNAM ) THEN

            CALL FMTCSRC( CSOURC( S ), NCHARS, BUFFER, L2 )
            WRITE( FDEV,93000 ) BLANK5 // BUFFER( 1:L2 )
            PSMKID = SMKID

        END IF

        LNAM = PNAM

        !.......  Write warning if source is controled and factor is 1.0
        IF( FACTOR .EQ. 1. ) THEN
            WRITE( MESG,93380 ) CTYPE, EMIS_IN, EMIS_OUT
            WRITE( FDEV,93000 ) TRIM( MESG )

        !.......  Otherwise, write standard control message
        ELSE
            WRITE( MESG,93400 ) CTYPE, EMIS_IN, EMIS_OUT, FACTOR
            WRITE( FDEV,93000 ) TRIM( MESG )

        END IF

        RETURN

    !*************  SUBPROGRAM FORMAT  STATEMENTS   **************************

    !.......   Formatted file I/O formats...... 93xxx

93000   FORMAT( A )

93380   FORMAT( 10X, A, ' Packet. Before: ', E11.4, ' After: ',         &
                E11.4, ' [tons/day]. WARNING: Control factor of 1.' )

93400   FORMAT( 10X, A, ' Packet. Before: ', E11.4, ' After: ',         &
                E11.4, ' [tons/day]. Factor:', E9.3 )

    END SUBROUTINE CONTROL_MESG

    !----------------------------------------------------------------------

END SUBROUTINE WCNTLREP

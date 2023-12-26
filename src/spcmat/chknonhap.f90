
SUBROUTINE CHKNONHAP( PNAM, EFLAG )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !    This subroutine compares the definition of a NONHAP* pollutant
    !    with the definition given in the GSPRO file.  The routine
    !    is called once per pollutant.  If the pollutant name is not
    !    a NONHAP pollutant or it is the same as the previous pollutant
    !    name, them the routine will exit immediately.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !    Version ??/???? by ???
    !    Version 11/2023 by CJC:  USE M3UTILIO and related changes
    !****************************************************************************/
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !    System
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

    !.......  MODULES for public variables
    !.......  This module contains the lists of unique source characteristics
    USE MODLISTS, ONLY: MXIDAT, ITNAMA, INVDVTS, INVDNAM,  ITKEEPA, ITVTSA

    !.......  This module contains the speciation profiles
    USE MODSPRO, ONLY: NSPDEF, NSPLST, SPCDEFPOL, SPCDEFLST

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   Subroutine arguments

    CHARACTER(*), INTENT  (IN) :: PNAM           ! pol name of interest
    LOGICAL     , INTENT (OUT) :: EFLAG          ! error flag

    !.......   Other local variables

    INTEGER         I, J, K      ! counters and indices
    INTEGER         IOS          ! allocate status
    INTEGER         CNT

    CHARACTER(256)   MESG                  ! message buffer

    CHARACTER(NAMLEN3)       :: PBUF      = ' '       ! tmp pollutant name
    CHARACTER(NAMLEN3), SAVE :: PREV_PNAM = ' '       ! PNAM from last call

    CHARACTER(16), PARAMETER :: PROGNAME = 'CHKNONHAP'     ! program name

    !***********************************************************************
    !   Begin body of subroutine CHKNONHAP

    J = INDEX( PNAM, 'NONHAP' )

    !.......  If current pollutant is not a NONHAP* pollutant, exit
    IF( J .LE. 0 .OR. PNAM .EQ. PREV_PNAM ) RETURN

    !.......  Store pollutant name to check in next iteration
    PREV_PNAM = PNAM

    !.......  Search for pollutant in list of available definitions
    I = INDEX1( PNAM, NSPDEF, SPCDEFPOL )

    !......  When pollutant found as having a definition in the
    !    speciation profile header, check it
    IF( I .GT. 0 ) THEN

        !......  Loop through the INVTABLE pollutants and ensure
        !    that all NONHAP contributors are also defined as part
        !    of this NONHAP* variable in the GSPRO file.
        CNT = 0
        DO J = 1, MXIDAT
            IF( INVDVTS( J ) .NE. 'N' ) THEN

                !.......  Search for data name in GSPRO definition
                K = INDEX1( INVDNAM(J), NSPLST(I), SPCDEFLST(1,I) )

                !.......  If data name not found in definition, then
                !    give an error.  If it is found count it.
                IF( K .LE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )    &
                           'ERROR: Pollutant "'// TRIM( INVDNAM(J) )    &
                           // '" part of definition of '//    &
                           TRIM( PNAM ) // ' in INVTABLE but' //    &
                           CRLF() // BLANK10 // 'is missing from '//    &
                           'the definition in the GSPRO file.'
                    CALL M3MSG2( MESG )
                ELSE
                    CNT = CNT + 1
                END IF

            END IF
        END DO

        !.......  If the count is different, report the pollutants
        !    that are in the GSPRO def'n, but not the inven's
        IF( CNT .NE. NSPLST(I) ) THEN

            DO K = 1, NSPLST(I)
                PBUF = SPCDEFLST(K,I)
                J = INDEX1( PBUF, MXIDAT, ITNAMA )
                IF( J .GT. 0 ) THEN
                    IF( ITKEEPA( J ) .AND.    &
                        ITVTSA ( J ) .EQ. 'N' ) THEN
                        EFLAG = .TRUE.
                        MESG = 'ERROR: Pollutant "'// TRIM( PBUF )//    &
                           '" part of ' // TRIM( PNAM ) //    &
                           ' definition in GSPRO but '// CRLF()//    &
                           BLANK10 // 'not in INVTABLE.'
                        CALL M3MSG2( MESG )
                    END IF
                END IF
            END DO

        END IF

    !......  If no definition found, give a warning
    ELSE
        MESG = 'WARNING: No definition found for ' // TRIM( PNAM )    &
               // ' in the GSPRO file.  Spcmat will' //    &
               CRLF() // BLANK10 // 'not be able to check for '//    &
               'consistency between the inventory and the GSPRO '//    &
               'file.'
        CALL M3MSG2( MESG )

    END IF

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE CHKNONHAP

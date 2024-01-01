
REAL FUNCTION UNITFAC( UNIT1, UNIT2, LNUM )

    !***********************************************************************
    !  function body starts at line
    !
    !  DESCRIPTION:
    !      This function returns a units conversion factor based on the numerator
    !      (if any) of argument 1, and the numerator (if any) of argument two.
    !      If there are no units quotients, it treats the units field like a
    !      numerator.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created by M. Houyoux 10/99
    !       Version 11/2023 by CJC:  USE M3UTILIO, ".f90" source format, and 
    !       related changes
    !**************************************************************************
    !
    ! Project Title: EDSS Tools Library
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
    ! Last updated: %G
    !
    !***************************************************************************

    USE M3UTILIO

    IMPLICIT NONE

    !.......   INCLUDES

    INCLUDE 'EMCNST3.h90'

    !.......   SUBROUTINE ARGUMENTS

    CHARACTER(*), INTENT (IN) :: UNIT1    ! first unit
    CHARACTER(*), INTENT (IN) :: UNIT2    ! second unit
    LOGICAL     , INTENT (IN) :: LNUM     ! true: for numerators

    !.......   Other local variables
    INTEGER         I, K1, K2, L1, L2, LB1, LB2
    INTEGER         EXP1       ! exponent on first unit
    INTEGER         EXP2       ! exponent on second unit

    REAL            PRECNV     ! Conversion for exponential prefixes

    CHARACTER(NAMLEN3) :: BUF1
    CHARACTER(NAMLEN3) :: BUF2

    CHARACTER(16), PARAMETER :: PROGNAME = 'UNITFAC' ! program name

    !***********************************************************************
    !   begin body of function UNITFAC

    !.......  Initialize factor
    UNITFAC = 1.0

    !.......  Retrieve length of unit strings
    L1 = LEN_TRIM( UNIT1 )
    L2 = LEN_TRIM( UNIT2 )

    !.......  Cases where one or more of the provided units are blank
    IF( L1 .LE. 0 .OR. L2 .LE. 0 ) RETURN

    !.......  Determine positions of the divide-by symbol in the units
    K1 = INDEX( UNIT1, '/' )
    K2 = INDEX( UNIT2, '/' )

    !.......  Define buffers as numerators or denominators, depending on
    !         value of LNUM
    IF( LNUM ) THEN
        IF( K1 .LE. 0 ) K1 = L1 + 1
        IF( K2 .LE. 0 ) K2 = L2 + 1
        BUF1 = UNIT1( 1:K1-1 )
        BUF2 = UNIT2( 1:K2-1 )

    ELSE
        IF( K1 .LE. 0 ) K1 = 0
        IF( K2 .LE. 0 ) K2 = 0
        BUF1 = UNIT1( K1+1:L1 )
        BUF2 = UNIT2( K2+1:L2 )

    END IF

    LB1  = LEN_TRIM( BUF1 )
    LB2  = LEN_TRIM( BUF2 )

    CALL UNITMATCH( BUF1 )
    CALL UNITMATCH( BUF2 )

    !.......  Separate out any leading adjustments (e.g., 10E6)
    EXP1 = 0
    EXP2 = 0
    K1 = INDEX( BUF1, '10E' )
    K2 = INDEX( BUF2, '10E' )

    IF( K1 .GT. 0 ) THEN
        I = INDEX( BUF1, ' ' )
        IF( I .GT. 0 ) THEN
            EXP1 = STR2INT( BUF1( K1+3:I ) )
            BUF1 = BUF1( I+1:LB1 )
        END IF
    END IF

    IF( K2 .GT. 0 ) THEN
        I = INDEX( BUF2, ' ' )
        IF( I .GT. 0 ) THEN
            EXP2 = STR2INT( BUF2( K2+3:I ) )
            BUF1 = BUF2( I+1:LB2 )
        END IF
    END IF

    !.......  Set conversion factor by comparing units prefixes
    PRECNV = 10**( EXP1 - EXP2 )

    !.......  Set conversion factor by comparing numerators of main units
    SELECT CASE( BUF1 )

      CASE( 'moles' )

        SELECT CASE( BUF2 )
          CASE( 'moles' )
            UNITFAC = 1.
          CASE DEFAULT
            CALL CASE_NOT_FOUND( UNITFAC )
        END SELECT

      CASE( 'g' )

        SELECT CASE( BUF2 )
          CASE( 'kg' )
            UNITFAC = 1. / 1000.
          CASE( 'tons' )
            UNITFAC = GM2TON
          CASE( 'g' )
          CASE DEFAULT
            CALL CASE_NOT_FOUND( UNITFAC )
        END SELECT

      CASE( 'kg' )

        SELECT CASE( BUF2 )
          CASE( 'g' )
            UNITFAC = 1000.
          CASE( 'tons' )
            UNITFAC = GM2TON / 1000.
          CASE( 'kg' )
          CASE DEFAULT
            CALL CASE_NOT_FOUND( UNITFAC )
        END SELECT

      CASE( 'tons' )

        SELECT CASE( BUF2 )
          CASE( 'g' )
            UNITFAC = TON2GM
          CASE( 'kg' )
            UNITFAC = TON2GM / 1000.
          CASE( 'tons' )
          CASE DEFAULT
            CALL CASE_NOT_FOUND( UNITFAC )
        END SELECT

      CASE( 'yr' )

        SELECT CASE( BUF2 )
          CASE( 'yr' )
          CASE( 'day' )
            UNITFAC = 365.
          CASE( 'hr' )
            UNITFAC = 365. * 24.
          CASE( 'min' )
            UNITFAC = 365. * 24. * 60.
          CASE( 's' )
            UNITFAC = 365. * 24. * 3600.
          CASE DEFAULT
            CALL CASE_NOT_FOUND( UNITFAC )
        END SELECT

      CASE( 'day' )

        SELECT CASE( BUF2 )
          CASE( 'yr' )
            UNITFAC = 1. / 365.
          CASE( 'day' )
          CASE( 'hr' )
            UNITFAC = 24.
          CASE( 'min' )
            UNITFAC = 24. * 60.
          CASE( 's' )
            UNITFAC = 24. * 3600.
          CASE DEFAULT
            CALL CASE_NOT_FOUND( UNITFAC )
        END SELECT

      CASE( 'hr' )

        SELECT CASE( BUF2 )
          CASE( 'yr' )
            UNITFAC = 1. / ( 365. * 24. )
          CASE( 'day' )
            UNITFAC = 1. / 24.
          CASE( 'hr' )
          CASE( 'min' )
            UNITFAC = 60.
          CASE( 's' )
            UNITFAC = 3600.
          CASE DEFAULT
            CALL CASE_NOT_FOUND( UNITFAC )
        END SELECT

      CASE( 'min' )

        SELECT CASE( BUF2 )
          CASE( 'yr' )
            UNITFAC = 1. / ( 365. * 24. * 60. )
          CASE( 'day' )
            UNITFAC = 1. / ( 24. * 60. )
          CASE( 'hr' )
            UNITFAC = 1. / 60.
          CASE( 'min' )
          CASE( 's' )
            UNITFAC = 60.
          CASE DEFAULT
            CALL CASE_NOT_FOUND( UNITFAC )
        END SELECT

      CASE( 's' )

        SELECT CASE( BUF2 )
          CASE( 'yr' )
            UNITFAC = 1. / ( 365. * 24. * 3600. )
          CASE( 'day' )
            UNITFAC = 1. / ( 24. * 3600. )
          CASE( 'hr' )
            UNITFAC = 1. / 3600.
          CASE( 'min' )
            UNITFAC = 1. / 60.
          CASE( 's' )
          CASE DEFAULT
            CALL CASE_NOT_FOUND( UNITFAC )
        END SELECT

      CASE( 'km' )

        SELECT CASE( BUF2 )
          CASE( 'km' )
          CASE( 'miles' )
            UNITFAC = 0.6213712
          CASE( 'm' )
            UNITFAC = 1000.
          CASE( 'ft' )
            UNITFAC = 3280.84
          CASE DEFAULT
            CALL CASE_NOT_FOUND( UNITFAC )
        END SELECT

      CASE( 'miles' )
        SELECT CASE( BUF2 )
          CASE( 'km' )
            UNITFAC = 1.609344
          CASE( 'miles' )
          CASE( 'm' )
            UNITFAC = 1609.344
          CASE( 'ft' )
            UNITFAC = 5280.
          CASE DEFAULT
            CALL CASE_NOT_FOUND( UNITFAC )
        END SELECT

      CASE( 'm' )

        SELECT CASE( BUF2 )
          CASE( 'km' )
            UNITFAC = 0.001
          CASE( 'miles' )
            UNITFAC = 6.213712E-4
          CASE( 'm' )
          CASE( 'ft' )
            UNITFAC = 3.2808399
          CASE DEFAULT
            CALL CASE_NOT_FOUND( UNITFAC )
        END SELECT

      CASE( 'ft' )

        SELECT CASE( BUF2 )
          CASE( 'km' )
            UNITFAC = 0.0003048
          CASE( 'miles' )
            UNITFAC = 1.893939E-4
          CASE( 'm' )
            UNITFAC = 0.3048
          CASE( 'ft' )
          CASE DEFAULT
            CALL CASE_NOT_FOUND( UNITFAC )
        END SELECT

      CASE DEFAULT
        CALL CASE_NOT_FOUND( UNITFAC )

    END SELECT

    !.......  Inckude any prefix adjustments
    UNITFAC = UNITFAC * PRECNV

    RETURN


CONTAINS    !******************  INTERNAL SUBPROGRAMS  *****************************


    !.......  This subprogram checks for an internal error of the
    !         units conversion case not being programmed in the code
    SUBROUTINE CASE_NOT_FOUND( FACTOR )

        !.......  Subprogram arguments
        REAL, INTENT (OUT) :: FACTOR

        !.......   Local variables
        CHARACTER(300) MESG

        !----------------------------------------------------------------------

        MESG = 'NOTE: Units conversion case not found for: "' //        &
               TRIM( BUF1 ) // '" with "' // TRIM( BUF2 ) // '"'//      &
               CRLF() // BLANK10// 'No conversion will be made.'
        CALL M3MSG2( MESG )
        FACTOR = -1.

        RETURN

    END SUBROUTINE CASE_NOT_FOUND

END FUNCTION UNITFAC


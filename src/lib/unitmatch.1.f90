
SUBROUTINE UNITMATCH( UNITBUF )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine removes or adds trailing "s" to units to make them
!      consistent with the syntax all calling programs use for that unit.
!      If a fraction of units is given, both the numerator and
!      denominator are modified.
!
!  PRECONDITIONS REQUIRED:

!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created by M. Houyoux 11/2000
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
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
! Last updated: $Date$
!
!***************************************************************************
    USE M3UTILIO

    IMPLICIT NONE

!...........   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN OUT) :: UNITBUF

!...........   Local variables
    INTEGER                I, L2

    LOGICAL                FRACFLAG  ! true: units are a fraction

    CHARACTER(NAMLEN3) CNUM    ! tmp numerator
    CHARACTER(NAMLEN3) CDEN    ! tmp denominator (if any)
    CHARACTER(NAMLEN3) PREFIX1 ! tmp numerator prefix (if any)
    CHARACTER(NAMLEN3) PREFIX2 ! tmp denominator prefix (if any)

    CHARACTER(16) :: PROGNAME = 'UNITMATCH' ! program name

!***********************************************************************
!   begin body of subroutine UNITMATCH

!.........  Initialize prefixes
    PREFIX1 = ' '
    PREFIX2 = ' '

!.........  Separate out the numerator and denominator, if any
    I  = INDEX( UNITBUF, '/' )
    L2 = LEN_TRIM( UNITBUF )
    IF( I .GT. 0 ) THEN
        FRACFLAG = .TRUE.
        CNUM = UNITBUF( 1  :I-1 )
        CDEN = UNITBUF( I+1:L2  )
    ELSE
        FRACFLAG = .FALSE.
        CNUM = UNITBUF
        CDEN = ' '
    END IF

!.........  Separate out any leading adjustments (e.g., 10E6) from numerator
!.........  Make sure the numerator is left-justified
    I = INDEX( CNUM, '10E' )

    IF( I .GT. 0 ) THEN
        I = INDEX( CNUM, ' ' )

        IF( I .GT. 0 ) THEN
            PREFIX1 = CNUM( 1:I )
            CNUM    = UNITBUF( I+1:LEN_TRIM( CNUM ) )
            CNUM    = ADJUSTL( CNUM )
        END IF

    ELSE
        CNUM = ADJUSTL( CNUM )

    END IF

!.........  Separate out any leading adjustments (e.g., 10E6) from denominator
!.........  Make sure the numerator is left-justified
    IF( FRACFLAG ) THEN

        I = INDEX( CDEN, '10E' )

        IF( I .GT. 0 ) THEN
            I = INDEX( CDEN, ' ' )

            IF( I .GT. 0 ) THEN
                PREFIX1 = CDEN( 1:I )
                CDEN    = UNITBUF( I+1:LEN_TRIM( CDEN ) )
                CDEN    = ADJUSTL( CDEN )
            END IF

        ELSE
            CDEN = ADJUSTL( CDEN )

        END IF

    END IF

!.........  Reset name of units for numerator
    CALL RENAME_UNITS( CNUM )

!.........  Reset name of units for denominator
    IF( FRACFLAG ) CALL RENAME_UNITS( CDEN )

!.........  Piece together output units
    IF( PREFIX1 .NE. ' ' ) THEN
        L2 = LEN_TRIM( PREFIX1 )
        CNUM = PREFIX1( 1:L2 ) // ' ' // CNUM
    END IF

    IF( PREFIX2 .NE. ' ' ) THEN
        L2 = LEN_TRIM( PREFIX2 )
        CNUM = PREFIX2( 1:L2 ) // ' ' // CDEN
    END IF

    IF( FRACFLAG ) THEN
        L2 = LEN_TRIM( CNUM )
        UNITBUF = CNUM( 1:L2 ) // '/' // CDEN

    ELSE
        UNITBUF = CNUM

    END IF

    RETURN

CONTAINS

!.............  This internal subprogram changes the units to the names
!               to be consistent.
    SUBROUTINE RENAME_UNITS( BUFFER )

        CHARACTER(*), INTENT( IN OUT ) :: BUFFER

!.............  Convert units to SMOKE notation

        SELECT CASE( BUFFER )
          CASE( 'gms' )
            BUFFER = 'g'
          CASE( 'gm' )
            BUFFER = 'g'
          CASE( 'gram' )
            BUFFER = 'g'
          CASE( 'grams' )
            BUFFER = 'g'
          CASE( 'kilogram' )
            BUFFER = 'kg'
          CASE( 'kilograms' )
            BUFFER = 'kg'
          CASE( 'kgs' )
            BUFFER = 'kg'
          CASE( 'ton' )
            BUFFER = 'tons'
          CASE( 'meter' )
            BUFFER = 'm'
          CASE( 'meters' )
            BUFFER = 'm'
          CASE( 'feet' )
            BUFFER = 'ft'
          CASE( 'mile' )
            BUFFER = 'miles'
          CASE( 'mi' )
            BUFFER = 'miles'
          CASE( 'mole' )
            BUFFER = 'moles'
          CASE( 'gm mole' )
            BUFFER = 'moles'
          CASE( 'gm moles' )
            BUFFER = 'moles'
          CASE( 'gm-mole' )
            BUFFER = 'moles'
          CASE( 'gm-moles' )
            BUFFER = 'moles'
          CASE( 'g mole' )
            BUFFER = 'moles'
          CASE( 'g moles' )
            BUFFER = 'moles'
          CASE( 'g-mole' )
            BUFFER = 'moles'
          CASE( 'g-moles' )
            BUFFER = 'moles'
          CASE( 'years' )
            BUFFER = 'yr'
          CASE( 'year' )
            BUFFER = 'yr'
          CASE( 'yrs' )
            BUFFER = 'yr'
          CASE( 'h' )
            BUFFER = 'hr'
          CASE( 'hours' )
            BUFFER = 'hr'
          CASE( 'hour' )
            BUFFER = 'hr'
          CASE( 'hrs' )
            BUFFER = 'hr'
          CASE( 'dy' )
            BUFFER = 'day'
          CASE( 'dys' )
            BUFFER = 'day'
          CASE( 'days' )
            BUFFER = 'day'
          CASE( 'mins' )
            BUFFER = 'min'
          CASE( 'mns' )
            BUFFER = 'min'
          CASE( 'mn' )
            BUFFER = 'min'
          CASE( 'minute' )
            BUFFER = 'min'
          CASE( 'minutes' )
            BUFFER = 'min'
          CASE( 'secs' )
            BUFFER = 's'
          CASE( 'seconds' )
            BUFFER = 's'
          CASE( 'second' )
            BUFFER = 's'
          CASE( 'sec' )
            BUFFER = 's'
        END SELECT

    END SUBROUTINE RENAME_UNITS

END SUBROUTINE UNITMATCH


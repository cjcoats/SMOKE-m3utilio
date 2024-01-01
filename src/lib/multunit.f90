
CHARACTER(*) FUNCTION MULTUNIT( UNIT1, UNIT2 )

    !***********************************************************************
    !  function body starts at line
    !
    !  DESCRIPTION:
    !      This function combines the two units by comparing the numerators
    !      and denominators and cancelling any of the same units. It works
    !      only for one string in the numerator and in the denominator, but
    !      could be enhanced to be more sophisticated.
    !
    !  PRECONDITIONS REQUIRED:

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

    !......   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: UNIT1        ! first unit
    CHARACTER(*), INTENT (IN) :: UNIT2        ! second unit

    !......   Other local variables
    INTEGER        K1, K2, L1, L2, LD1, LD2, LN1, LN2

    CHARACTER(NAMLEN3) :: DEN1 = ' '
    CHARACTER(NAMLEN3) :: DEN2 = ' '
    CHARACTER(NAMLEN3) :: NUM1 = ' '
    CHARACTER(NAMLEN3) :: NUM2 = ' '

    !***********************************************************************
    !   begin body of function MULTUNIT

    !......  Retrieve length of unit strings
    L1 = LEN_TRIM( UNIT1 )
    L2 = LEN_TRIM( UNIT2 )

    !......  Cases where one or more of the provided units are blank
    IF( L1 .LE. 0 .AND. L2 .LE. 0 ) THEN

        MULTUNIT = 'unitless'
        RETURN

    ELSE IF( L1 .LE. 0 ) THEN

        MULTUNIT = UNIT2
        RETURN

    ELSE IF( L2 .LE. 0 ) THEN

        MULTUNIT = UNIT1
        RETURN

    END IF

    !......  Determine positions of the divide-by symbol in the units
    K1 = INDEX( UNIT1, '/' )
    K2 = INDEX( UNIT2, '/' )

    IF( K1 .LE. 0 ) K1 = L1 + 1
    IF( K2 .LE. 0 ) K2 = L2 + 1

    !......  Define numerators, denominators, and their lengths
    NUM1 = UNIT1( 1 : K1-1 )
    IF( K1 .LT. L1 ) DEN1 = UNIT1( K1+1 : L1 )

    LN1  = LEN_TRIM( NUM1 )
    LD1  = LEN_TRIM( DEN1 )

    NUM2 = UNIT2( 1 : K2-1 )
    IF( K2 .LT. L2 ) DEN2 = UNIT2( K2+1 : L2 )

    LN2  = LEN_TRIM( NUM2 )
    LD2  = LEN_TRIM( DEN2 )

    !......  Set units by comparing numerators and denominators
    IF( ( NUM1 .EQ. DEN2 .AND. NUM2 .EQ. DEN1 ) .OR.&
        ( NUM1 .EQ. DEN1 .AND. NUM2 .EQ. DEN2 )      ) THEN

        MULTUNIT = 'unitless'

    ELSE IF( NUM1 .EQ. DEN2 ) THEN

        MULTUNIT = NUM2( 1:LN2 )// '/'// DEN1

    ELSE IF( NUM2 .EQ. DEN1 ) THEN

        MULTUNIT = NUM1( 1:LN1 )// '/'// DEN2

    ELSE IF( NUM1 .EQ. DEN1 ) THEN

        MULTUNIT = UNIT2

    ELSE IF( NUM2 .EQ. DEN2 ) THEN

        MULTUNIT = UNIT1

    ELSE IF( NUM1 .EQ. '1' .AND. NUM2 .EQ. '1' ) THEN

        MULTUNIT = '1/'// DEN1( 1:LD1 )// '*'// DEN2

    ELSE IF( NUM1 .EQ. '1' ) THEN

        MULTUNIT = NUM2( 1:LN2 )// '/'// DEN1( 1:LD1 )// '*'// DEN2

    ELSE IF( NUM2 .EQ. '1' ) THEN

        MULTUNIT = NUM1( 1:LN1 )// '/'// DEN1( 1:LD1 )// '*'// DEN2

    ELSE IF( DEN1 .EQ. ' ' .AND. DEN2 .EQ. ' ' ) THEN

        MULTUNIT = NUM1( 1:LN1 )// '*'// NUM2( 1:LN2 )

    ELSE IF( DEN1 .EQ. ' ' ) THEN

        MULTUNIT = NUM1( 1:LN1 )// '*'// NUM2( 1:LN2 )// '/'// DEN2

    ELSE IF( DEN2 .EQ. ' ' ) THEN

        MULTUNIT = NUM1( 1:LN1 )// '*'// NUM2( 1:LN2 )// '/'// DEN1

    ELSE

        MULTUNIT = NUM1( 1:LN1 )// '*'// NUM2( 1:LN2 )// '/'// DEN1( 1:LD1 )// '*'// DEN2

    END IF

    RETURN

END FUNCTION MULTUNIT


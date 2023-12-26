
SUBROUTINE BLDCSRC( CFIP , STR2, STR3, STR4,&
&                    STR5, STR6, STR7, CCOD , CSRC )

    !***********************************************************************
    !  subroutine body starts at line 97
    !
    !  DESCRIPTION:
    !      This subroutine combines the source characteristics into a single
    !      string.  It right justifies the fields and maintains the correct
    !      segment lengths.  Note that if the results from this function are
    !      used in the calling program with concatenation, the programmer should
    !      first use LEN_TRIM on the function result before doing the concatenation
    !
    !  PRECONDITIONS REQUIRED:
    !      Values for at least CFIP, PLANT, CHAR1, and CCOD, with other CHAR* values
    !      optional.  The width of the input strings *must* be consistent with
    !      the right-justified width desired for output, which is defined by
    !      the start and end positions of ARBEGL3/ARENDL3 for area sources,
    !      MBBEGL3/MBENDL3 for mobile sources, and PTBEGL3/PTENDL3 for point
    !      sources.  These arrays are defined in EMSTRG3.h90
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !     Subroutines: Models-3 subroutines
    !     Functions: Models-3 functions
    !
    !  REVISION  HISTORY:
    !       Created by M. Houyoux 10/98
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
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

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT(IN OUT) :: CFIP      ! string FIPS code
    CHARACTER(*), INTENT(IN OUT) :: STR2      ! ar:SCC, mb:SCC, pt: plant
    CHARACTER(*), INTENT(IN OUT) :: STR3      ! ar:blk, mb:VTYPE, pt: plt char 1
    CHARACTER(*), INTENT(IN OUT) :: STR4      ! ar:blk, mb:LINK, pt: plt char 2
    CHARACTER(*), INTENT(IN OUT) :: STR5      ! ar:blk, mb:blk, pt: plt char 3
    CHARACTER(*), INTENT(IN OUT) :: STR6      ! ar:blk, mb:blk, pt: plt char 4
    CHARACTER(*), INTENT(IN OUT) :: STR7      ! ar:blk, mb:blk, pt: plt char 5
    CHARACTER(*), INTENT(IN OUT ):: CCOD      ! string of int postn of pollutant
    CHARACTER(*), INTENT(OUT)    :: CSRC      ! concatenated result

    !.......   Other local variables
    INTEGER         L
    INTEGER         OUTLEN

    CHARACTER(300)  MESG
    CHARACTER(300)  STRING

    !***********************************************************************
    !   begin body of subroutine BLDCSRC

    !.....  Determine allocated length of string used for output
    OUTLEN = LEN( CSRC )

    !.....  Interpret -9 as blank
    IF( CFIP .EQ. '-9' ) CFIP = ' '
    IF( STR2 .EQ. '-9' ) STR2 = ' '
    IF( STR3 .EQ. '-9' ) STR3 = ' '
    IF( STR4 .EQ. '-9' ) STR4 = ' '
    IF( STR5 .EQ. '-9' ) STR5 = ' '
    IF( STR6 .EQ. '-9' ) STR6 = ' '
    IF( STR7 .EQ. '-9' ) STR7 = ' '
    IF( CCOD .EQ. '-9' ) CCOD = ' '

    !.....  Store right-justified entries in the output string

    STRING = ' '
    STRING = ADJUSTR( CFIP ) // &
             ADJUSTR( STR2 ) // &
             ADJUSTR( STR3 ) // &
             ADJUSTR( STR4 ) // &
             ADJUSTR( STR5 ) // &
             ADJUSTR( STR6 ) // &
             ADJUSTR( STR7 )
    STRING( POLPOS3:ALLLEN3 ) = ADJUSTR( CCOD )

    L = LEN_TRIM( STRING )

    IF( L .GT. OUTLEN ) THEN
        MESG= 'INTERNAL ERROR: Function BLDCSRC has been modified improperly'
        CALL M3MSG2( MESG )
        CALL M3EXIT( 'BLDCSRC', 0, 0, ' ', 2 )

    ELSE
        CSRC = STRING( 1:OUTLEN )

    ENDIF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

END SUBROUTINE BLDCSRC


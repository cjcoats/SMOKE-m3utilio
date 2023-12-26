
SUBROUTINE MRGUNITS

    !***********************************************************************
    !  subroutine MRGUNITS body starts at line 93
    !
    !  DESCRIPTION:
    !      The purpose of this subroutine
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 2/99 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***********************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !        System
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
    !***********************************************************************
    USE M3UTILIO

    !.....  MODULES for public variables
    !.....  This module contains the major data structure and control flags
    USE MODMERGE, ONLY: BIOGFAC, BIOTFAC, GRDFAC, TOTFAC,       &
                        BIOUNIT, GRDUNIT, TOTUNIT,              &
                        SPCUNIT, TOUNITS, EMIDX,                &
                        SFLAG, TFLAG, BFLAG, NIPPA, NMSPC, NUNITS

    IMPLICIT NONE

    !.....   INCLUDES:

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.....   EXTERNAL FUNCTIONS and their descriptions:

    CHARACTER(16), EXTERNAL :: MULTUNIT
    REAL         , EXTERNAL :: UNITFAC

    !.....   Other local variables

    INTEGER         IOS      ! tmp I/O status
    INTEGER         J, L, V        ! counter

    REAL            FAC1, FAC2

    CHARACTER(300)  BUFFER       ! text buffer
    CHARACTER(300)  MESG     ! message buffer

    CHARACTER(NAMLEN3) GRDUNIT_I       ! initialized gridded outputs units
    CHARACTER(NAMLEN3) GDEN_I      ! initialized gridded denominator
    CHARACTER(NAMLEN3) GNUM_I      ! initialized gridded numerator
    CHARACTER(NAMLEN3) GDEN        ! work gridded denominator
    CHARACTER(NAMLEN3) GNUM        ! work gridded numerator
    CHARACTER(NAMLEN3) GRDBUF      ! work gridded output units
    CHARACTER(NAMLEN3) GRDENV      ! gridded output units from envrmt

    CHARACTER(NAMLEN3) BIOUNIT_T       ! initialized biog totals units
    CHARACTER(NAMLEN3) TOTUNIT_I       ! initialized output totals units
    CHARACTER(NAMLEN3) TDEN_I      ! initialized totals denominator
    CHARACTER(NAMLEN3) TNUM_I      ! initialized totals numerator
    CHARACTER(NAMLEN3) TDEN        ! work totals denominator
    CHARACTER(NAMLEN3) TNUM        ! work totals numerator
    CHARACTER(NAMLEN3) TOTBUF      ! work output totals  units
    CHARACTER(NAMLEN3) TOTENV      ! output totals units from envrmt

    CHARACTER(16), PARAMETER :: PROGNAME = 'MRGUNITS'     ! program name

    !***********************************************************************
    !   begin body of subroutine MRGUNITS

    !.....  Allocate memory for units conversion factors and units.

    !.....  If using speciation, allocate arrays to number of species;
    !       otherwise, use total number of activities and pollutants.
    !       This approach assumes that each species has the same units
    !       regardless of which pollutant it comes from, i.e.
    !       ALD2 from VOC must have the same units as ALD2 from ACETALD
    IF( SFLAG ) THEN
        NUNITS = NMSPC
    ELSE
        NUNITS = NIPPA
    END IF

    ALLOCATE( GRDFAC( NUNITS ),    &
              TOTFAC( NUNITS ),    &
             GRDUNIT( NUNITS ),    &
             TOTUNIT( NUNITS ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TOTUNIT', PROGNAME )

    !.....  Initialize all
    BIOGFAC = 1.
    BIOTFAC = 1.
    GRDFAC  = 1.       ! array
    TOTFAC  = 1.       ! array
    GRDUNIT = ' '      ! array
    TOTUNIT = ' '      ! array

    !.....  Retrieve variables for setting the output units for gridded and
    !       country/state/county total data
    BUFFER = 'Units for output gridded emissions'
    CALL ENVSTR( 'MRG_GRDOUT_UNIT', BUFFER, ' ', GRDENV, IOS)
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_GRDOUT_UNIT"', 2 )
    END IF

    BUFFER = 'Units for output state/county total emissions'
    CALL ENVSTR( 'MRG_TOTOUT_UNIT', BUFFER, ' ', TOTENV, IOS)
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_TOTOUT_UNIT"', 2 )
    END IF

    !.....  Loop through pollutants, and create output units accordingly
    !.....  For speciation, use the first unit for the speciation units from
    !       a given pollutants speciation factors
    DO V = 1, NUNITS

        !.....  Initialize the output units for speciation
        IF( SFLAG ) THEN

            CALL UNITMATCH( SPCUNIT( V ) )

            !.....  Get pollutant position for this species
            !       NOTE: If multiple pollutants contribute to this species,
            !             we will only be using the last one in master list;
            !             this shouldn't be a problem unless the original
            !             pollutants have different units
            J = EMIDX( V )
            CALL UNITMATCH( TOUNITS( J ) )

            GRDUNIT_I = MULTUNIT( SPCUNIT( V ), TOUNITS( J ) )
            IF( TFLAG ) THEN
                TOTUNIT_I = MULTUNIT( GRDUNIT_I, 'hr/day' )
            ELSE
                TOTUNIT_I = GRDUNIT_I
            END IF

        !.....  For non-speciation, inventory pollutants or activities could
        !       have a variety of units depending on temporalized emissions or
        !       inventory emissions and activities.  Make sure that if the temporal
        !       resolution is per hour, that the TOTUNIT is still set as per day.
        ELSE
            J = V
            CALL UNITMATCH( TOUNITS( J ) )
            GRDUNIT_I = TOUNITS( J )

            L = INDEX( TOUNITS( J ), 'hr' )
            IF( L .GT. 0 ) THEN
                TOTUNIT_I = MULTUNIT( TOUNITS( J ), 'hr/day' )
            ELSE
                TOTUNIT_I = TOUNITS( J )
            END IF

        END IF

        !.....  Set the trial units
        GRDBUF = GRDUNIT_I
        TOTBUF = TOTUNIT_I
        IF( GRDENV .NE. ' ' ) GRDBUF = GRDENV
        IF( TOTENV .NE. ' ' ) TOTBUF = TOTENV

        !.....  Set the numerators and denominators
        L = INDEX( GRDUNIT_I, '/' )
        GNUM_I = ADJUSTL( GRDUNIT_I(   1:L-1     ) )
        GDEN_I = ADJUSTL( GRDUNIT_I( L+1:NAMLEN3 ) )

        L = INDEX( TOTUNIT_I, '/' )
        TNUM_I = ADJUSTL( TOTUNIT_I(   1:L-1     ) )
        TDEN_I = ADJUSTL( TOTUNIT_I( L+1:NAMLEN3 ) )

        L = INDEX( GRDBUF, '/' )
        GNUM = ADJUSTL( GRDBUF(   1:L-1     ) )
        GDEN = ADJUSTL( GRDBUF( L+1:NAMLEN3 ) )

        L = INDEX( TOTBUF, '/' )
        TNUM = ADJUSTL( TOTBUF(   1:L-1     ) )
        TDEN = ADJUSTL( TOTBUF( L+1:NAMLEN3 ) )

        !.....  Get factor for the numerators for the gridded outputs...
        IF( SFLAG ) THEN
            FAC1 = UNITFAC( SPCUNIT( V ), GRDBUF, .TRUE. )      ! speciation
        ELSE
            FAC1 = UNITFAC( TOUNITS( J ), GRDBUF, .TRUE. )
        END IF

        !.....  Get factor for the denominators for the gridded outputs
        FAC2 = UNITFAC( TOUNITS( J ), GRDBUF, .FALSE. )

        !.....  In case E.V. setting was bogus rebuild output units based
        !       on which factors were valid
        !.....  Also set negative factors (from unknown conversions) to 1.
        CALL CORRECT_UNITS( GNUM_I, GDEN_I, GNUM, GDEN, GRDBUF )

        !.....  Set factors for gridded outputs
        GRDFAC( V ) = FAC1 / FAC2

        !.....  Get conversion factor for the numerators for totals
        IF( SFLAG ) THEN
            FAC1 = UNITFAC( SPCUNIT( V ), TOTBUF, .TRUE. )          ! speciation
        ELSE
            FAC1 = UNITFAC( TOUNITS( J ), TOTBUF, .TRUE. )
        END IF

        !.....  Get factors for the denominators for the totals.  Note that
        !       the hourly data are output as daily totals.
        FAC2 = UNITFAC( TOTUNIT_I, TOTBUF, .FALSE. )

        !.....  In case E.V. setting was bogus rebuild output units based
        !       on which factors were valid
        !.....  Also set negative factors (from unknown conversions) to 1.
        CALL CORRECT_UNITS( TNUM_I, TDEN_I, TNUM, TDEN, TOTBUF )

        !.....  Set factors for totaled outputs
        TOTFAC( V ) = FAC1 / FAC2

        !.....  Set the output units per pollutant/activity
        GRDUNIT( V ) = GRDBUF
        TOTUNIT( V ) = TOTBUF

    END DO

    !.....  If biogenics, then get the factor needed for converting gridded
    !       outputs and totals outputs.
    IF( BFLAG ) THEN

        !.....  Set the trial units. NOTE - this could be too simplistic.
        GRDBUF = GRDUNIT( 1 )
        TOTBUF = TOTUNIT( 1 )
        IF( GRDENV .NE. ' ' ) GRDBUF = GRDENV
        IF( TOTENV .NE. ' ' ) TOTBUF = TOTENV

        !.....  Get factor for the numerators for the gridded outputs...
        CALL UNITMATCH( BIOUNIT )
        FAC1 = UNITFAC( BIOUNIT, GRDBUF, .TRUE. )

        !.....  Get factor for the denominators for the gridded outputs...
        FAC2 = UNITFAC( BIOUNIT, GRDBUF, .FALSE. )

        BIOGFAC = FAC1 / FAC2

        BIOUNIT_T = MULTUNIT( BIOUNIT, 'hr/day' )

        !.....  Get factor for the numerators for the output totals...
        FAC1 = UNITFAC( BIOUNIT_T, TOTBUF, .TRUE. )          ! speciation

        !.....  Get factors for the denominators for the totals.  Note that
        !       the hourly data are output as daily totals.
        FAC2 = UNITFAC( BIOUNIT_T, TOTBUF, .FALSE. )

        IF ( FAC1 .LT. 0. ) FAC1 = 1.
        IF ( FAC2 .LT. 0. ) FAC2 = 1.

        !.....  Set factors for gridded outputs
        BIOTFAC = FAC1 / FAC2

    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.....   Internal buffering formats.....94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

    !******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

    !.....  This internal subprogram corrects badly formatted units
    SUBROUTINE CORRECT_UNITS( NUM_I, DEN_I, NUM, DEN, OUTUNIT )

        !.....  Subprogram arguments
        CHARACTER(NAMLEN3), INTENT(IN   ) :: NUM_I
        CHARACTER(NAMLEN3), INTENT(IN   ) :: DEN_I
        CHARACTER(NAMLEN3), INTENT(INOUT) :: NUM
        CHARACTER(NAMLEN3), INTENT(INOUT) :: DEN
        CHARACTER(NAMLEN3), INTENT(  OUT) :: OUTUNIT

        !.....  Local variables
        INTEGER L1,L2

        !----------------------------------------------------------------------

        IF( FAC1 .LT. 0. ) THEN
            NUM = NUM_I
            FAC1 = 1.
        END IF
        IF( FAC2 .LT. 0. ) THEN
            DEN = DEN_I
            FAC2 = 1.
        END IF

        OUTUNIT =  TRIM( NUM ) // '/' // DEN
        
        RETURN

    END SUBROUTINE CORRECT_UNITS

END SUBROUTINE MRGUNITS

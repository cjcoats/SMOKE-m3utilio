
        SUBROUTINE MRGUNITS

C***********************************************************************
C  subroutine MRGUNITS body starts at line 93
C
C  DESCRIPTION:
C      The purpose of this subroutine
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 2/99 by M. Houyoux
C       Version 11/2023 by CJC:  USE M3UTILIO and related changes
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$
C
C***********************************************************************
        USE M3UTILIO

C.........  MODULES for public variables
C.........  This module contains the major data structure and control flags
        USE MODMERGE, ONLY: BIOGFAC, BIOTFAC, GRDFAC, TOTFAC,
     &                      BIOUNIT, GRDUNIT, TOTUNIT,
     &                      SPCUNIT, TOUNITS, EMIDX,
     &                      SFLAG, TFLAG, BFLAG, NIPPA, NMSPC, NUNITS

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   EXTERNAL FUNCTIONS and their descriptions:

        CHARACTER(16), EXTERNAL :: MULTUNIT
        REAL         , EXTERNAL :: UNITFAC

C...........   Other local variables

        INTEGER         IOS      ! tmp I/O status
        INTEGER         J, L, V        ! counter

        REAL            FAC1, FAC2

        CHARACTER(300)  BUFFER   ! text buffer
        CHARACTER(300)  MESG     ! message buffer

        CHARACTER(NAMLEN3) GRDUNIT_I   ! initialized gridded outputs units
        CHARACTER(NAMLEN3) GDEN_I      ! initialized gridded denominator
        CHARACTER(NAMLEN3) GNUM_I      ! initialized gridded numerator
        CHARACTER(NAMLEN3) GDEN        ! work gridded denominator
        CHARACTER(NAMLEN3) GNUM        ! work gridded numerator
        CHARACTER(NAMLEN3) GRDBUF      ! work gridded output units
        CHARACTER(NAMLEN3) GRDENV      ! gridded output units from envrmt

        CHARACTER(NAMLEN3) BIOUNIT_T   ! initialized biog totals units
        CHARACTER(NAMLEN3) TOTUNIT_I   ! initialized output totals units
        CHARACTER(NAMLEN3) TDEN_I      ! initialized totals denominator
        CHARACTER(NAMLEN3) TNUM_I      ! initialized totals numerator
        CHARACTER(NAMLEN3) TDEN        ! work totals denominator
        CHARACTER(NAMLEN3) TNUM        ! work totals numerator
        CHARACTER(NAMLEN3) TOTBUF      ! work output totals  units
        CHARACTER(NAMLEN3) TOTENV      ! output totals units from envrmt

        CHARACTER(16) :: PROGNAME = 'MRGUNITS' ! program name

C***********************************************************************
C   begin body of subroutine MRGUNITS

C.........  Allocate memory for units conversion factors and units.

C.........  If using speciation, allocate arrays to number of species;
C           otherwise, use total number of activities and pollutants.
C           This approach assumes that each species has the same units
C           regardless of which pollutant it comes from, i.e.
C           ALD2 from VOC must have the same units as ALD2 from ACETALD
        IF( SFLAG ) THEN
            NUNITS = NMSPC
        ELSE
            NUNITS = NIPPA
        END IF

        ALLOCATE( GRDFAC( NUNITS ),
     &            TOTFAC( NUNITS ),
     &           GRDUNIT( NUNITS ),
     &           TOTUNIT( NUNITS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TOTUNIT', PROGNAME )

C.........  Initialize all
        BIOGFAC = 1.
        BIOTFAC = 1.
        GRDFAC  = 1.   ! array
        TOTFAC  = 1.   ! array
        GRDUNIT = ' '  ! array
        TOTUNIT = ' '  ! array

C.........  Retrieve variables for setting the output units for gridded and
C           country/state/county total data
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

C.........  Loop through pollutants, and create output units accordingly
C.........  For speciation, use the first unit for the speciation units from
C           a given pollutants speciation factors
        DO V = 1, NUNITS

C.............  Initialize the output units for speciation
            IF( SFLAG ) THEN

                CALL UNITMATCH( SPCUNIT( V ) )

C.................  Get pollutant position for this species
C                   NOTE: If multiple pollutants contribute to this species,
C                         we will only be using the last one in master list;
C                         this shouldn't be a problem unless the original
C                         pollutants have different units
                J = EMIDX( V )
                CALL UNITMATCH( TOUNITS( J ) )

                GRDUNIT_I = MULTUNIT( SPCUNIT( V ), TOUNITS( J ) )
                IF( TFLAG ) THEN
                    TOTUNIT_I = MULTUNIT( GRDUNIT_I, 'hr/day' )
                ELSE
                    TOTUNIT_I = GRDUNIT_I
                END IF

C.............  For non-speciation, inventory pollutants or activities could
C           have a variety of units depending on temporalized emissions or
C           inventory emissions and activities.  Make sure that if the temporal
C           resolution is per hour, that the TOTUNIT is still set as per day.
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

C.............  Set the trial units
            GRDBUF = GRDUNIT_I
            TOTBUF = TOTUNIT_I
            IF( GRDENV .NE. ' ' ) GRDBUF = GRDENV
            IF( TOTENV .NE. ' ' ) TOTBUF = TOTENV

C.............  Set the numerators and denominators
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

C.............  Get factor for the numerators for the gridded outputs...
            IF( SFLAG ) THEN
                FAC1 = UNITFAC( SPCUNIT( V ), GRDBUF, .TRUE. )  ! speciation
            ELSE
                FAC1 = UNITFAC( TOUNITS( J ), GRDBUF, .TRUE. )
            END IF

C.............  Get factor for the denominators for the gridded outputs
            FAC2 = UNITFAC( TOUNITS( J ), GRDBUF, .FALSE. )

C.............  In case E.V. setting was bogus rebuild output units based
C               on which factors were valid
C.............  Also set negative factors (from unknown conversions) to 1.
            CALL CORRECT_UNITS( GNUM_I, GDEN_I, GNUM, GDEN, GRDBUF )

C.............  Set factors for gridded outputs
            GRDFAC( V ) = FAC1 / FAC2

C.............  Get conversion factor for the numerators for totals
            IF( SFLAG ) THEN
                FAC1 = UNITFAC( SPCUNIT( V ), TOTBUF, .TRUE. )  ! speciation
            ELSE
                FAC1 = UNITFAC( TOUNITS( J ), TOTBUF, .TRUE. )
            END IF

C.............  Get factors for the denominators for the totals.  Note that
C               the hourly data are output as daily totals.
            FAC2 = UNITFAC( TOTUNIT_I, TOTBUF, .FALSE. )

C.............  In case E.V. setting was bogus rebuild output units based
C               on which factors were valid
C.............  Also set negative factors (from unknown conversions) to 1.
            CALL CORRECT_UNITS( TNUM_I, TDEN_I, TNUM, TDEN, TOTBUF )

C.............  Set factors for totaled outputs
            TOTFAC( V ) = FAC1 / FAC2

C.............  Set the output units per pollutant/activity
            GRDUNIT( V ) = GRDBUF
            TOTUNIT( V ) = TOTBUF

        END DO

C.........  If biogenics, then get the factor needed for converting gridded
C           outputs and totals outputs.
        IF( BFLAG ) THEN

C.............  Set the trial units. NOTE - this could be too simplistic.
            GRDBUF = GRDUNIT( 1 )
            TOTBUF = TOTUNIT( 1 )
            IF( GRDENV .NE. ' ' ) GRDBUF = GRDENV
            IF( TOTENV .NE. ' ' ) TOTBUF = TOTENV

C.............  Get factor for the numerators for the gridded outputs...
            CALL UNITMATCH( BIOUNIT )
            FAC1 = UNITFAC( BIOUNIT, GRDBUF, .TRUE. )

C.............  Get factor for the denominators for the gridded outputs...
            FAC2 = UNITFAC( BIOUNIT, GRDBUF, .FALSE. )

            BIOGFAC = FAC1 / FAC2

            BIOUNIT_T = MULTUNIT( BIOUNIT, 'hr/day' )

C.............  Get factor for the numerators for the output totals...
            FAC1 = UNITFAC( BIOUNIT_T, TOTBUF, .TRUE. )  ! speciation

C.............  Get factors for the denominators for the totals.  Note that
C               the hourly data are output as daily totals.
            FAC2 = UNITFAC( BIOUNIT_T, TOTBUF, .FALSE. )

            IF ( FAC1 .LT. 0. ) FAC1 = 1.
            IF ( FAC2 .LT. 0. ) FAC2 = 1.

C.............  Set factors for gridded outputs
            BIOTFAC = FAC1 / FAC2

        END IF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats.............94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

C******************  INTERNAL SUBPROGRAMS  *****************************

        CONTAINS

C.............  This internal subprogram corrects badly formatted units
            SUBROUTINE CORRECT_UNITS( NUM_I, DEN_I, NUM, DEN, OUTUNIT )

C.............  Subprogram arguments
            CHARACTER(NAMLEN3) NUM_I
            CHARACTER(NAMLEN3) DEN_I
            CHARACTER(NAMLEN3) NUM
            CHARACTER(NAMLEN3) DEN
            CHARACTER(NAMLEN3) OUTUNIT

C.............  Local variables
            INTEGER L1,L2

C----------------------------------------------------------------------

            IF( FAC1 .LT. 0. ) THEN
                NUM = NUM_I
                FAC1 = 1.
            END IF
            IF( FAC2 .LT. 0. ) THEN
                DEN = DEN_I
                FAC2 = 1.
            END IF

            OUTUNIT =  TRIM( NUM ) // '/' // DEN

            END SUBROUTINE CORRECT_UNITS

        END SUBROUTINE MRGUNITS

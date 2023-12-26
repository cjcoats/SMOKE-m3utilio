
SUBROUTINE REPUNITS( RCNT )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      The REPUNITS routine is reponsible for generating column header
!      units and data conversion factors. Does not create conversion factors
!      for cell area normalization (included in gridding factor array) or
!      population normalization (BINPOPDIV).
!
!  PRECONDITIONS REQUIRED:
!      From previous subroutines, we should have indices defined for
!      which columns are for output for each species and pollutant and
!      the various combinations.
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 8/2000 by M Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!***********************************************************************
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
!***********************************************************************

!...........   MODULES for public variables
!.........  This module contains Smkreport-specific settings
    USE MODREPRT, ONLY: RPT_, AFLAG, OUTUNIT, UCNVFAC, MXINDAT,&
    &                    ALLRPT, ASCDATA, ALLUSET

!.........  This module contains report arrays for each output bin
    USE MODREPBN, ONLY: NSVARS, TODOUT, TOSOUT, SPCTOINV,&
    &                    SPCTOTPR, SLUNIT, SSUNIT

!.........  This module contains the temporal profile tables
    USE MODTMPRL, ONLY: NTPDAT, TPUNIT

!.........  This module contains the information about the source category
    USE MODINFO, ONLY: NIPPA, EAUNIT

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   EXTERNAL FUNCTIONS
    CHARACTER(16), EXTERNAL :: MULTUNIT

!...........   SUBROUTINE ARGUMENTS
    INTEGER, INTENT (IN) :: RCNT    ! report number

!...........   Local allocatable variables
    CHARACTER(IOULEN3) :: LOCUNIT( NIPPA + NTPDAT )

!...........   Other local variables
    INTEGER         E, J, V   ! counters and indices

    INTEGER         IOS               ! i/o status
    INTEGER         NDATA             ! number of data columns
    INTEGER         NV                ! tmp no. variables

    LOGICAL      :: FIRSTIME = .TRUE.  ! true: first time routine called
    LOGICAL      :: SFLAG    = .FALSE. ! true: speciation applies to rpt

    CHARACTER(IOULEN3) TMPUNIT     !  tmp units buffer
    CHARACTER(300)     MESG        !  message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'REPUNITS' ! program name

!***********************************************************************
!   begin body of subroutine REPUNITS

!.........  Report-specific local settings
    IF( AFLAG ) ALLRPT( RCNT )%NUMDATA = ASCDATA
    NDATA = ALLRPT( RCNT )%NUMDATA
    RPT_  = ALLRPT( RCNT )

    SFLAG = ( ALLRPT( RCNT )%USESLMAT .OR.&
    &          ALLRPT( RCNT )%USESSMAT      )

!.........  Set starting units for current report
    NV = NIPPA + NTPDAT
    LOCUNIT = ' '   ! array

    LOCUNIT( 1:NIPPA ) = EAUNIT( 1:NIPPA )   ! array

    IF( RPT_%USEHOUR ) THEN
        LOCUNIT( NIPPA+1:NV ) = TPUNIT( 1:NTPDAT )  ! array
    END IF

!.........  Set up input units and unit conversion factors...

!.........  Allocate memory for needed arrays
    IF( ALLOCATED( OUTUNIT ) ) DEALLOCATE( OUTUNIT, UCNVFAC )
    ALLOCATE( OUTUNIT( NDATA ),&
    &          UCNVFAC( NDATA ), STAT=IOS )
    CALL CHECKMEM( IOS, 'UCNVFAC', PROGNAME )

    OUTUNIT = ' '   ! array
    UCNVFAC = 1.    ! array

!.........  If current report has speciation, loop through species and update
!           units arrays
    DO V = 1, NSVARS

!.............  Set index to data arrays based on speciation status
        E = SPCTOINV( V )
        IF ( RPT_%USEHOUR ) E = SPCTOTPR( V )

        IF ( E .LE. 0 ) CYCLE

!.............  If current variable is a speciated variable and is output for
!               this report
        IF( TOSOUT( V,RCNT )%AGG .GT. 0 ) THEN

!.................  If using mole-speciation matrix
            IF( RPT_%USESLMAT ) THEN
                CALL UNITMATCH( SLUNIT( V ) )
                TMPUNIT = MULTUNIT( SLUNIT( V ), LOCUNIT( E ) )

!.................  If using mass-speciation matrix
            ELSE IF( RPT_%USESSMAT ) THEN
                CALL UNITMATCH( SSUNIT( V ) )
                TMPUNIT = MULTUNIT( SSUNIT( V ), LOCUNIT( E ) )

            END IF

!.................  Set units and conversion factors for appropriate columns
            CALL UPDATE_OUTUNIT( NDATA, TOSOUT( V,RCNT )%SPC,&
            &                     TMPUNIT, OUTUNIT, UCNVFAC )
            CALL UPDATE_OUTUNIT( NDATA, TOSOUT( V,RCNT )%ETPSPC,&
            &                     TMPUNIT, OUTUNIT, UCNVFAC )
            CALL UPDATE_OUTUNIT( NDATA, TOSOUT( V,RCNT )%PRCSPC,&
            &                     TMPUNIT, OUTUNIT, UCNVFAC )
            CALL UPDATE_OUTUNIT( NDATA, TOSOUT( V,RCNT )%SUMETP,&
            &                     TMPUNIT, OUTUNIT, UCNVFAC )
            CALL UPDATE_OUTUNIT( NDATA, TOSOUT( V,RCNT )%SUMPOL,&
            &                     TMPUNIT, OUTUNIT, UCNVFAC )
            CALL UPDATE_OUTUNIT( NDATA, TOSOUT( V,RCNT )%SUMSPC,&
            &                     TMPUNIT, OUTUNIT, UCNVFAC )

        END IF

    END DO

!.........  Now loop through pol/act/e-type and update units arrays
    DO E = 1, NV

!.............  If current variable is a pol/act/e-type and is used for this
!               report
        IF( TODOUT( E,RCNT )%AGG .GT. 0 ) THEN

            TMPUNIT = LOCUNIT( E )

!.................  Set units and conversion factors for appropriate columns
            CALL UPDATE_OUTUNIT( NDATA, TODOUT( E,RCNT )%ETP,&
            &                     TMPUNIT, OUTUNIT, UCNVFAC )
            CALL UPDATE_OUTUNIT( NDATA, TODOUT( E,RCNT )%DAT,&
            &                     TMPUNIT, OUTUNIT, UCNVFAC )

        END IF

    END DO      ! Done loop for setting input units\

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I10, :, 1X ) )

!******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

!.............  This internal function updates the units labels

    SUBROUTINE UPDATE_OUTUNIT( NREC, OUTCOL, I_UNIT,&
    &                           UNITS, FACTOR )

!.............  External functions
        REAL     UNITFAC
        EXTERNAL UNITFAC

!.............  Subprogram arguments
        INTEGER     , INTENT (IN) :: NREC    ! number of available colums
        INTEGER     , INTENT (IN) :: OUTCOL  ! specific column to update
        CHARACTER(*), INTENT (IN) :: I_UNIT ! input unit string
        CHARACTER(*), INTENT(OUT) :: UNITS( NREC )  ! string of units
        REAL        , INTENT(OUT) :: FACTOR( NREC ) ! conversion factors

!.............  Local subprogram variables
        INTEGER      L, LM1, LP1, L1, L2      ! indices

        INTEGER      UIDX    ! tmp unit array index

        REAL         FAC1    ! factor for converting units numerator
        REAL         FAC2    ! factor for converting units denominator

        CHARACTER(IOULEN3) DEN_I     ! tmp input denominator
        CHARACTER(IOULEN3) DEN_O     ! tmp output denominator
        CHARACTER(IOULEN3) O_UNIT    ! tmp output unit
        CHARACTER(IOULEN3) NUM_I     ! tmp input numerator
        CHARACTER(IOULEN3) NUM_O     ! tmp output numerator
        CHARACTER(IOULEN3) T_UNIT    ! tmp unit

!----------------------------------------------------------------------

!.............  Return immediately if current record is not an output column
        IF( OUTCOL .LE. 0 ) RETURN

!.............  Check if output data records is less than or equal to maximum
!               input records for indexing purposes
!.............  If output number is greater, then no SELECT DATA statement
!               was used, so units, if specified, should be retrieved from
!               the first entry in the specified output units.
        IF( NREC .GT. MXINDAT ) THEN
            UIDX = 1

!.............  If output records are the same as the input records, then
!               the units specified for each data column can be used.
        ELSE
            UIDX = OUTCOL

        END IF

!.............  Set temporary units from input units
        T_UNIT = I_UNIT

!.............  Set output units and conversion factor, if output units are
!               set by the report configuration file
        IF( ALLUSET( UIDX, RCNT ) .NE. ' ' ) THEN

            O_UNIT = ALLUSET( UIDX, RCNT )

!.................  Set the numerators and denominators...
!.................  Make sure if no denominator is given that there won't be
!                   a problem
!.................  For the input units:
            L2 = LEN_TRIM( T_UNIT )
            L  = INDEX( T_UNIT, '/' )
            LM1 = L - 1
            LP1 = L + 1
            IF( L .LE. 0 ) THEN
                LM1 = L2
                LP1 = 0
            END IF

            NUM_I = ADJUSTL( T_UNIT( 1:LM1 ) )
            IF( LP1 .GT. 0 ) DEN_I = ADJUSTL( T_UNIT( LP1:L2 ) )

!.................  If input denominator is hourly, but reporting is not
!                   hourly, then sum per day.  Change denominator accordingly.
            IF( DEN_I .EQ. 'hr'  .AND.&
            &          .NOT. ALLRPT( RCNT )%BYHOUR ) THEN
                DEN_I = 'day'
                L2 = LEN_TRIM( NUM_I )
                T_UNIT = NUM_I( 1:L2 ) // '/' // DEN_I
            END IF

!.................  For the output units:
            L2 = LEN_TRIM( O_UNIT )
            L  = INDEX( O_UNIT, '/' )
            LM1 = L - 1
            LP1 = L + 1
            IF( L .LE. 0 ) THEN
                LM1 = L2
                LP1 = 0
            END IF

            NUM_O = ADJUSTL( O_UNIT( 1:LM1 ) )
            IF( LP1 .GT. 0 ) DEN_O = ADJUSTL( O_UNIT( LP1:L2 ) )

!.................  Get factor for the numerator and denominator
            FAC1 = UNITFAC( T_UNIT, O_UNIT, .TRUE. )
            FAC2 = UNITFAC( T_UNIT, O_UNIT, .FALSE. )

            IF( FAC1 .LT. 0. ) THEN
                NUM_O = NUM_I
                FAC1 = 1.
            ENDIF
            IF( FAC2 .LT. 0. ) THEN
                DEN_O = DEN_I
                FAC2 = 1.
            ENDIF

!.................  Set the final output units
            UNITS ( OUTCOL ) = TRIM( NUM_O ) // '/' // TRIM( DEN_O )
            FACTOR( OUTCOL ) = FAC1 / FAC2

!.............  If output units not set, then use input units
!.............  Set conversion factor to 1
        ELSE

            UNITS ( OUTCOL ) = T_UNIT
            FACTOR( OUTCOL ) = 1.

        END IF

!.............  Update output units with population normalization
        IF( RPT_%NORMPOP ) THEN
            UNITS( OUTCOL )= MULTUNIT( UNITS( OUTCOL ), '1/person' )
        END IF

!.................  Update output units with cell area normalization
        IF( RPT_%NORMCELL ) THEN
            UNITS( OUTCOL ) = MULTUNIT( UNITS( OUTCOL ), '1/m^2' )
        END IF


        RETURN

    END SUBROUTINE UPDATE_OUTUNIT

END SUBROUTINE REPUNITS


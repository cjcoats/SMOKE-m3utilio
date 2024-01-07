
SUBROUTINE GETHDR( MXDATA, CFLAG, YFLAG, DFLAG,    &
                   LINE, ICC, INY, NPOA, EOS )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine interprets the header lines from the inventory files
    !
    !  PRECONDITIONS REQUIRED:
    !      Line must be read in and left-justified. Valid pollutant and activity
    !      list must be populated. Country codes list must be populated. The
    !      values of ICC, INY, and NPOA must be initialized.
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       created by M. Houyoux (2/2000)
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

    !.......   MODULES for public variables
    !.......  This module contains the lists of unique inventory information
    USE MODLISTS, ONLY: MXIDAT, INVDUNT, INVDCNV, INVDNAM,    &
                        ITCASA, ITNAMA, NINVTBL, ITKEEPA

    !.......  This module contains the arrays for state and county summaries
    USE MODSTCY, ONLY: NCOUNTRY, CTRYNAM, CTRYCOD

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: DATPOS, TMPNAM

    IMPLICIT NONE

    !.......   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    !.......   NOTE that NDROP and EDROP are not used at present
    INTEGER     , INTENT(IN   ) :: MXDATA      ! max no data variables allowed
    LOGICAL     , INTENT(IN   ) :: CFLAG       ! true: country header expected
    LOGICAL     , INTENT(IN   ) :: YFLAG       ! true: year header expected
    LOGICAL     , INTENT(IN   ) :: DFLAG       ! true: data names header expected
    CHARACTER(*), INTENT(INOUT) :: LINE      ! full record in string
    INTEGER     , INTENT(  OUT) :: ICC         ! country code
    INTEGER     , INTENT(  OUT) :: INY         ! inventory year
    INTEGER     , INTENT(  OUT) :: NPOA        ! no. data names
    INTEGER     , INTENT(  OUT) :: EOS         ! error status

    !.......   EXTERNAL FUNCTIONS and their descriptions:
    INTEGER, EXTERNAL :: GETNLIST
    REAL,    EXTERNAL :: UNITFAC
    LOGICAL, EXTERNAL :: USEEXPGEO

    !.......   Local allocatable arrays
    CHARACTER(NAMLEN3), ALLOCATABLE :: UNITS( : )

    !.......   Other local variables
    INTEGER         I, J, K, L, L1, L2, V      ! counters and indices

    INTEGER         COD           !  tmp data variable position in valid list
    INTEGER         IOS           !  i/o status
    INTEGER         NFINAL        !  final count of all "kept" pollutants
    INTEGER      :: NUNIT = 0     !  i/o status

    LOGICAL, SAVE:: ACT_FLAG  = .FALSE.     ! true: activities in data file
    LOGICAL, SAVE:: FIRSTIME  = .TRUE.      ! first time subroutine is called
    LOGICAL, SAVE:: UNITFLAG  = .FALSE.     ! true: units have already been set

    REAL            DFAC        !  denominator units conversion factor
    REAL            NFAC        !  numerator units conversion factor

    CHARACTER       CBUF        !  single char buffer
    CHARACTER(20)   CNTRY       !  country name
    CHARACTER(300)  MESG        !  message buffer
    CHARACTER(600)  BUFFER      !  tmp line buffer, upper case

    CHARACTER(NAMLEN3) CVAR          ! tmp variable name
    CHARACTER(NAMLEN3) DBUF          ! tmp units denominator
    CHARACTER(NAMLEN3) NBUF          ! tmp units numerator
    CHARACTER(NAMLEN3) UBUF          ! tmp units

    CHARACTER(16), PARAMETER :: PROGNAME = 'GETHDR'     ! Program name

    !***********************************************************************
    !   begin body of subroutine GETHDR

    !.......  First time the routine is called, initialize values
    IF( FIRSTIME ) THEN

        IF( CFLAG ) ICC  = -1
        IF( YFLAG ) INY  = 0
        IF( DFLAG ) NPOA = 0
        FIRSTIME = .FALSE.

    END IF

    !.......  Initialize i/o status
    EOS = 0

    !.......  Scan for header lines
    IF( LINE( 1:1 ) .EQ. CINVHDR ) THEN
        L1 = INDEX( LINE, '=' )
        L2 = LEN_TRIM( LINE )

        !.......  Convert to upper case
        BUFFER = LINE
        CALL UPCASE( BUFFER )

        !.......  Check for inventory type
        IF ( BUFFER(2:5) .EQ. 'TYPE' ) THEN

        !.......  Try to find activity in name of data, otherwise, assume emissions
            I = INDEX( BUFFER, 'ACTIVITY' )
            IF( I .GT. 0 ) THEN
                ACT_FLAG = .TRUE.
            END IF

        !.......  Check for country name
        ELSE IF ( BUFFER(2:8) .EQ. 'COUNTRY' ) THEN

            !.......  Skip country header if using expanded geographic codes
            IF ( USEEXPGEO() ) RETURN

            IF( L1 < 1 ) L1 = 8

            CNTRY = ADJUSTL( LINE( L1+1:L2 ) )
            I   = INDEX1( CNTRY, NCOUNTRY, CTRYNAM )

            IF( NCOUNTRY .LE. 0 ) THEN
                EOS = 2
                MESG = 'INTERNAL ERROR: Valid country list is uninitialized.'
                CALL M3MSG2( MESG )

            ELSE IF ( I .LE. 0 ) THEN
                EOS = 1
                MESG = 'ERROR: Country name "' // TRIM( CNTRY ) //    &
                       '" from inventory is not in country names file.'
                CALL M3MSG2( MESG )

            END IF

            ICC   = STR2INT( CTRYCOD( I ) ) / 100000

        !.......  Check for inventory year
        ELSE IF ( BUFFER(2:5) .EQ. 'YEAR' ) THEN

            IF( L1 < 1 ) L1 = 5

            INY = STR2INT( BUFFER( L1+1:L2 ) )

            IF ( INY .LT. 1971 .OR. INY .GT. MXIYEAR ) THEN
                WRITE( MESG, 94010 ) 'WARNING: Inventory year', INY,    &
                       'is outside of expected range.'
                CALL M3MSG2( MESG )
            END IF

        !.......  Check for units field to be used and compare to valid units
        !         from master input list
        !.......  Compute conversion factors for adjusting emissions in reader
        !         routines
        ELSE IF ( BUFFER(2:6) .EQ. 'UNITS' ) THEN         ! read in units

        !.......  Get the number of units fields and allocate memory
            BUFFER = BUFFER( 7:L2 )
            L = LEN_TRIM( BUFFER )
            NUNIT = GETNLIST( L, BUFFER )

            IF( UNITFLAG ) THEN
                MESG = 'ERROR: UNITS header encountered again in '//    &
                       'input. This is not allowed.'
                CALL M3MSG2( MESG )
                EOS = 5
                RETURN

            ELSE IF( .NOT. ALLOCATED( DATPOS ) ) THEN
                MESG = 'ERROR: UNITS header must be after DATA or POLID field.'
                CALL M3MSG2( MESG )
                EOS = 5
                RETURN

            ELSE IF( .NOT. ACT_FLAG ) THEN
                MESG = 'ERROR: UNITS header only valid for activity data.'
                CALL M3MSG2( MESG )
                EOS = 5
                RETURN

            ELSE IF( NUNIT .NE. NPOA ) THEN

                WRITE( MESG,94010 ) 'ERROR: number of UNITS ' //    &
                       'header entries (', NUNIT, ') is not ' //    &
                       CRLF() // BLANK10 // 'consistent with ' //    &
                       'number of DATA or POLID entries (', NPOA,    &
                       ')'
                CALL M3MSG2( MESG )
                EOS = 5
                RETURN

            END IF

            !.......  Allocate memory for local units field
            ALLOCATE( UNITS( NUNIT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'UNITS', PROGNAME )

            !.......  Parse the header line into the unit names
            CALL PARSLINE( LINE( 7:L2 ), NUNIT, UNITS )

            !.......  Post-process units and store in final arrays
            K = 0
            DO V = 1, NPOA

                !.......  Get position in master pollutants/activities list
                J  = DATPOS( V )

                !.......  Convert input units to valid SMOKE output units
                CALL UNITMATCH( UNITS( V ) )

                !.......  Set conversion factors for numerator and denominator,
                !                       if any
                NFAC = UNITFAC( UNITS( V ), INVDUNT( J ), .TRUE.  )
                DFAC = UNITFAC( UNITS( V ), INVDUNT( J ), .FALSE. )

                IF( NFAC .LT. 0. ) NFAC = 1.
                IF( DFAC .LT. 0. ) DFAC = 1.

                !.......  Error if units are not the same and the conversion
                !                       factor is equal to one
                IF( UNITS( V ) .NE. INVDUNT( J ) .AND.    &
                    NFAC .EQ. 1. .AND. DFAC .EQ. 1.      ) THEN

                    MESG = 'ERROR: Units conversion could not be '//    &
                           'made - use different input units.'
                    CALL M3MSG2( MESG )
                    EOS = 5

                    RETURN

                END IF

                !.......  Store current values of conversion factor
                INVDCNV( J ) = NFAC * DFAC

            END DO                   ! No. variables

            DEALLOCATE( UNITS )

            UNITFLAG = .TRUE.

        ELSE IF ( BUFFER(2:6) .EQ. 'POLID' .OR.    &
                  BUFFER(2:6) .EQ. 'DATA '       ) THEN         ! read in data names

            !....... Deallocate names for pollutant, if must
            IF( ALLOCATED( TMPNAM ) ) DEALLOCATE( TMPNAM,DATPOS )

            !.......  Allocate memory for current file for reading pol names
            !                   and storing positions in master list
            BUFFER = BUFFER( 7:L2 )
            L = LEN_TRIM( BUFFER )
            NPOA = GETNLIST( L, BUFFER )

            ALLOCATE( TMPNAM( NPOA ),    &
                      DATPOS( NPOA ), STAT=IOS )
            CALL CHECKMEM( IOS, 'DATPOS', PROGNAME )
            TMPNAM = ' '
            DATPOS = 0

            !.......  Set special value of error status when the number of
            !                   data variables attempted is greater than max allowed
            IF( NPOA .GT. MXDATA .AND. DFLAG ) THEN
                EOS = 4
            END IF

            !.......  Parse the header line into the pollutant names
            CALL PARSLINE( BUFFER, NPOA, TMPNAM )

            !.......  Store the position in master list of each pollutant
            !.......  Write error if pollutant is not found.
            NFINAL = 0
            DO V = 1, NPOA

                CVAR = TMPNAM( V )

                !.......  Look for variable in SMOKE variable names
                COD = INDEX1( CVAR, MXIDAT, INVDNAM )
                IF( COD .LE. 0 ) THEN

                    !......  Look for variable in Inventory Pollutant codes
                    COD = INDEX1( CVAR, NINVTBL, ITCASA )

                    !.......   If pollutant is not "kept", then take it
                    !                          out of the count and the list
                    IF( COD .GT. 0 ) THEN
                        IF( .NOT. ITKEEPA( COD ) ) CYCLE
                        COD = INDEX1( ITNAMA(COD), MXIDAT, INVDNAM )
                        DATPOS( V ) = COD

                        NFINAL = NFINAL + 1
                        TMPNAM( NFINAL ) = CVAR

                    !.......  If not found in list of names or codes, then error
                    ELSE
                        EOS = 1
                        MESG = 'ERROR: Data variable "' // TRIM( CVAR )//  &
                               '" not in master data variable list'
                        CALL M3MSG2( MESG )
                    END IF

                !.......  Variable found in SMOKE names
                ELSE
                    DATPOS( V ) = COD
                    NFINAL = NFINAL + 1

                END IF

            END DO

            !.......  Reset NPOA with final count that drops unkept variables
            NPOA = NFINAL

        END IF

        !.......  If this is a header line, return here without checking that
        !         all of the valid header fields are set properly
        RETURN

    END IF

    EOS = -1

    !.......  If the line is not a header line, make sure that all of the
    !         important header lines have been read in...

    !.......  Check for country header
    IF( CFLAG .AND. ICC .LT. 0 ) THEN
        EOS = 1
        ICC = 11             ! to turn off error message
        MESG = 'ERROR: Country name was not set with ' //    &
               '#COUNTRY header before first data line.'
        CALL M3MSG2( MESG )
    END IF

    !.......  Check for inventory year header
    IF( YFLAG .AND. INY .LE. 0 ) THEN
        EOS = 1
        INY = 1           ! to turn off error message
        MESG = 'ERROR: Inventory year was not set with ' //    &
               '#YEAR header before first data line.'
        CALL M3MSG2( MESG )
    END IF

    !.......  Check for pollutant names header
    !.......  This check is coordinated with initialized value
    IF( DFLAG .AND. NPOA .EQ. 0 ) THEN
        EOS = 1
        NPOA = -1           ! to turn off error message
        MESG = 'ERROR: Data variable names were not set with ' //    &
               '#POLID or #DATA header before first data line.'
        CALL M3MSG2( MESG )
    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE GETHDR


SUBROUTINE WRMRGREP( JDATE, JTIME, NIDX )

!***********************************************************************
!  subroutine WRMRGREP body starts at line
!
!  DESCRIPTION:
!      The purpose of this subroutine is to write state and county totals
!      with headers and formatted nicely by column width.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 8/99 by M. Houyoux
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
!****************************************************************************
    USE M3UTILIO

!.........  MODULES for public variables
!.........  This module contains the major data structure and control flags
    USE MODMERGE, ONLY: SDATE, STIME, EDATE, ETIME, TSTEP,&
    &                    AFLAG, BFLAG, MFLAG, PFLAG, XFLAG, SFLAG,&
    &                    ANMSPC, BNMSPC, MNMSPC, PNMSPC, NMSPC,&
    &                    AEMNAM, BEMNAM, MEMNAM, PEMNAM,&
    &                    ANIPOL, MNIPPA, PNIPOL, NIPPA,&
    &                    AEINAM, MEANAM, PEINAM,&
    &                    AEBCNY, BEBCNY, MEBCNY, PEBCNY, TEBCNY,&
    &                    AEBSTA, BEBSTA, MEBSTA, PEBSTA, TEBSTA,&
    &                    ARDEV,  BRDEV,  MRDEV,  PRDEV,  TRDEV,&
    &                    AUFLAG,         MUFLAG, PUFLAG, TUFLAG,&
    &                    ARFLAG,         MRFLAG, PRFLAG, TRFLAG,&
    &                    AECCNY,         MECCNY, PECCNY, TECCNY,&
    &                    AECSTA,         MECSTA, PECSTA, TECSTA,&
    &                    VGRPCNT, LREPSTA, LREPCTL, LREPCNY,&
    &                    LAVEDAY, EMNAM, EANAM, TOTUNIT,&
    &                    SIINDEX, SPINDEX

!.........  This module contains the arrays for state and county summaries
    USE MODSTCY, ONLY: NCOUNTY, NSTATE, STATNAM, CNTYNAM, CNTYCOD

    IMPLICIT NONE

!...........   INCLUDES:

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters\

!...........   Subroutine arguments
    INTEGER, INTENT (IN) :: JDATE  ! julian date  (YYYYDDD)
    INTEGER, INTENT (IN) :: JTIME  ! time (HHMMSS)
    INTEGER, INTENT (IN) :: NIDX   ! group index

!...........   Local allocatable arrays
    LOGICAL, ALLOCATABLE :: LUPDATE( : ) ! true: units/names not yet updated

    CHARACTER(NAMLEN3), ALLOCATABLE, SAVE ::  NAMES( : )
    CHARACTER(NAMLEN3), ALLOCATABLE, SAVE :: ANAMES( : )
    CHARACTER(NAMLEN3), ALLOCATABLE, SAVE :: BNAMES( : )
    CHARACTER(NAMLEN3), ALLOCATABLE, SAVE :: MNAMES( : )
    CHARACTER(NAMLEN3), ALLOCATABLE, SAVE :: PNAMES( : )

    CHARACTER(NAMLEN3), ALLOCATABLE, SAVE ::  UNITS( : )
    CHARACTER(NAMLEN3), ALLOCATABLE, SAVE :: AUNITS( : )
    CHARACTER(NAMLEN3), ALLOCATABLE, SAVE :: BUNITS( : )
    CHARACTER(NAMLEN3), ALLOCATABLE, SAVE :: MUNITS( : )
    CHARACTER(NAMLEN3), ALLOCATABLE, SAVE :: PUNITS( : )

!...........   Local group counts
    INTEGER, SAVE :: ACNT = 0       ! area source output vars count
    INTEGER, SAVE :: BCNT = 0       ! biogenic source output vars count
    INTEGER, SAVE :: MCNT = 0       ! mobile source output vars count
    INTEGER, SAVE :: PCNT = 0       ! point source output vars count
    INTEGER, SAVE :: TCNT = 0       ! total source output vars count

!...........   Other local variables

    INTEGER          C, F, I, J, K, L, L2, N, V  ! counters and indices
    INTEGER          IOS            ! i/o status
    INTEGER       :: KA = 0         ! tmp search indices
    INTEGER       :: KB = 0         ! tmp search indices
    INTEGER       :: KM = 0         ! tmp search indices
    INTEGER       :: KP = 0         ! tmp search indices
    INTEGER, SAVE :: MAXCYWID       ! max width of county names
    INTEGER, SAVE :: MAXSTWID       ! max width of state names
    INTEGER, SAVE :: NC             ! tmp no. counties in domain
    INTEGER, SAVE :: NS             ! tmp no. states in domain
    INTEGER, SAVE :: PGRP   = 0     ! group from previous call
    INTEGER, SAVE :: PDATE          ! date for computing PTIME
    INTEGER, SAVE :: PTIME          ! start time of period
    INTEGER, SAVE :: REPTIME        ! report time
    INTEGER, SAVE :: NVPGP          ! no. variables per group

    REAL   , SAVE :: UNITFAC        ! units conversion factor for reports

    LOGICAL, SAVE :: FIRSTIME= .TRUE. ! true: first time routine called
    LOGICAL, SAVE :: WARNON  = .TRUE. ! true: report warning

    CHARACTER(3000)    DATFMT     ! format for data
    CHARACTER(3000)    HDRFMT     ! format for header
    CHARACTER(3000)    HEADER     ! header for output files
    CHARACTER(3000)    LINFLD     ! line of dashes
    CHARACTER(300)     MESG       ! message buffer
    CHARACTER(NAMLEN3) SBUF       ! tmp pol or species name
    CHARACTER(NAMLEN3) CBUF       ! tmp units field

    CHARACTER(16) :: PROGNAME = 'WRMRGREP' ! program name

!***********************************************************************
!   begin body of subroutine WRMRGREP

!.........    If first time the routine is called
    IF( FIRSTIME ) THEN

        NC = NCOUNTY
        NS = NSTATE

!.............  Get time for output of reports
        MESG = 'Report time for country, state, and county totals'
        REPTIME = ENVINT( 'SMK_REPORT_TIME', MESG, 230000, IOS )

!............. NOTE that this would be a good place for an environment variable
!              that permits users to select the write format of the data

!.............  Get the maximum width for the state names
        MAXSTWID = 0
        DO I = 1, NS
            L = LEN_TRIM( STATNAM( I ) )
            MAXSTWID = MAX( MAXSTWID, L )
        END DO

!.............  Get the maximum width for the county names
        MAXCYWID = 0
        DO I = 1, NC
            L = LEN_TRIM( CNTYNAM( I ) )
            MAXCYWID = MAX( MAXCYWID, L )
        END DO

        PDATE = SDATE
        PTIME = STIME

!.............  Allocate memory for the units and names based on the
!               largest group size
!.............  Allocate memory for the logical update names/units flag
        J = MAX( NIPPA, NMSPC )
        I = MAXVAL( VGRPCNT )
        ALLOCATE( NAMES( I ),&
        &          UNITS( I ),&
        &         ANAMES( I ),&
        &         AUNITS( I ),&
        &         BNAMES( I ),&
        &         BUNITS( I ),&
        &         MNAMES( I ),&
        &         MUNITS( I ),&
        &         PNAMES( I ),&
        &         PUNITS( I ),&
        &         LUPDATE( J ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LUPDATE', PROGNAME )

        FIRSTIME = .FALSE.

    END IF

!.........  If current group is a new group...
    IF( NIDX .NE. PGRP ) THEN

!.............  Initialize units and names
        UNITS  = ' '  ! array
        AUNITS = ' '  ! array
        BUNITS = ' '  ! array
        MUNITS = ' '  ! array
        PUNITS = ' '  ! array
        ANAMES = ' '  ! array
        BNAMES = ' '  ! array
        MNAMES = ' '  ! array
        PNAMES = ' '  ! array

!.............  Initialize logical update states
        LUPDATE = .TRUE.   ! array

!.............  Set no. variables per group. Recall that this number can include
!               multiple species-pollutant combos (e.g., multiple processes for
!               mobile sources).
        NVPGP = VGRPCNT( NIDX )

!.............  Create units labels from variable units for current group
        ACNT = 0
        BCNT = 0
        MCNT = 0
        PCNT = 0
        TCNT = 0
        DO V = 1, NVPGP

!.................  Get indices
            I = SIINDEX( V,NIDX )  ! index to EANAM
            J = I
            SBUF = EANAM( J )

            IF( SFLAG ) THEN
                J = SPINDEX( V,NIDX )  ! index to EMNAM
                SBUF = EMNAM( J )
            END IF

!.................  If this species or pollutant has not already been
!                   encountered, update the source-category-specific arrays
!                   so that they will be in the same order as the global lists
            IF( LUPDATE( J ) ) THEN

                LUPDATE( J ) = .FALSE.

!.....................  Set names and units for totals output
                IF( SFLAG ) THEN
                    L = LEN_TRIM( TOTUNIT( J ) )
                    CBUF = '[' // TOTUNIT( J )( 1:L ) // ']'
                ELSE
                    L = LEN_TRIM( TOTUNIT( I ) )
                    CBUF = '[' // TOTUNIT( I )( 1:L ) // ']'
                END IF

                TCNT = TCNT + 1
                NAMES( TCNT ) = SBUF
                UNITS( TCNT ) = CBUF

!.....................  Get indices for source categories, depending on
!                       speciation or not
                IF( SFLAG ) THEN
                    IF( AFLAG ) KA = INDEX1( SBUF, ANMSPC, AEMNAM )
                    IF( BFLAG ) KB = INDEX1( SBUF, BNMSPC, BEMNAM )
                    IF( MFLAG ) KM = INDEX1( SBUF, MNMSPC, MEMNAM )
                    IF( PFLAG ) KP = INDEX1( SBUF, PNMSPC, PEMNAM )
                ELSE
                    IF( AFLAG ) KA = INDEX1( SBUF, ANIPOL, AEINAM )
                    KB = 0
                    IF( MFLAG ) KM = INDEX1( SBUF, MNIPPA, MEANAM )
                    IF( PFLAG ) KP = INDEX1( SBUF, PNIPOL, PEINAM )
                END IF

!.....................  Set names and units for area sources
                IF( KA .GT. 0 ) THEN
                    ACNT = ACNT + 1
                    ANAMES( ACNT ) = SBUF
                    AUNITS( ACNT ) = CBUF
                END IF

!.....................  Set names and units for biogenics
                IF( KB .GT. 0 ) THEN
                    BCNT = BCNT + 1
                    BNAMES( BCNT ) = SBUF
                    BUNITS( BCNT ) = CBUF
                END IF

!.....................  Set names and units for mobile sources
                IF( KM .GT. 0 ) THEN
                    MCNT = MCNT + 1
                    MNAMES( MCNT ) = SBUF
                    MUNITS( MCNT ) = CBUF
                END IF

!.....................  Set names and units for point sources
                IF( KP .GT. 0 ) THEN
                    PCNT = PCNT + 1
                    PNAMES( PCNT ) = SBUF
                    PUNITS( PCNT ) = CBUF
                END IF

            END IF  ! End if global units already defined or not

        END DO      ! End loop on group

        PGRP = NIDX

    END IF   ! End of processing for current group


!.........  Do not report if this time is not appropriate
    IF( JTIME .NE. REPTIME .AND.&
    &  ( JDATE .NE. EDATE .OR. JTIME .NE. ETIME ) ) RETURN

!.............  If required, create and write state totals, either controlled
!               or uncontrolled, depending on sector and which are controlled.
    IF ( LREPSTA ) THEN

!.............  Controlled area sources
        IF( ( AUFLAG .OR. ARFLAG ) .AND. LREPCTL ) THEN
            CALL CREATE_HEADER( 'Controlled area' )
            CALL CREATE_STATE( NC, NS, ACNT, AECCNY, AECSTA )
            CALL WRITE_STA( ARDEV, NS, ACNT, ANAMES, AUNITS, AECSTA)

!.................  Update state totals
            IF( XFLAG ) THEN
                CALL TOT_UPDATE( NS, ACNT, ANAMES, AECSTA, TECSTA)
            END IF

!.............  Uncontrolled area sources
        ELSE IF ( AFLAG ) THEN

            CALL CREATE_HEADER( 'Area' )
            CALL CREATE_STATE( NC, NS, ACNT, AEBCNY, AEBSTA )
            CALL WRITE_STA( ARDEV, NS, ACNT, ANAMES, AUNITS, AEBSTA)

!.....................  Update state totals
            IF( XFLAG ) THEN
                CALL TOT_UPDATE( NS, ACNT, ANAMES, AEBSTA, TEBSTA)
            END IF

        END IF

!.............  Biogenic sources (speciated by definition)
        IF( BFLAG ) THEN

            CALL CREATE_HEADER( 'Biogenic' )
            CALL CREATE_STATE( NC, NS, BCNT, BEBCNY, BEBSTA )
            CALL WRITE_STA( BRDEV, NS, BCNT, BNAMES, BUNITS, BEBSTA)

!.................  Update state totals
            IF( XFLAG ) THEN
                CALL TOT_UPDATE( NS, BCNT, BNAMES, BEBSTA, TEBSTA )
            END IF

        END IF

!.............  Controlled mobile sources
        IF( ( MUFLAG .OR. MRFLAG ) .AND. LREPCTL ) THEN
            CALL CREATE_HEADER( 'Controlled mobile' )
            CALL CREATE_STATE( NC, NS, MCNT, MECCNY, MECSTA )
            CALL WRITE_STA( MRDEV, NS, MCNT, MNAMES, MUNITS, MECSTA)

!.................  Update state totals
            IF( XFLAG ) THEN
                CALL TOT_UPDATE( NS, MCNT, MNAMES, MECSTA, TECSTA)
            END IF

!.............  Uncontrolled mobile sources
        ELSE IF( MFLAG ) THEN

            CALL CREATE_HEADER( 'Mobile' )
            CALL CREATE_STATE( NC, NS, MCNT, MEBCNY, MEBSTA )
            CALL WRITE_STA( MRDEV, NS, MCNT, MNAMES, MUNITS, MEBSTA)

!.................  Update state totals
            IF( XFLAG ) THEN
                CALL TOT_UPDATE( NS, MCNT, MNAMES, MEBSTA, TEBSTA)
            END IF

        END IF

!.............  Controlled point sources
        IF( ( PUFLAG .OR. PRFLAG ) .AND. LREPCTL ) THEN
            CALL CREATE_HEADER( 'Controlled point' )
            CALL CREATE_STATE( NC, NS, PCNT, PECCNY, PECSTA )
            CALL WRITE_STA( PRDEV, NS, PCNT, PNAMES, PUNITS, PECSTA)

!.................  Update state totals
            IF( XFLAG ) THEN
                CALL TOT_UPDATE( NS, PCNT, PNAMES, PECSTA, TECSTA)
            END IF

!.............  Uncontrolled point sources
        ELSE IF( PFLAG ) THEN

            CALL CREATE_HEADER( 'Point' )
            CALL CREATE_STATE( NC, NS, PCNT, PEBCNY, PEBSTA )
            CALL WRITE_STA( PRDEV, NS, PCNT, PNAMES, PUNITS, PEBSTA)

!.................  Update state totals
            IF( XFLAG ) THEN
                CALL TOT_UPDATE( NS, PCNT, PNAMES, PEBSTA, TEBSTA )
            END IF
        END IF

!.............  Controlled total sources
        IF( ( TUFLAG .OR. TRFLAG ) .AND. LREPCTL ) THEN
            CALL CREATE_HEADER( 'Controlled total' )
            CALL CREATE_STATE( NC, NS, TCNT, TECCNY, TECSTA )
            CALL WRITE_STA( TRDEV, NS, TCNT, NAMES, UNITS, TECSTA)

!.............  Uncontrolled total sources
        ELSE IF( XFLAG ) THEN

            CALL CREATE_HEADER( 'Total' )
            CALL WRITE_STA( TRDEV, NS, TCNT, NAMES, UNITS, TEBSTA )

        END IF

    END IF

!.........  If required, write county totals
    IF( LREPCNY ) THEN

!.............  Give warning if not already given and if reporting of controlled
!               emissions was selected
        IF ( LREPCTL .AND. WARNON ) THEN
            WARNON = .FALSE.

            MESG = 'WARNING: County reports do not support '//&
            &       'including controlled emissions, but this '//&
            &       'combination was requested.'
            CALL M3MESG( MESG )

        END IF

!.............  Area sources
        IF( AFLAG ) THEN

            CALL CREATE_HEADER( 'Area' )
            CALL WRITE_CNY( ARDEV, NC, ACNT, ANAMES, AUNITS, AEBCNY)

!.................  Update state totals
            IF( XFLAG ) THEN
                CALL TOT_UPDATE( NC, ACNT, ANAMES, AEBCNY, TEBCNY)
            END IF
        END IF

!.............  Biogenc sources (speciated by definition)
        IF( BFLAG ) THEN

            CALL CREATE_HEADER( 'Biogenic' )
            CALL WRITE_CNY( BRDEV, NC, BCNT, BNAMES, BUNITS, BEBCNY)

!.................  Update state totals
            IF( XFLAG ) THEN
                CALL TOT_UPDATE( NC, BCNT, BNAMES, BEBCNY, TEBCNY )
            END IF

        END IF

!.............  Mobile sources
        IF( MFLAG ) THEN

            CALL CREATE_HEADER( 'Mobile' )
            CALL WRITE_CNY( MRDEV, NC, MCNT, MNAMES, MUNITS, MEBCNY)

!.................  Update state totals
            IF( XFLAG ) THEN
                CALL TOT_UPDATE( NC, MCNT, MNAMES, MEBCNY, TEBCNY)
            END IF

        END IF

!.............  Point sources
        IF( PFLAG ) THEN

            CALL CREATE_HEADER( 'Point' )
            CALL WRITE_CNY( PRDEV, NC, PCNT, PNAMES, PUNITS, PEBCNY)

!.................  Update state totals
            IF( XFLAG ) THEN
                CALL TOT_UPDATE( NC, PCNT, PNAMES, PEBCNY, TEBCNY)
            END IF

        END IF

!.............  Combined sources
        IF( XFLAG ) THEN

            CALL CREATE_HEADER( 'Total' )
            CALL WRITE_CNY( TRDEV, NC, TCNT, NAMES, UNITS, TEBCNY )

        END IF

    END IF

!.........  Intialize summed emissions to zero
    CALL INITSTCY

!........  After output, increment date and time one step to set the start
!          of the next period
    PDATE = JDATE
    PTIME = JTIME
    CALL NEXTIME( PDATE, PTIME, TSTEP )  ! advance one

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats.............94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

!*****************  INTERNAL SUBPROGRAMS   *****************************

CONTAINS

    SUBROUTINE CREATE_STATE( NC, NS, NDIM, CY_EMIS, ST_EMIS )

!.............  Subprogram arguments
        INTEGER     , INTENT (IN) :: NC
        INTEGER     , INTENT (IN) :: NS
        INTEGER     , INTENT (IN) :: NDIM
        REAL        , INTENT (IN) :: CY_EMIS( NC, NDIM )
        REAL        , INTENT(OUT) :: ST_EMIS( NS, NDIM )

!.............  Local variables
        INTEGER  I, J, N

        CHARACTER(FIPLEN3) PSTA, STA

!..............................................................................

        PSTA = ' '
        N = 0
        DO I = 1, NC


            STA = CNTYCOD( I )
            STA( FIPLEN3-2:FIPLEN3 ) = '000'
            IF( STA .NE. PSTA ) THEN
                N = N + 1
                PSTA = STA
            END IF

            DO J = 1, NDIM
                ST_EMIS( N,J ) = ST_EMIS( N,J ) + CY_EMIS( I,J )
            END DO

        END DO

        RETURN

    END SUBROUTINE CREATE_STATE

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

    SUBROUTINE CREATE_HEADER( CATNAME )

!.............  MODULES for public variables
!.............  This module contains the global variables for the 3-d grid
        USE MODGRID, ONLY: GRDNM

!.............  Subprogram arguments
        CHARACTER(*), INTENT (IN) :: CATNAME

!.............  Local variables
        INTEGER   K1, K2, L, LD1, LD2

        CHARACTER(10) TYPENAM
        CHARACTER(15) DATANAM

!..............................................................................

        HEADER = '# ' // CATNAME // ' source'
        L  = LEN_TRIM( HEADER )

        DATANAM = ' average'
        IF( LAVEDAY .AND. ( AFLAG .OR. PFLAG ) )&
        &    DATANAM = ' average day'
        LD1 = LEN_TRIM( DATANAM )

        TYPENAM = ' inventory'
        IF( SFLAG ) TYPENAM = ' speciated'

        HEADER = HEADER( 1:L ) // DATANAM( 1:LD1 ) //&
        &         TYPENAM // ' emissions'
        L = LEN_TRIM( HEADER )

        IF( JDATE .NE. 0 ) THEN

            K1  = WKDAY( PDATE )
            K2  = WKDAY( JDATE )
            LD1 = LEN_TRIM( DAYS( K1 ) )
            LD2 = LEN_TRIM( DAYS( K2 ) )

            WRITE( HEADER,94010 ) HEADER( 1:L ) // ' from ' //&
            &   CRLF() // '#' // BLANK5(1:4) //&
            &   DAYS( K1 )( 1:LD1 ) // ' ' // MMDDYY( PDATE ) //&
            &   ' at', PTIME, 'to' // CRLF() // '#' // BLANK5(1:4) //&
            &   DAYS( K2 )( 1:LD2 ) // ' '// MMDDYY( JDATE ) //&
            &   ' at', JTIME
            L = LEN_TRIM( HEADER )

        END IF

        HEADER = HEADER( 1:L ) // ' within grid ' // GRDNM

        RETURN

!.......................... FORMAT STATEMENTS ................................

94010   FORMAT( 10( A, :, I8, :, 1X ) )

    END SUBROUTINE CREATE_HEADER

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

    SUBROUTINE CREATE_FORMATS( NDIM, MAXCOL1, INNAMS, INUNIT,&
    &                           WIDTHS, OUTNAMS, OUTUNIT  )

!.............  Subprogram arguments
        INTEGER     , INTENT (IN)     :: NDIM
        INTEGER     , INTENT (IN)     :: MAXCOL1
        CHARACTER(*), INTENT (IN)     :: INNAMS ( NDIM )
        CHARACTER(*), INTENT (IN)     :: INUNIT ( NDIM )
        INTEGER     , INTENT (IN OUT) :: WIDTHS ( 0:NDIM )
        CHARACTER(*), INTENT (IN OUT) :: OUTNAMS( NDIM )
        CHARACTER(*), INTENT (IN OUT) :: OUTUNIT( NDIM )

!.............  Local parameters
        INTEGER, PARAMETER :: EFMTWID = 10  ! minimum width of fields
        INTEGER, PARAMETER :: EFMTDEC = 4   ! number of decimal places

!.............  Local variables
        INTEGER       I1, I2, J, L, L1, L2
        CHARACTER(30)  :: SPACE = ' '
        CHARACTER(3000) :: TMPFMT

!.............................................................................

!.............  Adjust maximum width of numbers in case width of variable
!               names or units is greater than width of numbers
!.............  Also, move the position of the names and units so that output
!               strings will look right-justified in the file.
        DO J = 1, NDIM
            L1 = LEN_TRIM( INNAMS( J ) )
            WIDTHS ( J ) = MAX( WIDTHS( J ), L1, EFMTWID )
            L2 = LEN_TRIM( INUNIT( J ) )
            WIDTHS ( J ) = MAX( WIDTHS( J ), L2 )

            I1 = WIDTHS( J ) - L1
            IF( I1 .GT. 0 ) THEN
                WRITE( OUTNAMS( J ), '(A,A)' ) SPACE( 1:I1 ),&
                &                               INNAMS( J )( 1:L1 )
            ELSE
                WRITE( OUTNAMS( J ), '(A)' ) INNAMS( J )( 1:L1 )
            END IF

            I2 = WIDTHS( J ) - L2
            IF( I2 .GT. 0 ) THEN
                WRITE( OUTUNIT( J ), '(A,A)' ) SPACE( 1:I2 ),&
                &                               INUNIT( J )( 1:L2 )
            ELSE
                WRITE( OUTUNIT( J ), '(A)' ) INUNIT( J )( 1:L2 )
            END IF

        END DO

!.............  Create format statement for output of header
        WIDTHS( 0 ) = MAXCOL1
        WRITE( HDRFMT, '( "(A",A)' ) ',";"'
        DO J = 1, NDIM
            TMPFMT = HDRFMT
            WRITE( HDRFMT, '(A, ",1X,A",I2.2,A)' )&
            &       TRIM( TMPFMT ), WIDTHS( J ), ',";"'
        END DO
        TMPFMT = HDRFMT
        WRITE( HDRFMT, '(A)' ) TRIM( TMPFMT ) // ')'

!.............  Create format statement for output of emissions
        WRITE( DATFMT, '( "(A",I2.2,A)' ) WIDTHS( 0 ), ',";"'
        DO J = 1, NDIM
            TMPFMT = DATFMT
            WRITE( DATFMT, '(A, ",1X,E",I2.2,".",I1,A)' )&
            &       TRIM( TMPFMT ), WIDTHS( J ), EFMTDEC, ',";"'
        END DO
        TMPFMT = DATFMT
        L = LEN_TRIM( TMPFMT )
        WRITE( DATFMT, '(A)' ) TMPFMT( 1:L ) // ')'

        RETURN

    END SUBROUTINE CREATE_FORMATS

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

    SUBROUTINE WRITE_STA( FDEV, NS, NDIM, VNAMES, INUNIT,&
    &                      ST_EMIS )

!.............  Subprogram arguments
        INTEGER     , INTENT (IN) :: FDEV
        INTEGER     , INTENT (IN) :: NS
        INTEGER     , INTENT (IN) :: NDIM
        CHARACTER(*), INTENT (IN) :: VNAMES ( NDIM )
        CHARACTER(*), INTENT (IN) :: INUNIT ( NDIM )
        REAL        , INTENT (IN) :: ST_EMIS( NS, NDIM )

!.............  Arrays allocated by subprogram argument
        INTEGER            MAXWID ( 0:NDIM )
        CHARACTER(NAMLEN3) OUTNAMS( NDIM )
        CHARACTER(NAMLEN3) OUTUNIT( NDIM )

!.............  Local variables
        INTEGER       I, J, L, L2

        REAL          VAL

        CHARACTER(20) :: HDRBUF  = '#'
        CHARACTER(20) :: STLABEL = '# State'
        CHARACTER(30)    BUFFER

!..............................................................................

!.............  Get maximum width of numbers and
        DO J = 1, NDIM

            VAL = MAXVAL( ST_EMIS( 1:NS, J ) )
            WRITE( BUFFER, '(F30.1)' ) VAL
            BUFFER = ADJUSTL( BUFFER )
            MAXWID( J ) = LEN_TRIM( BUFFER )

        END DO

!.............  Rearrange labels and units to be in order of master list, which
!               is the order that the emission values themselves will be in

!.............  Get column labels and formats
        CALL CREATE_FORMATS( NDIM, 6+MAXSTWID, VNAMES, INUNIT,&
        &                     MAXWID, OUTNAMS, OUTUNIT )

!.............  Create line format
!            L2 = SUM( MAXWID ) + NDIM
!            LINFLD = REPEAT( '-', L2 )

!.............  Write header for state totals
        WRITE( FDEV, '(A)' ) '# '
        WRITE( FDEV, '(A)' ) HEADER( 1:LEN_TRIM( HEADER ) )

!.............  Write line
!            WRITE( FDEV, '(A)' ) LINFLD( 1:L2 )

!.............  Write units for columns
        WRITE( FDEV, HDRFMT ) ADJUSTL( HDRBUF),&
        &                      ( OUTUNIT(J), J=1,NDIM )

!.............  Write column labels
        WRITE( FDEV, HDRFMT ) ADJUSTL( STLABEL ),&
        &                    ( OUTNAMS( J ), J=1, NDIM )

!.............  Write state total emissions
        DO I = 1, NSTATE

!.................  Build output format depending on data values
            CALL DYNAMIC_FORMATS( NSTATE, NDIM, I, ST_EMIS,&
            &                      MAXWID(1), DATFMT )

!.................  Write out state name and converted emissions
            WRITE( FDEV, DATFMT ) STATNAM( I ),&
            &                      ( ST_EMIS( I,J ), J=1, NDIM )
        END DO

        RETURN

    END SUBROUTINE WRITE_STA

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

    SUBROUTINE WRITE_CNY( FDEV, NC, NDIM, VNAMES, INUNIT,&
    &                      CY_EMIS )

!.............  Subprogram arguments
        INTEGER     , INTENT (IN) :: FDEV
        INTEGER     , INTENT (IN) :: NC
        INTEGER     , INTENT (IN) :: NDIM
        CHARACTER(*), INTENT (IN) :: VNAMES ( NDIM )
        CHARACTER(*), INTENT (IN) :: INUNIT ( NDIM )
        REAL        , INTENT (IN) :: CY_EMIS( NC, NDIM )

!.............  Arrays allocated by subprogram argument
        INTEGER            MAXWID ( 0:NDIM )
        CHARACTER(NAMLEN3) OUTNAMS( NDIM )
        CHARACTER(NAMLEN3) OUTUNIT( NDIM )

!.............  Local variables
        INTEGER       I, J, L, L1, L2, N

        REAL          VAL

        CHARACTER(FIPLEN3+8) CDATFIP
        CHARACTER(FIPLEN3) PSTA, STA
        CHARACTER(60) :: HDRBUF  = '#'
        CHARACTER(60) :: STLABEL = '# County'
        CHARACTER(30)    BUFFER

!..............................................................................

!.............  Get maximum width of numbers
        DO J = 1, NDIM

            VAL = MAXVAL( CY_EMIS( 1:NC, J ) )
            WRITE( BUFFER, '(F30.1)' ) VAL
            BUFFER = ADJUSTL( BUFFER )
            MAXWID( J ) = LEN_TRIM( BUFFER )

        END DO

!.............  Get column labels and formats
        CALL CREATE_FORMATS( NDIM, 26+MAXSTWID+MAXCYWID, VNAMES,&
        &                     INUNIT, MAXWID, OUTNAMS, OUTUNIT )

!.............  Create line format
!            L2 = SUM( MAXWID ) + NDIM
!            LINFLD = REPEAT( '-', L2 )

!.............  Write header for county totals
        WRITE( FDEV, '(A)' ) '# '
        WRITE( FDEV, '(A)' ) HEADER( 1:LEN_TRIM( HEADER ) )

!.............  Write units for columns
        WRITE( FDEV, HDRFMT ) ADJUSTL( HDRBUF),&
        &                      ( OUTUNIT(J), J=1,NDIM )

!.............  Write column labels
        WRITE( FDEV, HDRFMT ) ADJUSTL( STLABEL ),&
        &                      ( OUTNAMS( J ), J=1, NDIM )

!.............  Write line
!            WRITE( FDEV, '(A)' ) LINFLD( 1:L2 )

!.............  Write county total emissions
        PSTA = ' '
        N = 0
        DO I = 1, NC

            STA = CNTYCOD( I )
            STA( FIPLEN3-2:FIPLEN3 ) = '000'
            IF( STA .NE. PSTA ) THEN
                N = N + 1
                PSTA = STA
            END IF

!.................  Write out county name and converted emissions
            WRITE( CDATFIP, '(I7.7,1X,A)' ) JDATE, CNTYCOD( I )

!.................  Build output format depending on data values
            CALL DYNAMIC_FORMATS( NC, NDIM, I, CY_EMIS,&
            &                      MAXWID(1), DATFMT )

            WRITE( FDEV,DATFMT ) CDATFIP // ' '// STATNAM(N) //&
            &                     CNTYNAM(I),&
            &                     ( CY_EMIS( I,J ), J=1, NDIM )

        END DO

        RETURN

    END SUBROUTINE WRITE_CNY

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

    SUBROUTINE TOT_UPDATE( NU, NDIM, VNAMES, INEMIS, TOTEMIS )

!.............  Subprogram arguments
        INTEGER     , INTENT    (IN) :: NU               ! no. summing units
        INTEGER     , INTENT    (IN) :: NDIM           ! no. pols or species
        CHARACTER(*), INTENT    (IN) :: VNAMES ( NDIM )   ! pol or spc names
        REAL        , INTENT    (IN) :: INEMIS ( NU, NDIM ) ! component emis
        REAL        , INTENT(IN OUT) :: TOTEMIS( NU, * )    ! total emis

!.............  Local variables
        INTEGER       I, K, V

!..............................................................................

!.............  Abort function call if only one source category
        IF( .NOT. XFLAG ) RETURN

!.............  Loop through pollutants or species, and search for
!               local names in master lists to get index to totals array
        DO V = 1, NDIM

            IF( SFLAG ) THEN
                K = INDEX1( VNAMES( V ), NMSPC, EMNAM )
            ELSE
                K = INDEX1( VNAMES( V ), NIPPA, EANAM )
            END IF

!.................  If pollutant or species is found, then add component
!                   emissions to the total
            IF( K .GT. 0 ) THEN

                DO I = 1, NU

                    TOTEMIS( I,K ) = TOTEMIS( I,K ) + INEMIS( I,V )

                END DO

            END IF

        END DO

        RETURN

    END SUBROUTINE TOT_UPDATE

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

    SUBROUTINE DYNAMIC_FORMATS( N1, N2, STCNT, EMIS,&
    &                            WIDTH, DATFMT )

        INTEGER     , INTENT (IN)  :: N1
        INTEGER     , INTENT (IN)  :: N2
        INTEGER     , INTENT (IN)  :: STCNT      ! counter for state/county index
        REAL        , INTENT (IN)  :: EMIS( N1, N2 )  ! state/count emission
        INTEGER     , INTENT (IN)  :: WIDTH( N2 )
        CHARACTER(*), INTENT (OUT) :: DATFMT

!.............  Local variables
        INTEGER  I
        CHARACTER(5) :: FMT

!..............................................................................

!.............  Initialize format array
        DATFMT = '(A'

!............. Determine significant figures that we want to report.
!............. Use a minimum of 5 significant figures, and a minimum
!              of 1 decimal place.  If value is < 0.1, then use exponential
!              with 4 decimal places
        DO I = 1, N2

            IF ( EMIS( STCNT,I ) .EQ. 0. ) THEN
                WRITE( FMT, '(A,I2.2,A)' ) 'F',WIDTH(I),'.0'

            ELSE IF( EMIS( STCNT,I ) .GE. 1000. ) THEN
                WRITE( FMT, '(A,I2.2,A)' ) 'F',WIDTH(I),'.1'

            ELSE IF ( EMIS( STCNT,I ) .GE. 100. ) THEN
                WRITE( FMT, '(A,I2.2,A)' ) 'F',WIDTH(I),'.2'

            ELSE IF ( EMIS( STCNT,I ) .GE. 10. ) THEN
                WRITE( FMT, '(A,I2.2,A)' ) 'F',WIDTH(I),'.3'

            ELSE IF ( EMIS( STCNT,I ) .GE. 0.1 ) THEN
                WRITE( FMT, '(A,I2.2,A)' ) 'F',WIDTH(I),'.4'

            ELSE IF ( EMIS( STCNT,I ) .LT. 0.1 ) THEN
                WRITE( FMT, '(A,I2.2,A)' ) 'E',WIDTH(I),'.4'
            END IF

            DATFMT = TRIM( DATFMT ) // ',"; "' // FMT

        END DO

        DATFMT = TRIM( DATFMT ) // ')'

        RETURN

    END SUBROUTINE DYNAMIC_FORMATS

END SUBROUTINE WRMRGREP


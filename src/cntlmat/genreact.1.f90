
SUBROUTINE GENREACT( PYEAR, ENAME, RPOL, USEPOL )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine processes the reactivity data ...
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 3/99 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!************************************************************************
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
!*************************************************************************
    USE M3UTILIO

!.........  MODULES for public variables
!.........  This module contains the inventory arrays
    USE MODSOURC, ONLY: CSOURC, CSCC, CIFIP, CISIC, CLINK

!.........  This module contains the control packet data and control matrices
    USE MODCNTRL, ONLY: RMTXMASS, RMTXMOLE, PCRIDX, PCRREPEM,&
    &                    PCRPRJFC, PCRMKTPN, PCRCSCC, PCRSPROF,&
    &                    RPTDEV, EMREPREA, CSPFREA, PRJFCREA,&
    &                    MKTPNREA, CSCCREA

!.........  This module contains the speciation profiles
    USE MODSPRO, ONLY: NSPFUL, NSPROF, SPROFN, SPECID, MOLEFACT,&
    &                   MASSFACT, INPRF, MXSPFUL, SPCNAMES,&
    &                   IDXSPRO, IDXSSPEC, NSPECIES

!.........  This module contains the information about the source category
    USE MODINFO, ONLY: NIPOL, MXCHRS, EINAM, NSRC, CATDESC, BYEAR,&
    &                   NCHARS, SC_BEGP, SC_ENDP, JSCC, CATEGORY, CRL

    IMPLICIT NONE

!...........   INCLUDES

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
    INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions


!...........   SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN) :: PYEAR  ! projection year for reactivity
    CHARACTER(*), INTENT (IN) :: ENAME  ! emission inventory file name
    CHARACTER(NAMLEN3), INTENT (IN) :: RPOL! pol for procssing reactvity
    LOGICAL     , INTENT (IN) :: USEPOL( NIPOL ) ! true: pol valid for pkt

!...........   Local parameters
    INTEGER,       PARAMETER :: NHEADER  = 15
    CHARACTER(16), PARAMETER :: PROGNAME = 'GENREACT' ! program name
    CHARACTER(15), PARAMETER :: HEADERS( NHEADER ) =&
    &                        ( / 'Source         ',&
    &                            'Region         ',&
    &                            'Plant          ',&
    &                            'Char1          ',&
    &                            'Char2          ',&
    &                            'Char3          ',&
    &                            'Char4          ',&
    &                            'Link           ',&
    &                            'Base SCC       ',&
    &                            'Base Emis      ',&
    &                            'New Base Emis  ',&
    &                            'Proj Factor    ',&
    &                            'New SCC        ',&
    &                            'New Spc Profile',&
    &                            'Mkt Pen Rate   '  / )

!...........   Local arrays allocated by subroutine arguments
    INTEGER          MCWID( MXCHRS )      !  max characteristic width
    INTEGER          HINDX( NHEADER )     !  header label index
    INTEGER          HCWID( NHEADER )     !  header label widths
    LOGICAL          LF   ( MXCHRS )      !  true: column should be output
    CHARACTER(20)    CHARS( MXCHRS )      !  source fields for output

!.........  Allocatable arrays
    INTEGER, ALLOCATABLE :: ISREA( : )   ! reactivity control data table index
    REAL   , ALLOCATABLE :: EMIS ( : )   ! base inventory emissions

    CHARACTER(NAMLEN3), ALLOCATABLE :: MASSONAM( : ) ! mass output species names
    CHARACTER(NAMLEN3), ALLOCATABLE :: MOLEONAM( : ) ! mole output species names

!...........   Logical names and unit numbers
    INTEGER, SAVE :: PDEV         !  speciation profiles unit no.
    INTEGER       :: SDEV         !  supplement file unit no.

    CHARACTER(16) :: RNAME = 'IOAPI_DAT' ! logical name for reading pols
    CHARACTER(16)    SNAME   ! logical name for mass-based react. cntrl mat
    CHARACTER(16)    LNAME   ! logical name for mole-based react. cntrl mat

!...........   Other local variables

    INTEGER          I, J, K, L, N, P, S, V     ! counters and indices

    INTEGER          IDUM        ! dummy integer
    INTEGER          IOS         ! i/o error status
    INTEGER          ITBL        ! position in full table of current profile
    INTEGER          NC          ! tmp number of src chars
    INTEGER          NCOLS       ! no. of output columns in report
    INTEGER          NCOUT       ! no. of src chars to put in report
    INTEGER          NMSPC       ! number of model species
    INTEGER          NSREAC      ! number of srcs w/ reactivity controls
    INTEGER          NTBL        ! number of species in current profile
    INTEGER          RDEV          ! Report unit number

    REAL             EMREP       ! tmp replacement emissions
    REAL             YFAC        ! tmp yr conversion factor

    CHARACTER(SPNLEN3) SPCODE ! tmp speciation profile code

    LOGICAL       :: EFLAG    = .FALSE.  ! true: error has occurred
    LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! true: first call to subroutine
    LOGICAL       :: LFLAG    = .FALSE.  ! true: link will be included in report
    LOGICAL, SAVE :: LAVEDAY  = .FALSE.  ! true: use average day emissions
    LOGICAL, SAVE :: MASSOUT  = .FALSE.  ! true: output mass-based spc facs
    LOGICAL, SAVE :: MOLEOUT  = .FALSE.  ! true: output mole-based spc facs
    LOGICAL, SAVE :: PFLAG    = .FALSE.  ! true: point source processing

    CHARACTER(4)     OUTTYPE             ! speciation output type
    CHARACTER(NAMLEN3) VNAM              ! tmp average day var name
    CHARACTER(256)   BUFFER              ! string buffer for building output fmt
    CHARACTER(256)   HDRSTR              ! string for part of header line
    CHARACTER(256)   MESG                ! message buffer
    CHARACTER(256)   OUTFMT              ! output format for report


!***********************************************************************
!   begin body of subroutine GENREACT

    IF( FIRSTIME ) THEN

!.............  Get environment variables that control subroutine behavior
!.............  Retrieve the type of speciation outputs (mass,mole,or all)
        MESG = 'Type of speciation outputs to create'
        CALL ENVSTR( 'SPEC_OUTPUT', MESG, 'ALL', OUTTYPE, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SPEC_OUTPUT"', 2 )
        END IF

!.............  Get environment variables that control program behavior
        MESG = 'Use annual or average day emissions'
        LAVEDAY = ENVYN( 'SMK_AVEDAY_YN', MESG, .FALSE., IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK_AVEDAY_YN"', 2 )
        END IF

!.........  Open reports file
        RPTDEV( 3 ) = PROMPTFFILE(&
        &               'Enter logical name for REACTIVITY REPORT',&
        &              .FALSE., .TRUE., CRL // 'REACREP', PROGNAME )
        RDEV = RPTDEV( 3 )

!.............  Set flags that depend on the value of OUTTYPE
        CALL UPCASE( OUTTYPE )
        MASSOUT = ( INDEX( OUTTYPE, 'ALL' ) .GT. 0 )
        MOLEOUT = ( INDEX( OUTTYPE, 'ALL' ) .GT. 0 )
        MASSOUT = ( INDEX( OUTTYPE, 'MASS' ) .GT. 0 .OR. MASSOUT )
        MOLEOUT = ( INDEX( OUTTYPE, 'MOLE' ) .GT. 0 .OR. MOLEOUT )

        MESG = 'Speciation profiles needed to process ' //&
        &       'reactivity packet...'
        CALL M3MSG2( MESG )

        PDEV = PROMPTFFILE(&
        &       'Enter logical name for SPECIATION PROFILES file',&
        &       .TRUE., .TRUE., 'GSPRO', PROGNAME )

        CALL DSCSPROF( PDEV, NIPOL, EINAM )

!...........  Determine if source category is point sources, or other
        PFLAG = ASSOCIATED( CISIC )

!.............  Allocate memory for arrays of speciation tables and unique lists
!               using the maximum number of profile table entires per pollutant,
!               MXSPFUL, which is from module MODSPRO
        ALLOCATE( INPRF( MXSPFUL ),&
        &         SPECID( MXSPFUL ),&
        &       MOLEFACT( MXSPFUL ),&
        &       MASSFACT( MXSPFUL ),&
        &          ISREA( NSRC ),&
        &           EMIS( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMIS', PROGNAME )
        ISREA = 0
        EMIS  = 0.

!...........  Write output report header
        IF( RDEV .GT. 0 ) THEN
            WRITE(RDEV,93000) 'Processed as '// CATDESC// ' sources'
            WRITE(RDEV,93390) 'Base inventory year ', BYEAR
            IF ( PYEAR .GT. 0 ) THEN
                WRITE(RDEV,93390) 'Reactivity packet projected ' //&
                &                  'year ', PYEAR
            END IF

            WRITE( RDEV, 93000 )&
            &        'Controls applied with /REACTIVITY/ packet '//&
            &        'for pollutant "' // TRIM( RPOL ) // '".'

            IF( LAVEDAY ) THEN
                WRITE(RDEV,93000)'Average day data basis in report'
            ELSE
                WRITE(RDEV,93000)'Annual total data basis in report'
            END IF

        END IF

        FIRSTIME = .FALSE.

    END IF

!.........  For pollutant in subroutine argument...

!.........  Determine which reactivity packet goes to each source
    CALL ASGNCNTL( NSRC, 1, 'REACTIVITY', USEPOL, RPOL,&
    &               IDUM, ISREA )

! NOTE:
!.........  Please note that the indexing for MCWID, HCWID, and HINDX
!           is confusing. The main confusing aspect is that the
!           +1's are added to index past the "Source" column, which
!           is not a part of the MCWID (which is just for the
!           source characteristics).  The +2's are added to index
!           past the "Source" column and to the column after the
!           last NCOUT column (last column of source chars).

!.........  Initialize valid columns
    LF = .FALSE.  ! array
    DO I = 1, NCHARS
        LF( I ) = .TRUE.
    END DO

!.........  Count up the number of sources that have reactivity data
    NSREAC = 0
    N      = 0
    MCWID   = 0  ! array
    DO S = 1, NSRC

        IF( ISREA( S ) .GT. 0 ) THEN
            NSREAC = NSREAC + 1

!............... Determine maximum width of report columns for source info
            CALL PARSCSRC( CSOURC( S ), MXCHRS, SC_BEGP, SC_ENDP,&
            &               LF, NC, CHARS )

            DO J = 1, NC
                MCWID( J ) = MAX( LEN_TRIM( CHARS(J) ), MCWID(J) )
                IF( MCWID( J ) .GT. 0 ) N = MAX( N, J )
            END DO

        END IF

    END DO

!........  If SCC is not included in source chars, then add it to report
    IF( JSCC .EQ. 0 ) THEN
        N = N + 1
        MCWID( N ) = SCCLEN3
    END IF

    NCOUT = N

!.........  Determine source characteristic header indices
    HINDX = 0    ! array
    SELECT CASE( CATEGORY )

!........  Area source
      CASE( 'AREA' )
        NCOLS = 3
        HINDX( 1 ) = 1
        HINDX( 2 ) = 2
        HINDX( 3 ) = 9

!........  Mobile source
      CASE( 'MOBILE' )
        NCOLS = 3
        HINDX( 1 ) = 1
        HINDX( 2 ) = 2
        HINDX( 3 ) = 9

!...........  Determine if link will be included in reports
        IF( ALLOCATED( CLINK ) ) THEN
            DO S = 1, NSRC
                IF( CLINK( S ) .NE. ' ' ) THEN
                    LFLAG = .TRUE.
                    EXIT
                END IF
            END DO
        END IF
        IF( LFLAG ) THEN
            HINDX( 4 ) = 8
            NCOLS = NCOLS + 1
        END IF

!........  Point source
      CASE( 'POINT' )
        NCOLS = 1
        HINDX( 1 ) = 1
        DO J = 1, NCOUT
            HINDX( J+1 ) = J+1
            IF( JSCC .EQ. J ) HINDX( J+1 ) = 9  ! Ensure SCC is labeled as such
            NCOLS = NCOLS + 1
        END DO

!............  If SCC is not included in source chars, then add it to report
        IF( JSCC .EQ. 0 ) HINDX( NCOUT ) = 9

    END SELECT

!.........  Reset max characteristics widths based on headers
    HCWID( 1 ) = 6
    DO J = 1, NCOUT
        HCWID( 1+J ) = LEN_TRIM( HEADERS( HINDX(1+J) ))
        MCWID( J ) = MAX( MCWID(J), HCWID( 1+J ) )
        HCWID( 1+J ) = MCWID( J )
    END DO

!.........  Create remaining header widths - must be consistent with format
!           statments. Maximum is 15 because of HEADERS string size.
    NCOLS = NCOLS + 6
    HCWID( NCOUT + 2 ) = MAX( LEN_TRIM( HEADERS( 10 ) ), 10 )
    HCWID( NCOUT + 3 ) = MAX( LEN_TRIM( HEADERS( 11 ) ), 10 )
    HCWID( NCOUT + 4 ) = MAX( LEN_TRIM( HEADERS( 12 ) ), 7 )
    HCWID( NCOUT + 5 ) = MAX( LEN_TRIM( HEADERS( 13 ) ), SCCLEN3 )
    HCWID( NCOUT + 6 ) = MAX( LEN_TRIM( HEADERS( 14 ) ), SPNLEN3 )
    HCWID( NCOUT + 7 ) = MAX( LEN_TRIM( HEADERS( 15 ) ), 10 )

!.........  Set remaining header indices
    I = 9
    DO J = NCOUT+2, NCOLS
        I = I + 1
        HINDX( J ) = I
    END DO

!.........  Write column headers
    HDRSTR = ' ' // HEADERS( HINDX( 1 ) )( 1:HCWID( 1 ) ) // ';'
    DO J = 2, NCOLS-1
        BUFFER = HDRSTR
        HDRSTR = TRIM( BUFFER ) // ' ' //&
        &         HEADERS( HINDX( J ) )( 1:HCWID( J ) ) // ';'
    END DO
    BUFFER = HDRSTR
    HDRSTR = TRIM( BUFFER ) // ' ' //&
    &         HEADERS( HINDX( NCOLS ) )( 1:HCWID( NCOLS ) )

    WRITE( RDEV, 93000 ) TRIM( HDRSTR )

!.........  Write units headers
    HDRSTR = '       ;'
    DO I = 1, NCOUT
        BUFFER = HDRSTR
        HDRSTR = TRIM( BUFFER ) // REPEAT( ' ',MCWID( I )+1 ) // ';'
    END DO
    BUFFER = HDRSTR
    HDRSTR = TRIM( BUFFER ) // ' [tons/day];'//&
    &         '    [tons/day];            ;           ;' //&
    &         '                ; [frac/year]'

    WRITE( RDEV, 93000 ) TRIM( HDRSTR )
    WRITE( RDEV, 93000 ) REPEAT( '-', LEN_TRIM( HDRSTR ) )

!.........  Develop format for output report. Use buffer to prevent run-time
!           error with Portland Group compilers.
    OUTFMT = '(I7,"; "'
    DO J = 1, NCOUT
        BUFFER = OUTFMT
        WRITE( OUTFMT, '(A,I2.2,A)' ) TRIM( BUFFER ) // ',A',&
        &       MCWID( J ), ',"; "'
    END DO
    BUFFER = OUTFMT
    WRITE( OUTFMT,'(7(A,:,I2.2))' ) TRIM( BUFFER ) // ' E',&
    &       HCWID( NCOUT+2 ), '.4,"; ",E',&
    &       HCWID( NCOUT+3 ), '.4,"; ",F',&
    &       HCWID( NCOUT+4 ), '.4,"; ",A',&
    &       HCWID( NCOUT+5 ), '  ,"; ",A',&
    &       HCWID( NCOUT+6 ), '  ,"; ",F',&
    &       HCWID( NCOUT+7 ), '.4)'

!.........  Read speciation profiles file

    MESG = 'Reading SPECIATION PROFILES file for ' // TRIM( RPOL )
    CALL M3MSG2( MESG )

    CALL RDSPROF( PDEV, RPOL, NMSPC )

!.........  Ensure that profile(s) exist for this pollutant
    IF( NSPFUL .LE. 0 ) THEN
        MESG = 'No speciation profiles found for pollutant "' //&
        &       TRIM( RPOL ) // '"! '
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Abridge profiles so that there is an array of unique profiles
    V = INDEX1( RPOL, NIPOL, EINAM )
    CALL PROCSPRO( NMSPC, SPCNAMES( 1,V ) )

!.........  Allocate memory for the compressed reactivity matrix.  Use data
!           structures for point sources, but this routine can be used for area
!           sources or mobile sources as well.
!.........  Allocate memory for names of output variables

    ALLOCATE(   PCRIDX( NSREAC ),&
    &          PCRREPEM( NSREAC ),&
    &          PCRPRJFC( NSREAC ),&
    &          PCRMKTPN( NSREAC ),&
    &           PCRCSCC( NSREAC ),&
    &          PCRSPROF( NSREAC ),&
    &          MASSONAM( NMSPC ),&
    &          MOLEONAM( NMSPC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'MOLEONAM', PROGNAME )

    IF( MASSOUT ) THEN
        ALLOCATE( RMTXMASS( NSREAC, NMSPC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RMTXMASS', PROGNAME )
        RMTXMASS = 0.
    END IF

    IF( MOLEOUT ) THEN
        ALLOCATE( RMTXMOLE( NSREAC, NMSPC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RMTXMOLE', PROGNAME )
        RMTXMOLE = 0.
    END IF

!.........  Read emissions for current pollutant from the inventory file
!.........  Note that average day read will not work if RPOL is
!           more than 13 characters - must use BLDENAMS routine to
!           do this correctly.
    IF( LAVEDAY ) THEN
        VNAM = AVEDAYRT // RPOL( 1:MIN( LEN_TRIM( RPOL ), 13 ) )
    ELSE
        VNAM = RPOL
    END IF

    CALL RDMAPPOL( NSRC, 1, 1, VNAM, EMIS )

!.........  Loop through all sources and store reactivity information for
!           those that have it
    YFAC = YR2DAY( BYEAR )  ! for loop below
    N = 0
    DO S = 1, NSRC

        K = ISREA( S )       ! index to reactivity data tables

        IF( K .GT. 0 ) THEN

!................. Adjust annual emissions values
            IF( .NOT. LAVEDAY ) EMIS( S ) = EMIS( S ) * YFAC

!.................  For storing inventory emissions, use the value from the
!                   reactivity table only if it is greater than zero. Otherwise,
!                   store the original base-year emissions from the inventory.
            IF( EMREPREA( K ) .GT. 0 ) THEN
                EMREP = EMREPREA( K )
            ELSE
                EMREP = EMIS( S )
            END IF

!.................  Make sure speciation profile is valid for the current source
!                   and pollutant
            SPCODE = ADJUSTR( CSPFREA( K ) )
            V = MAX( FINDC( SPCODE, NSPROF, SPROFN ), 0 )

            IF( V .EQ. 0 ) THEN

                CALL WRITE_WARNING
                CYCLE  ! To next iteration

!.................  Get indices to full speciation table
            ELSE

                ITBL = IDXSPRO ( V )
                NTBL = NSPECIES( V )

            END IF

!.................  Increment counter of reactivity sources
            N = N + 1

!.................  Store various parameters in reactivity matrix output
!                   structures, but double-check for array overflow
            IF( N .LE. NSREAC ) THEN

                PCRIDX  ( N ) = S
                PCRREPEM( N ) = EMREP
                PCRPRJFC( N ) = PRJFCREA( K )
                PCRMKTPN( N ) = MKTPNREA( K )
                PCRCSCC ( N ) = CSCCREA ( K )
                PCRSPROF( N ) = SPCODE

!.....................  Store speciation factors based on speciation profiles
!                       and indices retrieved for SPCODE
                I = ITBL - 1
                DO P = 1, NTBL   ! Loop over species for this profile

                    I = I + 1
                    J = IDXSSPEC( V,P )

                    IF( MASSOUT ) THEN
                        RMTXMASS( N,J )= MASSFACT( I )
                    END IF

                    IF( MOLEOUT ) THEN
                        RMTXMOLE( N,J )= MOLEFACT( I )
                    END IF

                END DO

            END IF  ! End check for overflow

!.................  Format source characteristic information for report
            CALL PARSCSRC( CSOURC( S ), MXCHRS, SC_BEGP, SC_ENDP,&
            &               LF, NC, CHARS )

!..............  Write report entry: source number, source characteristics,
!                current base-year emissions, replacement base-year emissions,
!                projection factor, future-year SCC, future-year profile number,
!                market penetration of new SCC/speciation.
            WRITE( RDEV, OUTFMT ) S,&
            &     ( CHARS( J )( 1:MCWID(J) ), J=1,NCOUT ), EMIS( S ),&
            &       EMREP, PRJFCREA( K ), CSCCREA ( K ), SPCODE,&
            &       MKTPNREA( K )

        END IF      ! End check if reactivity applies to this source

    END DO

!.........  Reset number of reactivity sources in case any were dropped
    NSREAC = N

    IF( NSREAC .LE. 0 ) THEN
        MESG = 'No sources assigned reactivity controls.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Find position of pollutant in list
    V = INDEX1( RPOL, NIPOL, EINAM )

!.........  Set up for and open output reactivity matrices for current pollutant

    CALL OPENRMAT( ENAME, RPOL, MASSOUT, MOLEOUT,&
    &               BYEAR, PYEAR, NSREAC, NMSPC, SPCNAMES( 1,V ),&
    &               SDEV, SNAME, LNAME, MASSONAM, MOLEONAM )

!.........  Write reactivity matrices for current pollutant
    IF( MASSOUT ) THEN
        CALL WRRMAT( NSREAC, NMSPC, SDEV, SNAME, PCRIDX, PCRREPEM,&
        &             PCRPRJFC, PCRMKTPN, RMTXMASS,&
        &             PCRCSCC, PCRSPROF, MASSONAM )

    END IF

    IF( MOLEOUT ) THEN
        CALL WRRMAT( NSREAC, NMSPC, SDEV, LNAME, PCRIDX, PCRREPEM,&
        &             PCRPRJFC, PCRMKTPN, RMTXMOLE,&
        &             PCRCSCC, PCRSPROF, MOLEONAM )

    END IF

!.........  Deallocate memory used to generate reactivity matrices
    DEALLOCATE( INPRF, SPECID, MOLEFACT, MASSFACT )

!        DEALLOCATE( PCRIDX, PCRREPEM, PCRPRJFC, PCRMKTPN, PCRCSCC )
!     &              PCRSPROF )

!        IF( ALLOCATED( RMTXMASS ) ) DEALLOCATE( RMTXMASS )
!        IF( ALLOCATED( RMTXMOLE ) ) DEALLOCATE( RMTXMOLE )

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

93390 FORMAT( A, I4.4 )


!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

!******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

!.............  This internal subprogram writes a warning message for
!               a bad speciation profile for either point sources or
!               other sources.
    SUBROUTINE WRITE_WARNING

!.............  Local variables
        INTEGER                L2

        CHARACTER(300)     BUFFER   ! source information buffer
        CHARACTER(300)     MESG     ! message buffer
        CHARACTER(FIPLEN3) CFIP     ! tmp (character) FIPS code
        CHARACTER(SRCLEN3) CSRC     ! tmp source chars string
        CHARACTER(SCCLEN3) TSCC     ! tmp 10-digit SCC

!----------------------------------------------------------------------

!.............  Retrieve SCC for this source
        TSCC = CSCC( S )   ! SCC

!.............  Generate point source part of message
        IF( PFLAG ) THEN

            CSRC = CSOURC( S )   ! source characteristics
            CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )

            BUFFER = BUFFER( 1:L2 ) // CRLF() // BLANK10 //&
            &         'SCC: ' // TSCC // ' POL: ' // RPOL

!.............  Generate area or mobile source part of message
        ELSE

            CFIP = CIFIP( S )

            BUFFER = 'FIP: ' // CFIP // ' SCC: ' // TSCC //&
            &         ' POL: ' // RPOL&
&

        END IF

        MESG = 'WARNING: Speciation profile "' // SPCODE //&
        &       '" is not in profiles, but it was assigned'//&
        &       CRLF() // BLANK10 // 'to source:' //&
        &       CRLF() // BLANK10 // TRIM( BUFFER )//&
        &       CRLF() // BLANK5 //&
        &       'Source will be excluded from reactivity matrix.'
        CALL M3MESG( MESG )

        RETURN

    END SUBROUTINE WRITE_WARNING

END SUBROUTINE GENREACT

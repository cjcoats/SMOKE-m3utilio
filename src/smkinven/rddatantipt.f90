
SUBROUTINE RDDATAORLPT( LINE, READDATA, READPOL, IYEAR, DESC,    &
                        ERPTYP, SRCTYP, HT, DM, TK, FL, VL, SIC,    &
                        MACT, NAICS, CTYPE, LAT, LON, UTMZ,    &
                        NEID, CORS, BLID, FUGHT, FUGAR,    &
                        EXTORL, HDRFLAG, EFLAG )

    !***********************************************************************
    !  subroutine body starts at line 156
    !
    !  DESCRIPTION:
    !      This subroutine processes a line from an ORL format point-source inventory
    !      file and returns the inventory data values.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created by C. Seppanen (01/03) based on rddataidapt.f
    !       Version June 2016 by Carlie Coats:  add fugitive-emissions properties
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

    !.......   MODULES for public variables
    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: NEM, NDY, NEF, NCE, NRE, NC1, NC2

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*),       INTENT  (IN) :: LINE                      ! input line
    CHARACTER(*),       INTENT (OUT) :: READDATA( 1,NPTPPOL3 )    ! array of data values
    CHARACTER(NAMLEN3), INTENT (OUT) :: READPOL( 1 )              ! array of pollutant names
    INTEGER,            INTENT (OUT) :: IYEAR                     ! inventory year
    CHARACTER(40),      INTENT (OUT) :: DESC                      ! plant description
    CHARACTER(ERPLEN3), INTENT (OUT) :: ERPTYP                    ! emissions release point type
    CHARACTER(STPLEN3), INTENT (OUT) :: SRCTYP                    ! source type code
    CHARACTER(16),      INTENT (OUT) :: HT                        ! stack height
    CHARACTER(16),      INTENT (OUT) :: DM                        ! stack diameter
    CHARACTER(16),      INTENT (OUT) :: TK                        ! exit temperature
    CHARACTER(16),      INTENT (OUT) :: FL                        ! flow rate
    CHARACTER(16),      INTENT (OUT) :: VL                        ! exit velocity
    CHARACTER(SICLEN3), INTENT (OUT) :: SIC                       ! SIC
    CHARACTER(MACLEN3), INTENT (OUT) :: MACT                      ! MACT code
    CHARACTER(NAILEN3), INTENT (OUT) :: NAICS                     ! NAICS code
    CHARACTER,          INTENT (OUT) :: CTYPE                     ! coordinate type
    CHARACTER(16),      INTENT (OUT) :: LAT                       ! stack latitude
    CHARACTER(16),      INTENT (OUT) :: LON                       ! stack longitude
    CHARACTER(2),       INTENT (OUT) :: UTMZ                      ! UTM zone
    CHARACTER(NEILEN3), INTENT (OUT) :: NEID                      ! NEI unique ID
    CHARACTER(ORSLEN3), INTENT (OUT) :: CORS                      ! DOE plant ID
    CHARACTER(BLRLEN3), INTENT (OUT) :: BLID                      ! boiler ID
    CHARACTER(16),      INTENT (OUT) :: FUGHT                     ! RELEASE_HEIGHT_FUGITIVE
    CHARACTER(16),      INTENT (OUT) :: FUGAR                     ! HORIZONTAL_AREA_FUGITIVE
    CHARACTER(EXTLEN3), INTENT (OUT) :: EXTORL                    ! additional ext vars
    LOGICAL,            INTENT (OUT) :: HDRFLAG                   ! true: line is a header line
    LOGICAL,            INTENT (OUT) :: EFLAG                     ! error flag

    !.......   Local parameters, indpendent
    INTEGER      , PARAMETER :: MXPOLFIL = 60      ! arbitrary maximum pollutants in file
    INTEGER      , PARAMETER :: NSEG     = 75          ! number of segments in line
    CHARACTER(16), PARAMETER :: PROGNAME = 'RDDATAORLPT'     ! Program name


    !.......   Other local variables
    INTEGER         I, L, L1, LL           ! counters and indices

    INTEGER, SAVE:: ICC         !  position of CNTRY in CTRYNAM
    INTEGER, SAVE:: INY         !  inventory year
    INTEGER         IOS         !  i/o status
    INTEGER, SAVE:: NPOL        !  number of pollutants in file
    LOGICAL      :: BLKFLAG  = .TRUE.      ! true when it is blank

    CHARACTER(40)      TMPSEG              ! tmp segments of line
    CHARACTER(40)      SEGMENT( NSEG )     ! segments of line
    CHARACTER(CASLEN3) TCAS                ! tmp cas number
    CHARACTER(300)     MESG                ! message buffer

    !***********************************************************************
    !   begin body of subroutine RDDATAORLPT

    !.......  Scan for header lines and check to ensure all are set
    !           properly
    CALL GETHDR( MXPOLFIL, .TRUE., .TRUE., .FALSE.,    &
                 LINE, ICC, INY, NPOL, IOS )

    !.......  Interpret error status
    IF( IOS == 4 ) THEN
        WRITE( MESG,94010 )    &
               'Maximum allowed data variables MXPOLFIL=', MXPOLFIL,    &
               CRLF() // BLANK10 // ' exceeded in input file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    ELSE IF( IOS > 0 ) THEN
        EFLAG = .TRUE.

    END IF

    !.......  If a header line was encountered, set flag and return
    IF( IOS >= 0 ) THEN
        HDRFLAG = .TRUE.
        IYEAR = INY
        RETURN
    ELSE
        HDRFLAG = .FALSE.
    END IF

    !.......  Separate line into segments
    CALL PARSLINE( LINE, NSEG, SEGMENT )

    !.......  Use the file format definition to parse the line into
    !           the various data fields
    DESC   = ADJUSTL( SEGMENT( 6 ) )       ! plant description
    ERPTYP = ADJUSTL( SEGMENT( 8 ) )       ! emissions release point type
    SRCTYP = ADJUSTL( SEGMENT( 9 ) )       ! source type code
    HT     = SEGMENT( 10 )                 ! stack height
    DM     = SEGMENT( 11 )                 ! stack diameter
    TK     = SEGMENT( 12 )                 ! exit temperature
    FL     = SEGMENT( 13 )                 ! flow rate
    VL     = SEGMENT( 14 )                 ! exit velocity
    SIC    = SEGMENT( 15 )                 ! SIC
    MACT   = ADJUSTL( SEGMENT( 16 ) )      ! MACT code
    NAICS  = ADJUSTL( SEGMENT( 17 ) )      ! NAICS code
    CTYPE  = ADJUSTL( SEGMENT( 18 ) )      ! coordinate type
    LON    = SEGMENT( 19 )                 ! stack longitude
    LAT    = SEGMENT( 20 )                 ! stack latitude
    UTMZ   = ADJUSTL( SEGMENT( 21 ) )      ! UTM zone

    READPOL ( 1     ) = SEGMENT( 22 )
    READDATA( 1,NEM ) = SEGMENT( 23 )     ! annual emissions
    READDATA( 1,NDY ) = SEGMENT( 24 )     ! average-day emissions
    READDATA( 1,NEF ) = ' '               ! emission factor
    READDATA( 1,NCE ) = SEGMENT( 25 )     ! control efficiency
    READDATA( 1,NRE ) = SEGMENT( 26 )     ! rule effectiveness
    READDATA( 1,NC1 ) = SEGMENT( 27 )     ! primary control equipment code
    READDATA( 1,NC2 ) = SEGMENT( 28 )     ! secondary control equipment code

    NEID   = ADJUSTL( SEGMENT( 29 ) )      ! NEI Unique ID
    IF( NEID == ' ' ) NEID = '-9'
    CORS   = ADJUSTL( SEGMENT( 30 ) )      ! DOE plant ID
    BLID   = ADJUSTL( SEGMENT( 31 ) )      ! boiler ID
    FUGAR  = SEGMENT( 38 )
    FUGHT  = SEGMENT( 39 )

    !.......  Read extended orl variables and store it as string
    EXTORL = ' '
    DO I = 32, 70
        IF( SEGMENT( I ) == ' ' ) THEN
            TMPSEG = ','
        ELSE
            TMPSEG = ',' // TRIM( SEGMENT( I ) )
            BLKFLAG = .FALSE.
        ENDIF

        EXTORL = TRIM( EXTORL ) // TRIM( TMPSEG )
    END DO

    IF( BLKFLAG ) EXTORL = ' '

    !.......  Return from subroutine
    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Formatted file I/O formats...... 93xxx

93000 FORMAT( A )

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

94120 FORMAT( I6.6 )

94125 FORMAT( I5 )

END SUBROUTINE RDDATAORLPT

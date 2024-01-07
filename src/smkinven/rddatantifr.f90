
SUBROUTINE RDDATAORLFR( LINE, READDATA, READPOL, NDATPERLN, &
                        IYEAR, DESC, SIC, MACT, CTYPE,      &
                        LAT, LON, HDRFLAG, EFLAG )

    !***********************************************************************
    !  DESCRIPTION:
    !      This subroutine processes a line from an ORL FIRE format wildfire
    !      point-source inventory file and returns the unique source characteristics.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !      Created by B.H. Baek (02/06) based on rddatantipt.f
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
    USE MODINFO, ONLY: NEM, NDY, NEF, NCE, NRE, NC1, NC2, TMPNAM

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*),       INTENT  (IN) :: LINE                      ! input line
    CHARACTER(*),       INTENT (OUT) :: READDATA( NDATPERLN,NPTPPOL3 )    ! array of data values
    CHARACTER(NAMLEN3), INTENT (OUT) :: READPOL( NDATPERLN )              ! array of pollutant names
    INTEGER,            INTENT(INOUT):: NDATPERLN                 ! number of data values per line
    INTEGER,            INTENT (OUT) :: IYEAR                     ! inventory year
    CHARACTER(40),      INTENT (OUT) :: DESC                      ! plant description
    CHARACTER(SICLEN3), INTENT (OUT) :: SIC                       ! Material burned code (stored in SIC)
    CHARACTER(MACLEN3), INTENT (OUT) :: MACT                      ! NFDRS code (stored in MACT)
    CHARACTER,          INTENT (OUT) :: CTYPE                     ! coordinate type
    CHARACTER(9),       INTENT (OUT) :: LAT                       ! stack latitude
    CHARACTER(9),       INTENT (OUT) :: LON                       ! stack longitude
    LOGICAL,            INTENT (OUT) :: HDRFLAG                   ! true: line is a header line
    LOGICAL,            INTENT (OUT) :: EFLAG                     ! error flag

    !.......   Local parameters, indpendent
    INTEGER, PARAMETER :: MXPOLFIL = 1000   ! arbitrary maximum pollutants in file
    INTEGER, PARAMETER :: NSEG     =   10   ! number of segments in line for format
    INTEGER, PARAMETER :: NEXTRA   =    4   ! number of extra non-data fields that need
                                            ! to be added as "pollutants)
    CHARACTER(16), PARAMETER :: PROGNAME = 'RDDATAORLFR'     ! Program name
    CHARACTER(NAMLEN3), PARAMETER :: FIREVNAM( NEXTRA ) =     &    ! fire variable names
                        (/ 'HEATCONTENT     ',    &
                           'HFLUX           ',    &
                           'ENDHOUR         ',    &
                           'BEGHOUR         '  /)

    !.......   Other local variables
    INTEGER         I, L           ! counters and indices

    INTEGER, SAVE:: ICC         !  position of CNTRY in CTRYNAM
    INTEGER, SAVE:: INY         !  inventory year
    INTEGER         IOS         !  i/o status
    INTEGER, SAVE:: NDAT = -1     !  number of data values as set by header

    LOGICAL, SAVE:: FIRSTDATA = .TRUE.     ! true: first time data row in encountered

    CHARACTER(40)      SEGMENT( NSEG )     ! segments of line
    CHARACTER(300)     MESG                ! message buffer

    !***********************************************************************
    !   begin body of subroutine RDDATAORLFR

    !.......  Scan for header lines and check to ensure all are set
    !           properly
    CALL GETHDR( MXPOLFIL, .TRUE., .TRUE., .TRUE., LINE, ICC, INY, NDAT, IOS )

    !.......  Interpret error status
    IF( IOS == 4 ) THEN
        WRITE( MESG,94010 )    &
               'Maximum allowed data variables MXPOLFIL=', MXPOLFIL,    &
               CRLF() // BLANK10 //' exceeded in input file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    ELSE IF( IOS > 0 ) THEN
        EFLAG = .TRUE.
    END IF

    !.......  If a header line was encountered, set flag and return
    IF( IOS >= 0 ) THEN
        HDRFLAG = .TRUE.
        IYEAR = INY
        IF( NDAT > 0 ) NDATPERLN = NDAT + NEXTRA
        RETURN
    ELSE
        HDRFLAG = .FALSE.
    END IF

    !.....  Give error if #DATA line has not defined the pollutants that
    !       are contained in the day-specific data file. NOTE: This code
    !       will not be reached until after the last header line)
    IF ( NDAT < 0 ) THEN
        WRITE( MESG, 94010 ) 'First data line reached in ORL '//    &
               'FIRE file with required #DATA header found.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    !.....  Otherwise, the READPOL array has been redfined and now can
    !       be populated with the TMPNAM array set by the GETHDR routine
    ELSE IF ( FIRSTDATA ) THEN

        FIRSTDATA = .FALSE.
        READPOL( 1:NEXTRA ) = FIREVNAM( 1:NEXTRA )
        READPOL( NEXTRA+1:NDATPERLN ) = TMPNAM( 1:NDAT )      ! array

    END IF

    !.......  Separate line into segments
    CALL PARSLINE( LINE, NSEG, SEGMENT )

    !.....  Use the file format definition to parse the line into
    !       the various data fields
    DESC   = SEGMENT( 5 )                  ! fire description
    LAT    = SEGMENT( 6 )                  ! fire latitude
    LON    = SEGMENT( 7 )                  ! fire longitude
    MACT   = SEGMENT( 8 )                  ! MACT code (NFDRSCODE field)
    SIC    = SEGMENT( 9 )                  ! SIC (MATBURNED field)
    CTYPE  = 'L'                           ! lat-lon coordinate type part of format

    !.......  Populate all of the data fields with dummy values
    DO I = 1, NDATPERLN
        READDATA( I,NEM ) = '0.0'       ! dummy annual emissions
        READDATA( I,NDY ) = '-9'        ! dummy average-day emissions
        READDATA( I,NEF ) = '-9'        ! dummy emission factor
        READDATA( I,NCE ) = '-9'        ! dummy control efficiency
        READDATA( I,NRE ) = '-9'        ! dummy rule effectiveness
        READDATA( I,NC1 ) = '-9'        ! dummy primary control equipment code
        READDATA( I,NC2 ) = '-9'        ! dummy secondary control equipment code
    END DO

    !.......  Populate heat content value
    L = LEN( READDATA( 1, NEM ) )
    READDATA( 1, NEM ) = SEGMENT( 10 )( 1:L )

    RETURN

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDDATAORLFR

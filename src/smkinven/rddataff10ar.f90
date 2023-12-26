
SUBROUTINE RDDATAFF10AR( LINE, READDATA, READPOL, IYEAR,    &
                         SRCTYP, SIC, SHAPE, EXTORL,        &
                         HDRFLAG, AVEFLAG, EFLAG )

    !***********************************************************************
    !  subroutine body starts at line 156
    !
    !  DESCRIPTION:
    !      This subroutine processes a line from an ORL format area-source inventory
    !      file and returns the inventory data values.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created by B.H. Baek (Aug, 2011)
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !**************************************************************************
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
    !***************************************************************************

    USE M3UTILIO

    !......   MODULES for public variables
    !......  This module contains the information about the source category
    USE MODINFO, ONLY: NEM, NDY, NEF, NCE, NRE, NRP, INV_MON

    !......  This module contains data for day- and hour-specific data
    USE MODDAYHR, ONLY: FF10INVFLAG

    IMPLICIT NONE

    !......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !......   SUBROUTINE ARGUMENTS
    CHARACTER(*),       INTENT  (IN) :: LINE                      ! input line
    CHARACTER(*),       INTENT (OUT) :: READDATA( 1,NARPPOL3 )    ! array of data values
    CHARACTER(IOVLEN3), INTENT (OUT) :: READPOL( 1 )              ! pollutant name
    INTEGER,            INTENT (OUT) :: IYEAR                     ! inventory year
    CHARACTER(STPLEN3), INTENT (OUT) :: SRCTYP                    ! source type code
    CHARACTER(SICLEN3), INTENT (OUT) :: SIC                       ! SIC
    CHARACTER(SHPLEN3), INTENT (OUT) :: SHAPE                     ! SHAPE_ID
    CHARACTER(EXTLEN3), INTENT (OUT) :: EXTORL                    ! additional ext vars
    LOGICAL,            INTENT (OUT) :: HDRFLAG                   ! true: line is a header line
    LOGICAL,            INTENT (OUT) :: AVEFLAG                   ! true: Aveday inv is processed
    LOGICAL,            INTENT (OUT) :: EFLAG                     ! error flag

    !......   EXTERNAL FUNCTIONS and their descriptions:
    LOGICAL, EXTERNAL :: CHKINT

    !......   Local parameters
    INTEGER      , PARAMETER :: MXDATFIL = 60               ! arbitrary max no. data variables
    INTEGER      , PARAMETER :: NSEG     = 60               ! number of segments in line
    CHARACTER(16), PARAMETER :: PROGNAME = 'RDDATAFF10AR'   ! Program name

    !......   Other local variables
    INTEGER         I            ! counters and indices
    INTEGER, SAVE:: ICC          !  position of CNTRY in CTRYNAM
    INTEGER, SAVE:: INY          !  inventory year
    INTEGER         IOS          !  i/o status
    INTEGER, SAVE:: NPOA         !  number of pollutants in file
    INTEGER      :: MDAYS        !  days of modeling inventory month
    INTEGER      :: LYEAR        !  Leap year (366 days per year)
    REAL         :: AVEINV       !  annual total estimate from monthly total VMT

    CHARACTER(25)      SEGMENT( NSEG )     ! segments of line
    CHARACTER(25)      TMPSEG              ! tmp segments of line
    CHARACTER(CASLEN3) TCAS                ! tmp cas number
    CHARACTER(300)     MESG                !  message buffer

    !***********************************************************************
    !   begin body of subroutine RDDATAFF10AR

    !......  Scan for header lines and check to ensure all are set
    !        properly (country and year required)
    CALL GETHDR( MXDATFIL, .TRUE., .TRUE., .FALSE., LINE, ICC, INY, NPOA, IOS )

    !......  Interpret error status
    IF( IOS == 4 ) THEN
        WRITE( MESG,94010 )   &
               'Maximum allowed data variables MXDATFIL=', MXDATFIL,    &
               CRLF() // BLANK10 // ' exceeded in input file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    ELSE IF( IOS > 0 ) THEN
        EFLAG = .TRUE.

    END IF

    !......  If a header line was encountered, set flag and return
    IF( IOS >= 0 ) THEN
        HDRFLAG = .TRUE.
        IYEAR = INY
        RETURN
    ELSE
        HDRFLAG = .FALSE.
    END IF

    !......  Separate line into segments
    CALL PARSLINE( LINE, NSEG, SEGMENT )

    !...... Return if the first line is a header line
    IF( .NOT. CHKINT( SEGMENT( 2 ) ) ) THEN
        HDRFLAG = .TRUE.
        RETURN
    END IF

    !......  Use the file format definition to parse the line into
    !        the various data fields
    READPOL ( 1     ) = ADJUSTL( SEGMENT( 8 ) )
    READDATA( 1,NEM ) = SEGMENT( 9 )
    READDATA( 1,NDY ) = ''
    READDATA( 1,NEF ) = '-9'
    READDATA( 1,NCE ) = '-9'
    READDATA( 1,NRE ) = '-9'
    READDATA( 1,NRP ) = '-9'

    SRCTYP = ' '       ! source type code
    SIC    = ' '
    SHAPE  = ADJUSTL( SEGMENT(  5 ) )
    EXTORL = ' '       ! extended orl (N/A)

    !......  Compute annual total based on monthly total
    AVEINV = 0.0
    AVEFLAG = .FALSE.
    IF( INV_MON > 0 ) THEN

        DO I = 1, 12
            IF( LEN_TRIM( SEGMENT( 20+I ) ) < 1 ) THEN
                SEGMENT( 20+I ) = '0.0'
            END IF
            AVEINV = AVEINV + STR2REAL( SEGMENT( 20+I ) )
        END DO

        IF( AVEINV > 0.0 ) AVEFLAG = .TRUE.

        IF( AVEFLAG ) THEN

            READDATA( 1,NEM ) = '0.0'
            READDATA( 1,NDY ) = SEGMENT( 20 + INV_MON )

            MDAYS = MON_DAYS( INV_MON )        ! day of months

            LYEAR = INT( 1 / YR2DAY ( INY ) )                    ! convert year to days
            IF( LYEAR > 365 .AND. INV_MON == 2 ) MDAYS = 29      ! leap year (feb = 29days)

            AVEINV = STR2REAL( READDATA(1,NDY) ) / MDAYS
            WRITE( READDATA( 1,NDY ), '( E15.10 )' ) AVEINV

            IF( AVEINV < 0.0 ) THEN
                MESG = 'ERROR: Can not process negative value'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        END IF

    END IF

    !......  Reset annual total inventory to zero for daily/hourly FF10 processing
    IF( FF10INVFLAG ) THEN
        READPOL ( 1     ) = SEGMENT( 9 )
        READDATA( 1,NEM ) = ''
        READDATA( 1,NDY ) = ''
    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************
    !......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDDATAFF10AR

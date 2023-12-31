
SUBROUTINE RDSRCFF10PT( LINE, CFIP, FCID, PTID, SKID, SGID, TSCC,   &
                        CORS, BLID, NPOLPERLN, HDRFLAG, EFLAG )

    !***********************************************************************
    !  subroutine body starts at line 156
    !
    !  DESCRIPTION:
    !      This subroutine processes a line from an FF10 format point-source inventory
    !      file and returns the unique source characteristics.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 10/2011 by Dongmei Yang based on rdsrcntipt.f
    !       Version 10/2016 by C. Coats:  USE M3UTILIO
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
    !......  This module contains the lists of unique inventory information
    USE MODLISTS, ONLY: UCASNKEP, NUNIQCAS, UNIQCAS

    !......  This module contains data for day- and hour-specific data
    USE MODDAYHR, ONLY: DAYINVFLAG, HRLINVFLAG, FF10INVFLAG

    IMPLICIT NONE

    !......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !......   EXTERNAL FUNCTIONS and their descriptions:

    LOGICAL, EXTERNAL :: CHKINT, USEEXPGEO

    !......   SUBROUTINE ARGUMENTS
    CHARACTER(*),       INTENT (IN) :: LINE          ! input line
    CHARACTER(FIPLEN3), INTENT(OUT) :: CFIP          ! fip code
    CHARACTER(PLTLEN3), INTENT(OUT) :: FCID          ! facility/plant ID
    CHARACTER(CHRLEN3), INTENT(OUT) :: PTID          ! point ID
    CHARACTER(CHRLEN3), INTENT(OUT) :: SKID          ! stack ID
    CHARACTER(CHRLEN3), INTENT(OUT) :: SGID          ! segment ID
    CHARACTER(SCCLEN3), INTENT(OUT) :: TSCC          ! scc code
    CHARACTER(ORSLEN3), INTENT(OUT) :: CORS          ! ORIS code
    CHARACTER(BLRLEN3), INTENT(OUT) :: BLID          ! boiler ID
    INTEGER,            INTENT(OUT) :: NPOLPERLN     ! no. pollutants per line
    LOGICAL,            INTENT(OUT) :: HDRFLAG       ! true: line is a header line
    LOGICAL,            INTENT(OUT) :: EFLAG         ! error flag

    !......   Local parameters, indpendent
    INTEGER, PARAMETER :: MXPOLFIL = 60      ! arbitrary maximum pollutants in file
    INTEGER, PARAMETER :: NSEG = 80          ! number of segments in line

    !......   Other local variables
    INTEGER         I, L1, L2           ! counters and indices

    INTEGER, SAVE:: ICC         !  position of CNTRY in CTRYNAM
    INTEGER, SAVE:: INY         !  inventory year
    INTEGER         IOS         !  i/o status
    INTEGER, SAVE:: NPOL        !  number of pollutants in file

    CHARACTER(50)      SEGMENT( NSEG )     ! segments of line
    CHARACTER(CASLEN3) TCAS                ! tmp cas number
    CHARACTER(300)     MESG                ! message buffer
    CHARACTER(500)     LBUF                ! line buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'RDSRFF10LPT'     ! Program name

    !***********************************************************************
    !   begin body of subroutine RDSRCORLPT

    !......  Scan for header lines and check to ensure all are set
    !        properly
    CALL GETHDR( MXPOLFIL, .NOT. USEEXPGEO(), .TRUE., .FALSE.,   &
                 LINE, ICC, INY, NPOL, IOS )

    !......  Interpret error status
    IF( IOS == 4 ) THEN
        WRITE( MESG,94010 )   &
               'Maximum allowed data variables MXPOLFIL=', MXPOLFIL,    &
               CRLF() // BLANK10 // ' exceeded in input file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    ELSE IF( IOS > 0 ) THEN
        EFLAG = .TRUE.

    END IF

    !......  If a header line was encountered, set flag and return
    IF( IOS >= 0 ) THEN
    !......  Determine whether processing daily/hourly inventories or not
        LBUF = LINE
        CALL UPCASE( LBUF )

        IF( FF10INVFLAG ) THEN
            L1 = INDEX( LBUF, 'FF10_DAILY_POINT' )
            L2 = INDEX( LBUF, 'FF10_HOURLY_POINT' )

            IF( .NOT. DAYINVFLAG .AND. L1  > 0 ) THEN
                MESG = 'ERROR: MUST set DAY_SPECIFIC_YN to Y '//   &
                     'to process daily FF10_DAILY_POINT inventory'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF( .NOT. HRLINVFLAG .AND. L2 > 0 ) THEN
                MESG = 'ERROR: MUST set HOUR_SPECIFIC_YN to Y '//   &
                     'to process hourly FF10_HOURLY_POINT inventory'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END IF

        HDRFLAG = .TRUE.
        RETURN
    ELSE
        HDRFLAG = .FALSE.
    END IF

    !......  Separate line into segments
    CALL PARSLINE( LINE, NSEG, SEGMENT )

    !......  Use the file format definition to parse the line into
    !        the various data fields
    CFIP = REPEAT( '0', FIPLEN3 )
    IF( USEEXPGEO() ) THEN
        CFIP(  1: 3 ) = ADJUSTR( SEGMENT( 1 )( 1:3 ) )
        CFIP(  4: 9 ) = ADJUSTR( SEGMENT( 2 )( 1:6 ) )
        CFIP( 10:12 ) = ADJUSTR( SEGMENT( 3 )( 1:3 ) )
    ELSE
        WRITE( CFIP( FIPEXPLEN3+1:FIPEXPLEN3+1 ), '(I1)' ) ICC      ! country code of FIPS
        CFIP( FIPEXPLEN3+2:FIPEXPLEN3+6 ) = ADJUSTR( SEGMENT( 2 )( 1:5 ) )      ! state/county code
    END IF

    !......  Replace blanks with zeros
    DO I = 1,FIPLEN3
        IF( CFIP( I:I ) == ' ' ) CFIP( I:I ) = '0'
    END DO

    FCID = ADJUSTL( SEGMENT( 4 ) )     ! facility/plant ID
    PTID = ADJUSTL( SEGMENT( 5 ) )     ! point ID
    SKID = ADJUSTL( SEGMENT( 6 ) )     ! stack ID
    SGID = ADJUSTL( SEGMENT( 7 ) )     ! segment ID
    CORS = ADJUSTL( SEGMENT( 42 ) )    ! DOE plant ID
    BLID = ADJUSTL( SEGMENT( 43 ) )    ! boiler ID

    !......  Determine number of pollutants for this line based on CAS number
    IF( FF10INVFLAG ) THEN
        TSCC = ADJUSTL( SEGMENT( 8 ) )                               ! SCC code
        TCAS = ADJUSTL( SEGMENT( 9 ) )
    ELSE
        TSCC = ADJUSTL( SEGMENT( 12 ) )                              ! SCC code
        TCAS = ADJUSTL( SEGMENT( 13 ) )
    END IF

    I = FINDC( TCAS, NUNIQCAS, UNIQCAS )
    IF( I < 1 ) THEN
        NPOLPERLN = 0
    ELSE
        NPOLPERLN = UCASNKEP( I )
    END IF

    !......  Return from subroutine
    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !......   Formatted file I/O formats...... 93xxx

93000 FORMAT( A )

    !......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

94120 FORMAT( I6.6 )

94125 FORMAT( I5 )

END SUBROUTINE RDSRCFF10PT

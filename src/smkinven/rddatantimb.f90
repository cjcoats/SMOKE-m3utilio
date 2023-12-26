
SUBROUTINE RDDATAORLMB( LINE, READDATA, READPOL, IYEAR,    &
                        SRCTYP, TSCC, EXTORL, HDRFLAG, EFLAG )

    !***********************************************************************
    !  subroutine body starts at line 156
    !
    !  DESCRIPTION:
    !      This subroutine processes a line from an ORL format mobile-source inventory
    !      file and returns the inventory data values.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !      Created by C. Seppanen (01/03) based on rdntiar.f
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
    USE MODINFO, ONLY: NEM, NDY

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*),       INTENT  (IN) :: LINE                      ! input line
    CHARACTER(*),       INTENT (OUT) :: READDATA( 1,NMBPPOL3 )    ! array of data values
    CHARACTER(NAMLEN3), INTENT (OUT) :: READPOL( 1 )              ! pollutant name
    INTEGER,            INTENT (OUT) :: IYEAR                     ! inventory year
    CHARACTER(STPLEN3), INTENT (OUT) :: SRCTYP                    ! source type code
    CHARACTER(SCCLEN3), INTENT (OUT) :: TSCC                      ! scc code
    CHARACTER(EXTLEN3), INTENT (OUT) :: EXTORL                    ! additional ext vars
    LOGICAL,            INTENT (OUT) :: HDRFLAG                   ! true: line is a header line
    LOGICAL,            INTENT (OUT) :: EFLAG                     ! error flag

    !.......   Local parameters
    INTEGER      , PARAMETER :: MXDATFIL = 60      ! arbitrary max no. data variables
    INTEGER      , PARAMETER :: NSEG     = 30          ! number of segments in line
    CHARACTER(16), PARAMETER :: PROGNAME = 'RDDATAORLMB'     ! Program name

    !.......   Other local variables
    INTEGER         I            ! counters and indices

    INTEGER, SAVE:: ICC          !  position of CNTRY in CTRYNAM
    INTEGER, SAVE:: INY          !  inventory year
    INTEGER         IOS          !  i/o status
    INTEGER, SAVE:: NPOA         !  number of pollutants in file

    LOGICAL      :: BLKFLAG  = .TRUE.      ! true when it is blank

    CHARACTER(25)      SEGMENT( NSEG )     ! segments of line
    CHARACTER(25)      TMPSEG              ! tmp segments of line
    CHARACTER(CASLEN3) TCAS                ! tmp cas number
    CHARACTER(300)     MESG                !  message buffer

    !***********************************************************************
    !   begin body of subroutine RDDATAORLMB

    !.......  Scan for header lines and check to ensure all are set
    !           properly (country and year required)
    CALL GETHDR( MXDATFIL, .TRUE., .TRUE., .FALSE., LINE, ICC, INY, NPOA, IOS )

    !.......  Interpret error status
    IF( IOS == 4 ) THEN
        WRITE( MESG,94010 )    &
               'Maximum allowed data variables MXDATFIL=', MXDATFIL,    &
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
    READPOL ( 1     ) = SEGMENT( 3 )
    READDATA( 1,NEM ) = SEGMENT( 4 )
    READDATA( 1,NDY ) = SEGMENT( 5 )

    TSCC   = SEGMENT( 2 )
    SRCTYP = ADJUSTL( SEGMENT( 6 ) )       ! source type code

    !.......  Read extended orl variables and store it as string
    EXTORL = ' '
    DO I = 7, 16
        IF( SEGMENT( I ) == ' ' ) THEN
            TMPSEG = ','
        ELSE
            TMPSEG = ',' // TRIM( SEGMENT( I ) )
            BLKFLAG = .FALSE.
        ENDIF

        EXTORL = TRIM( EXTORL ) // TRIM( TMPSEG )
    END DO

    IF( BLKFLAG ) EXTORL = ' '

    RETURN

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDDATAORLMB

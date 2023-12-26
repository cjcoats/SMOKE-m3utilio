
SUBROUTINE RDSRCORLAR( LINE, CFIP, TSCC, NPOLPERLN, HDRFLAG, EFLAG )

    !***********************************************************************
    !  subroutine body starts at line 156
    !
    !  DESCRIPTION:
    !      This subroutine processes a line from an ORL format area-source inventory
    !      file and returns the unique source characteristics.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created by C. Seppanen (01/03) based on rdntiar.f
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
    !.......  This module contains the lists of unique inventory information
    USE MODLISTS, ONLY: UCASNKEP, NUNIQCAS, UNIQCAS

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*),       INTENT (IN) :: LINE          ! input line
    CHARACTER(FIPLEN3), INTENT(OUT) :: CFIP          ! fip code
    CHARACTER(SCCLEN3), INTENT(OUT) :: TSCC          ! scc code
    INTEGER,            INTENT(OUT) :: NPOLPERLN     ! no. pollutants per line
    LOGICAL,            INTENT(OUT) :: HDRFLAG       ! true: line is a header line
    LOGICAL,            INTENT(OUT) :: EFLAG         ! error flag

    !.......   Local parameters
    INTEGER      , PARAMETER :: MXDATFIL = 60        ! arbitrary max no. data variables
    INTEGER      , PARAMETER :: NSEG     = 50        ! number of segments in line
    CHARACTER(16), PARAMETER :: PROGNAME = 'RDSRCORLAR'     ! Program name

    !.......   Other local variables
    INTEGER         I            ! counters and indices

    INTEGER, SAVE:: ICC          !  position of CNTRY in CTRYNAM
    INTEGER, SAVE:: INY          !  inventory year
    INTEGER         IOS          !  i/o status
    INTEGER, SAVE:: NPOA         !  number of pollutants in file

    CHARACTER(128)      SEGMENT( NSEG )     ! segments of line
    CHARACTER(CASLEN3) TCAS                ! tmp cas number
    CHARACTER(300)     MESG                !  message buffer

    !***********************************************************************
    !   begin body of subroutine RDSRCORLAR

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
        RETURN
    ELSE
        HDRFLAG = .FALSE.
    END IF

    !.......  Separate line into segments
    CALL PARSLINE( LINE, NSEG, SEGMENT )

    !.......  Use the file format definition to parse the line into
    !           the various data fields
    CFIP = REPEAT( '0', FIPLEN3 )
    WRITE( CFIP( FIPEXPLEN3+1:FIPEXPLEN3+1 ), '(I1)' ) ICC      ! country code of FIPS
    CFIP( FIPEXPLEN3+2:FIPEXPLEN3+6 ) = ADJUSTR( SEGMENT( 1 )( 1:5 ) )      ! state/county code

    !.......  Replace blanks with zeros
    DO I = 1,FIPLEN3
        IF( CFIP( I:I ) == ' ' ) CFIP( I:I ) = '0'
    END DO

    TSCC = SEGMENT( 2 )                  ! SCC code

    !.......  Determine number of pollutants for this line based on CAS number
    TCAS = ADJUSTL( SEGMENT( 3 ) )
    I = FINDC( TCAS, NUNIQCAS, UNIQCAS )
    IF( I < 1 ) THEN
        NPOLPERLN = 0
    ELSE
        NPOLPERLN = UCASNKEP( I )
    END IF

    RETURN

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDSRCORLAR

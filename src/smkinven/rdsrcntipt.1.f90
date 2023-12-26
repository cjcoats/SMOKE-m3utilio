
SUBROUTINE RDSRCORLPT( LINE, CFIP, FCID, PTID, SKID, SGID, TSCC,&
&                       CORS, BLID, NPOLPERLN, HDRFLAG, EFLAG )

!***********************************************************************
!  subroutine body starts at line 156
!
!  DESCRIPTION:
!      This subroutine processes a line from an ORL format point-source inventory
!      file and returns the unique source characteristics.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!      Created by C. Seppanen (01/03) based on rdidapt.f
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
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

!...........   MODULES for public variables
!.........  This module contains the lists of unique inventory information
    USE MODLISTS, ONLY: UCASNKEP, NUNIQCAS, UNIQCAS

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   SUBROUTINE ARGUMENTS
    CHARACTER(*),       INTENT (IN) :: LINE      ! input line
    CHARACTER(FIPLEN3), INTENT(OUT) :: CFIP      ! fip code
    CHARACTER(PLTLEN3), INTENT(OUT) :: FCID      ! facility/plant ID
    CHARACTER(CHRLEN3), INTENT(OUT) :: PTID      ! point ID
    CHARACTER(CHRLEN3), INTENT(OUT) :: SKID      ! stack ID
    CHARACTER(CHRLEN3), INTENT(OUT) :: SGID      ! segment ID
    CHARACTER(SCCLEN3), INTENT(OUT) :: TSCC      ! scc code
    CHARACTER(ORSLEN3), INTENT(OUT) :: CORS      ! DOE plant ID
    CHARACTER(BLRLEN3), INTENT(OUT) :: BLID      ! boiler ID
    INTEGER,            INTENT(OUT) :: NPOLPERLN ! no. pollutants per line
    LOGICAL,            INTENT(OUT) :: HDRFLAG   ! true: line is a header line
    LOGICAL,            INTENT(OUT) :: EFLAG     ! error flag

!...........   Local parameters, indpendent
    INTEGER, PARAMETER :: MXPOLFIL = 60  ! arbitrary maximum pollutants in file
    INTEGER, PARAMETER :: NSEG = 75      ! number of segments in line

!...........   Other local variables
    INTEGER         I       ! counters and indices

    INTEGER, SAVE:: ICC     !  position of CNTRY in CTRYNAM
    INTEGER, SAVE:: INY     !  inventory year
    INTEGER         IOS     !  i/o status
    INTEGER, SAVE:: NPOL    !  number of pollutants in file

    LOGICAL, SAVE:: FIRSTIME = .TRUE. ! true: first time routine is called

    CHARACTER(128) SEGMENT( NSEG ) ! segments of line
    CHARACTER(CASLEN3) TCAS            ! tmp cas number
    CHARACTER(300)     MESG            ! message buffer

    CHARACTER(16) :: PROGNAME = 'RDSRCORLPT' ! Program name

!***********************************************************************
!   begin body of subroutine RDSRCORLPT

!.........  Scan for header lines and check to ensure all are set
!           properly
    CALL GETHDR( MXPOLFIL, .TRUE., .TRUE., .FALSE.,&
    &             LINE, ICC, INY, NPOL, IOS )

!.........  Interpret error status
    IF( IOS == 4 ) THEN
        WRITE( MESG,94010 )&
        &       'Maximum allowed data variables ' //&
        &       '(MXPOLFIL=', MXPOLFIL, CRLF() // BLANK10 //&
        &       ') exceeded in input file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    ELSE IF( IOS > 0 ) THEN
        EFLAG = .TRUE.

    END IF

!.........  If a header line was encountered, set flag and return
    IF( IOS >= 0 ) THEN
        HDRFLAG = .TRUE.
        RETURN
    ELSE
        HDRFLAG = .FALSE.
    END IF

!.........  Separate line into segments
    CALL PARSLINE( LINE, NSEG, SEGMENT )

!.........  Use the file format definition to parse the line into
!           the various data fields
    CFIP = REPEAT( '0', FIPLEN3 )
    WRITE( CFIP( FIPEXPLEN3+1:FIPEXPLEN3+1 ), '(I1)' ) ICC  ! country code of FIPS
    CFIP( FIPEXPLEN3+2:FIPEXPLEN3+6 ) = ADJUSTR( SEGMENT( 1 )( 1:5 ) )  ! state/county code

!.........  Replace blanks with zeros
    DO I = 1,FIPLEN3
        IF( CFIP( I:I ) == ' ' ) CFIP( I:I ) = '0'
    END DO

    FCID = ADJUSTL( SEGMENT( 2 ) ) ! facility/plant ID
    PTID = ADJUSTL( SEGMENT( 3 ) ) ! point ID
    SKID = ADJUSTL( SEGMENT( 4 ) ) ! stack ID
    SGID = ADJUSTL( SEGMENT( 5 ) ) ! segment ID
    TSCC = ADJUSTL( SEGMENT( 7 ) ) ! scc code
    CORS = ADJUSTL( SEGMENT( 30 ) )  ! DOE plant ID
    BLID = ADJUSTL( SEGMENT( 31 ) )  ! boiler ID

!.........  Determine number of pollutants for this line based on CAS number
    TCAS = ADJUSTL( SEGMENT( 22 ) )
    I = FINDC( TCAS, NUNIQCAS, UNIQCAS )
    IF( I < 1 ) THEN
        NPOLPERLN = 0
    ELSE
        NPOLPERLN = UCASNKEP( I )
    END IF

!.........  Make sure routine knows it's been called already
    FIRSTIME = .FALSE.

!.........  Return from subroutine
    RETURN

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDSRCORLPT

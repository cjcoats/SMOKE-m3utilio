
SUBROUTINE RDSRCMEDSPT( LINE, CFIP, FCID, PTID, SKID, SGID, TSCC,    &
                       NPOLPERLN, HDRFLAG, EFLAG )

    !***********************************************************************
    !  subroutine body starts at line 156
    !
    !  DESCRIPTION:
    !      This subroutine processes a line from an MEDS format point-source inventory
    !      file and returns the unique source characteristics.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created on 12/2013 by B.H. Baek
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
    USE MODSOURC, ONLY: NMEDGAI, COABDST

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*),       INTENT (IN) :: LINE          ! input line
    CHARACTER(FIPLEN3), INTENT(OUT) :: CFIP          ! fip code
    CHARACTER(PLTLEN3), INTENT(OUT) :: FCID          ! facility ID
    CHARACTER(CHRLEN3), INTENT(OUT) :: PTID          ! point ID
    CHARACTER(CHRLEN3), INTENT(OUT) :: SKID          ! stack ID
    CHARACTER(CHRLEN3), INTENT(OUT) :: SGID          ! segment ID
    CHARACTER(SCCLEN3), INTENT(OUT) :: TSCC          ! scc code
    INTEGER,            INTENT(OUT) :: NPOLPERLN     ! no. pollutants per line
    LOGICAL,            INTENT(OUT) :: HDRFLAG       ! true: line is a header line
    LOGICAL,            INTENT(OUT) :: EFLAG         ! error flag

    !.......   Local parameters, indpendent
    INTEGER      , PARAMETER :: MXPOLFIL = 53      ! maximum pollutants in file
    CHARACTER(16), PARAMETER :: PROGNAME = 'RDSRCMEDSPT'     ! Program name

    !.......   Other local variables
    INTEGER         I           ! counters and indices

    INTEGER, SAVE:: ICC         !  position of CNTRY in CTRYNAM
    INTEGER, SAVE:: INY         !  inventory year
    INTEGER         IOS         !  i/o status
    INTEGER, SAVE:: NPOL        !  number of pollutants in file
    INTEGER         ROW, COL      ! tmp row and col

    CHARACTER( 3 )    :: STA = '006'       ! State code for CA (006)
    CHARACTER( 3 )       ARBN, CNTY
    CHARACTER(300)       MESG        !  message buffer
    CHARACTER(CHRLEN3)   GAI         !  GAI lookup code

    !***********************************************************************
    !   begin body of subroutine RDSRCMEDSPT

    !.......  Fixed no of pollutants in MEDS
    HDRFLAG = .FALSE.
    NPOLPERLN = 6           ! fixed no of poll (CO,NOx,SOx,TOG,PM,NH3) in MEDS

    !.......  Use the file format definition to parse the line into
    !           the various data fields
    GAI = ADJUSTL( LINE( 71:73 ) )      ! GAI lookup code

    IF( .NOT. ALLOCATED( COABDST ) ) THEN
        EFLAG = .TRUE.
        MESG='ERROR: MUST set IMPORT_MEDS_YN to Y to process'    &
           //' pregridded MEDS-formatted inventory'
        CALL M3MESG( MESG )
    END IF

    I = INDEX1( GAI, NMEDGAI, COABDST( :,1 ) )
    IF( I < 1 ) THEN
        ARBN = ADJUSTR( LINE( 68:70 ) )
        CALL PADZERO( ARBN )
        WRITE( CNTY, '(I3.3)' ) STR2INT( LINE( 57:58 ) )
        CFIP = ARBN // STA // CNTY // '000'
    ELSE
        CFIP = TRIM( COABDST( I,2 ) )        ! FIPS code
    ENDIF

    FCID = ADJUSTL( LINE( 43:51 ) )      ! platn/facility ID
    SKID = ADJUSTL( LINE( 52:56 ) )      ! stack ID
    PTID = ADJUSTL( LINE( 37:39 ) )      ! Column ID
    SGID = ADJUSTL( LINE( 40:42 ) )      ! Row ID
    TSCC = ADJUSTR( LINE(  9:22 ) )      ! EIC code
    CALL PADZERO( TSCC )

    RETURN

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDSRCMEDSPT

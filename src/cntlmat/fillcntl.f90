
SUBROUTINE FILLCNTL( PKTTYP, JTMAX, JXMAX, PKTINFO, JPOL, JT, JX )

    !***********************************************************************
    !  subroutine body starts at line
    !      This routine increments JX and JT and stores the control data tables
    !      and ungrouped control x-ref tables
    !
    !  DESCRIPTION:
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 3/99 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !****************************************************************************/
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

    !.......  MODULES for public variables
    !.......  This module is for cross reference tables
    USE MODXREF, ONLY: INDXTA, ISPTA, CSCCTA, CSRCTA, MPRNA, CMACTA, CISICA

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY

    IMPLICIT NONE

    !.......   INCLUDES:
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters
    INCLUDE 'CPKTDAT.h90'       !  control packet contents


    !.......   SUBROUTINE ARGUMENTS:
    INTEGER     , INTENT (IN) :: PKTTYP     ! packet type number
    INTEGER     , INTENT (IN) :: JTMAX      ! max allowed JT
    INTEGER     , INTENT (IN) :: JXMAX      ! max allowed JX
    TYPE( CPACKET ),INTENT(IN):: PKTINFO    ! packet information
    INTEGER     , INTENT (IN) :: JPOL       ! pollutant position
    INTEGER  , INTENT(IN OUT) :: JT         ! idx to control data tables
    INTEGER  , INTENT(IN OUT) :: JX         ! idx to ungrouped cntl x-ref tables

    !.......   Other local variables
    INTEGER         N                   !  counters and indices
    INTEGER         SIC                 !  tmp SIC

    CHARACTER(SCCLEN3) TMPSCC       !  tmp SCC
    CHARACTER(ALLLEN3) CSRCALL      !  buffer for source char, incl pol

    CHARACTER(16), PARAMETER :: PROGNAME = 'FILLCNTL'       ! program name

    !***********************************************************************
    !   begin body of subroutine FILLCNTL

    !.......  Increment data table counter
    JT = JT + 1

    !.......  Store packet information in temporary variables
    TMPSCC = PKTINFO%TSCC
    IF( PKTINFO%CSIC /= REPEAT( '0', SICLEN3 ) ) THEN
        TMPSCC = PKTINFO%CSIC
    END IF

    JX = JX + 1

    IF( JX .GT. JXMAX ) RETURN      ! to next iteration

    !.......  Store unsorted x-ref table entries
    INDXTA( JX ) = JX
    ISPTA ( JX ) = JPOL

    !.......  Parse the line of data into segments based on the rules
    !.......  Ensure that pollutant is in master list of pollutants or
    !               skip the pollutant-specific entry
    CSCCTA( JX ) = PKTINFO%TSCC
    CMACTA( JX ) = PKTINFO%CMCT
    CISICA( JX ) = PKTINFO%CSIC
    MPRNA ( JX ) = JT       ! Position in data table

    !.......  Store sorting criteria as right-justified in fields
    CSRCALL = ' '
    SELECT CASE( CATEGORY )

      CASE( 'AREA' )
        CALL BLDCSRC( PKTINFO%CFIP, PLTBLNK3, CHRBLNK3, &
                      CHRBLNK3, CHRBLNK3, CHRBLNK3,     &
                      CHRBLNK3, POLBLNK3, CSRCALL )

        CSRCTA( JX ) = CSRCALL( 1:SRCLEN3 ) // TMPSCC //&
                       PKTINFO%CMCT // PKTINFO%CPOS

      CASE( 'MOBILE' )
        CALL BLDCSRC( PKTINFO%CFIP, PLTBLNK3, CHRBLNK3, &
                      CHRBLNK3, CHRBLNK3, CHRBLNK3,     &
                      CHRBLNK3, POLBLNK3, CSRCALL )

        CSRCTA( JX ) = CSRCALL( 1:SRCLEN3 ) // TMPSCC //&
                       PKTINFO%CPOS

      CASE( 'POINT' )

        CALL BLDCSRC( PKTINFO%CFIP, PKTINFO%PLT, PKTINFO%CHAR1, &
                      PKTINFO%CHAR2, PKTINFO%CHAR3,             &
                      PKTINFO%CHAR4, PKTINFO%CHAR5, POLBLNK3,   &
                      CSRCALL )

        CSRCTA( JX ) = CSRCALL( 1:SRCLEN3 ) // TMPSCC // PKTINFO%CMCT // PKTINFO%CPOS

    END SELECT

    !.......  Now store the packet information in the packet tables...

    !.......  Double check for memory allocation
    IF( JT .GT. JTMAX ) RETURN

    !.......  Store control data table, depending on type of packet
    CALL FILLCDAT( PKTLIST( PKTTYP ), JT, PKTINFO )

    RETURN

END SUBROUTINE FILLCNTL

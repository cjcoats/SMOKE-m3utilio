
SUBROUTINE SETSRCDY( NSRC, SDATE, TZONES, LDAYSAV, MFLAG,&
&                     DAYBEGT, DAYENDT )

!***********************************************************************
!  subroutine SETSRCDY body starts at line < >
!
!  DESCRIPTION:
!
!  PRECONDITIONS REQUIRED:
!     Sets the start and end of the day for all sources.
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!       Version ??/???? by ???
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!***************************************************************************
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
!****************************************************************************

    USE M3UTILIO

!.........  MODULES for public variables
!.........  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY

    IMPLICIT NONE

!...........   SUBROUTINE ARGUMENTS
    INTEGER,      INTENT (IN) :: NSRC             ! no. sources
    INTEGER,      INTENT (IN) :: SDATE            ! start julian date in GMT
    INTEGER,      INTENT (IN) :: TZONES  ( NSRC ) ! time zones per source
    LOGICAL,      INTENT (IN) :: LDAYSAV ( NSRC ) ! true: use daylight time
    LOGICAL,      INTENT (IN) :: MFLAG            ! true: setting for MOVES. false: for MOBILE6
    INTEGER,      INTENT(OUT) :: DAYBEGT( NSRC )  ! start time of SDATE
    INTEGER,      INTENT(OUT) :: DAYENDT( NSRC )  ! end time of SDATE

!...........   Other local variables

    INTEGER       EDATE        ! end julian date of day for source
    INTEGER       ETIME        ! end time HHMMSS
    INTEGER       IOS          ! status of ENVINT
    INTEGER       JDATE        ! tmp Julian date
    INTEGER       S            ! source no.
    INTEGER       STIME        ! start time HHMMSS
    INTEGER       STIME_SET    ! start time HHMMSS when all srcs get same

    LOGICAL       :: DAYLIT = .FALSE. ! true: date is daylight savings
    LOGICAL, SAVE :: FIRSTIME = .TRUE.! true: first time routine called
    LOGICAL, SAVE :: UFLAG  = .FALSE. ! true: all srcs use same day start

    CHARACTER(300) MESG          ! message buffer

    CHARACTER(16) :: PNAME = 'SETSRCDY' ! program name

!***********************************************************************
!   begin body of subroutine SETSRCDY

!.........  Get environment variable for old-style processing in which
!           all of the sources use the same start and end of the day
    IF( FIRSTIME ) THEN

        MESG = 'Start time for using uniform start time'
        STIME_SET = ENVINT ( 'UNIFORM_STIME', MESG, -1, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PNAME,0,0, 'Bad env vble "UNIFORM_STIME"', 2 )
        END IF

        UFLAG = ( STIME_SET > 0 )

        IF( UFLAG ) THEN

            WRITE( MESG,94010 ) 'NOTE: A daily start time of ',&
            &       STIME_SET, 'is being used for all sources.'
            CALL M3MSG2( MESG )

        END IF

        FIRSTIME = .FALSE.

    END IF

!.........  Set start and end time for all sources and return when uniform
!           start time is being used
    IF( UFLAG ) THEN

        EDATE = SDATE
        ETIME = STIME_SET
        CALL NEXTIME( EDATE, ETIME, 230000 )

        DAYBEGT = STIME_SET   ! array
        DAYENDT = ETIME       ! array

        RETURN

    END IF

!.........  Determine if this date is in the range for daylight savings

    DAYLIT = ISDSTIME( SDATE )

!.........  Loop through sources and set start and end date in GMT for
!           each source.  Note that the time zone, adjusted for daylight
!           savings, is the same as the start hour of the day in GMT.0
    DO S = 1, NSRC

        STIME = TZONES( S ) * 10000
        JDATE = SDATE

!.............  If this date is during daylight savings, and if this
!               source is affected by daylight savings
        IF( DAYLIT .AND. LDAYSAV( S ) ) THEN

            CALL NEXTIME( JDATE, STIME, -10000 )

        END IF

!.............  Set start time to  6 A.M. local time for MOBILE6 processing
        IF( CATEGORY == 'MOBILE' .AND. .NOT. MFLAG ) THEN
            CALL NEXTIME( JDATE, STIME, 60000 )
        END IF

!.............  Store start time
        DAYBEGT( S ) = STIME

!.............  Compute end date and time
        EDATE = JDATE
        ETIME = STIME
        CALL NEXTIME( EDATE, ETIME, 230000 )

!.............  Store end time
        DAYENDT( S ) = ETIME

    END DO

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I9, :, 1X ) )

END SUBROUTINE SETSRCDY

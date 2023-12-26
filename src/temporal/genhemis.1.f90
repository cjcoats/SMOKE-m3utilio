
SUBROUTINE GENHEMIS( IGRP, NGRP, NGSZ, JDATE, JTIME, TZONE, DNAME,&
&                     HNAME, PNAME, NAMIN, NAMOUT, EAREAD2D, LDATE )

!***********************************************************************
!  subroutine body starts at line 173
!
!  DESCRIPTION:
!      This subroutine is responsible for applying the temporal profiles
!      to each source using the temporal cross-reference tables.
!      NOTE - Must use EMISV because want to keep EMIS constant for all of the
!      time steps.  This is important to support a different number of
!      day- or hour-specific sources per time step. Many of the IN/OUT
!      arrays are like that because of memory allocation
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created by M. Houyoux 1/99
!
!       Updated by B.H. Baek 5/06 (Re-normalize hourly factors for wildfires)
!
!       Version 07/2014 by C.Coats for  new GENTPRO CSV profiles and cross-references
!
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!************************************************************************
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
!.........  MODSOURC contains the inventory arrays
!........   MODXREF  contains the cross-reference tables
!........   MODTMPRL contains the temporal profile tables
!.........  MODDAYHR contains data for day- and hour-specific data
!.........  MODINFO contains the information about the source category

    USE MODSOURC, ONLY: TZONES, TPFLAG, CSOURC

    USE MODXREF,  ONLY: MDEX, WDEX, DDEX

    USE MODTMPRL, ONLY: NHOLIDAY, HOLJDATE, HOLALTDY, HRLFAC, HRLPROF,&
    &                    NMETPROF, METPROF, IPOL2D, HOUR_TPROF, LTFLAG

    USE MODDAYHR, ONLY: INDXD, INDXH, EMACD, EMACH, NDYSRC, NHRSRC,&
    &                    LDSPOA, LHSPOA, LHPROF,&
    &                    EMAC, EMACV, EMIST, EMFAC, TMAT

    USE MODINFO, ONLY: NSRC, CATEGORY, NIPPA, EACNV, NCHARS

    IMPLICIT NONE

!...........   INCLUDES

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   SUBROUTINE ARGUMENTS

    INTEGER     , INTENT (IN)    :: IGRP      ! pollutant group
    INTEGER     , INTENT (IN)    :: NGRP
    INTEGER     , INTENT (IN)    :: NGSZ      ! Number of pols+emis types
    INTEGER     , INTENT (IN)    :: JDATE     ! Julian date (YYYYDDD) in TZONE
    INTEGER     , INTENT (IN)    :: JTIME     ! Time (HHMMSS) in TZONE
    INTEGER     , INTENT (IN)    :: TZONE     ! Output time zone (typcly 0)
    CHARACTER(*), INTENT (IN)    :: DNAME     ! day-spec file name or NONE
    CHARACTER(*), INTENT (IN)    :: HNAME     ! hour-spec file name or NONE
    CHARACTER(*), INTENT (IN)    :: PNAME     ! met-profile file name or NONE
    CHARACTER(*), INTENT (IN)    :: NAMIN ( NGSZ )      ! inv pol names
    CHARACTER(*), INTENT (IN)    :: NAMOUT( NGSZ )      ! inv pol names
    CHARACTER(*), INTENT (IN)    :: EAREAD2D( NIPPA )   ! tmp inv pol names
    INTEGER     , INTENT (IN OUT):: LDATE                   ! reset previous

!...........   TMAT update variables

    INTEGER, SAVE :: MONTH ( 24, -23:23 )  ! time zone's month 1 ... 12
    INTEGER, SAVE :: DAYOW ( 24, -23:23 )  ! time zone's day-of-week    1 ... 7
    INTEGER, SAVE :: DAYOM ( 24, -23:23 )  ! time zone's day-of-month   1 ... 31

    REAL, ALLOCATABLE, SAVE :: STHOUR( : )         ! episode start hour
    REAL, ALLOCATABLE, SAVE :: EDHOUR( : )         ! episode end hour
    REAL, ALLOCATABLE       :: METVAL( : )         ! tmp met-based temporal factors
    REAL                    :: TMPHRLFAC( 24 )     ! tmp hourly factors

!...........   Other local variables

    INTEGER          C, H, I, II, IS, J, K, K1, K2, KK, L, L2, M, S, V !  indices and counters
    INTEGER          IHR
    INTEGER, SAVE :: TZMIN      ! minimum time zone in inventory
    INTEGER, SAVE :: TZMAX      ! maximum time zone in inventory
    INTEGER, SAVE :: WARNCNT=0  ! warning count
    INTEGER, SAVE :: MXWARN     ! max no. warnings

    INTEGER          MDAY       ! tmp emissions day of month
    INTEGER          WDAY       ! tmp emissions day of week (1=monday)
    INTEGER          HCORR      ! hour correction factor
    INTEGER          HOUR       ! hour of day (1 ... 24)
    INTEGER       :: ST = 0     ! resetting time for episode begin
    INTEGER       :: ED = 0     ! resetting time for episode end
    INTEGER          IOS        ! i/o status
    INTEGER, SAVE :: LTIME = -1 ! time used in previous subroutine call
    INTEGER          MON        ! tmp month number (1=Jan)
    INTEGER          IMET       ! met-based profile ID by source
    INTEGER          VIDX       ! tmp pollutant index for profile
    INTEGER          PIDX       ! tmp pollutant/activity index
    INTEGER          TDATE      ! date for computing time zones update arr
    INTEGER          TTIME      ! time for computing time zones update arr
    INTEGER          EPSBEG     ! tmp episode start hour
    INTEGER          EPSEND     ! tmp episode end hour

    REAL             UFAC            ! tmp units conversion factor
    REAL             TOT             ! tmp total value (denominator)
    REAL          :: NORMFAC   = 0.  ! normalizing factors for hourly factors
    REAL          :: SUMHRLFAC = 0.  ! tmp partial sum of hourly factors
    REAL          :: TOTHRLFAC = 0.  ! tmp sum of hourly factors

    LOGICAL, SAVE :: DFLAG              ! true: day-specific data
    LOGICAL, SAVE :: EFLAG    = .FALSE. ! true: error found
    LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first call to subrtn
    LOGICAL, SAVE :: FIRSTSTP = .TRUE.  ! true: first time step
    LOGICAL, SAVE :: HFLAG              ! true: hour-specific data
    LOGICAL, SAVE :: OUTMSG = .TRUE.    ! true: output message for new day
    LOGICAL       :: RDFLAG = .TRUE.    ! true: read dy data for this iter
    LOGICAL, SAVE :: FIREFLAG=.FALSE.   ! true: read wildfires only
    LOGICAL       :: RHFLAG = .TRUE.    ! true: read hr data for this iter
    LOGICAL, SAVE :: TMATCALC           ! true: need to calculate new TMAT
    LOGICAL, SAVE :: UFLAG  = .FALSE.   ! true: use src-spec hr profiles
    LOGICAL, SAVE :: WKEMSG = .FALSE.   ! true: wkend-profile msg written
    LOGICAL, SAVE :: ZONE4WM            ! true: src zone for week/mon temp prof

    CHARACTER(300)     BUFFER    ! source info buffer
    CHARACTER(300)     MESG      ! message buffer
    CHARACTER(SRCLEN3) CSRC      ! tmp source chars string
    CHARACTER(NAMLEN3) NAMBUF    ! variable name buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'GENHEMIS' ! program name

!***********************************************************************
!   begin body of subroutine GENHEMIS

!.........  If current date is less than previous date or if the date is the
!           same but the time is earlier, then we know that the calling
!           program has started over with a new set of pollutants.  So
!           reset the flag for first timestep.
    IF( JDATE .LT. LDATE .OR. FIRSTIME .OR.&
    &  ( JDATE .EQ. LDATE .AND. JTIME .LT. LTIME )) FIRSTSTP = .TRUE.

!.........  For the first time the subroutine is called,
    IF( FIRSTIME ) THEN

!.............  Check source category name
        IF( CATEGORY .NE. 'AREA'   .AND.&
        &    CATEGORY .NE. 'MOBILE' .AND.&
        &    CATEGORY .NE. 'POINT'        ) THEN

            MESG = 'INTERNAL ERROR: Source category ' // CATEGORY //&
            &       ' is not recognized by program ' // PROGNAME
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        END IF  ! End category selection

!.............  Retrieve environment variables
        MXWARN = ENVINT( WARNSET , 'Maximum warning messages', 100, IOS )
        IF ( IOS .GT. 0 ) THEN
            MESG = 'ERROR:  bad env vble "SMK_MAXWARNING"'
            CALL M3EXIT( PROGNAME, 0,0, MESG, 2 )
        END IF

        MESG = 'Assign weekly/monthly profiles using time zones'
        ZONE4WM = ENVYN ( 'ZONE4WM', MESG, .TRUE., IOS )
        IF ( IOS .GT. 0 ) THEN
            MESG = 'ERROR:  bad env vble "ZONE4WM"'
            CALL M3EXIT( PROGNAME, 0,0, MESG, 2 )
        END IF

!.............  Define the minimum and maximum time zones in the inventory
        TZMIN = MINVAL( TZONES )
        TZMAX = MAXVAL( TZONES )

!.............  Adjust TZMIN and TZMAX for possibility of daylight savings
        TZMIN = TZMIN - 1
        TZMAX = TZMAX + 1

!.............  Ouput hourly emissions in local time
        IF( LTFLAG ) THEN
            TZMIN = 0
            TZMAX = 0
        END IF

!.............  Set flags for daily and hourly data
        DFLAG = ( DNAME .NE. 'NONE' )
        HFLAG = ( HNAME .NE. 'NONE' )

        FIRSTIME = .FALSE.

!.............  Allocate memories for BEGHOUR and ENDHOUR
        ALLOCATE( STHOUR( NDYSRC ),&
        &          EDHOUR( NDYSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EDHOUR', PROGNAME )
        STHOUR = 0.0
        EDHOUR = 0.0

!.............  Define whether processing wildfire or not
        IF( INDEX1( 'ACRESBURNED',NIPPA,EAREAD2D ) > 0 ) THEN
            FIREFLAG = .TRUE.
            CALL M3MSG2( 'Processing Wildfire Emissions....' )
        END IF

    END IF  ! End of first time section

!.........  For new date...

    IF( JDATE .NE. LDATE ) THEN

!.............  Store month and day of week for this output date for
!.............  all source's time zones and for the *next* 24 hours in the
!               simulation
!.............  Use same loop for case where this feature is turned off
        DO I = TZMIN, TZMAX

            TDATE = JDATE
            TTIME = 0     ! Set for loop below from 1 to 24

!.................  Adjust time zone I based on output time zone and correct
!                   for negative value.
            K = I - TZONE
!                IF( K .LT. 0 ) K = 24 + K  !     (e.g., K= -1 -> K= 23)

!.................  Convert output time to local time of I time zone, adjusted
!                   by output time zone.
            CALL NEXTIME( TDATE, TTIME, -K * 10000 )

            DO H = 1, 24

!.....................  When using time zones to set monthly and weekly profiles
!.....................  NOTE - this is more correct
                IF( ZONE4WM ) THEN

                    CALL DAYMON( TDATE, MON, MDAY )   ! get month & day-of-month
                    WDAY   = WKDAY( TDATE )           ! get day-of-week

!.....................  When not using time zones as in previous emissions
!                       systems
                ELSE
                    CALL DAYMON( JDATE, MON, MDAY )
                    WDAY  = WKDAY( JDATE )
                    TDATE = JDATE                  ! set for holiday check

                END IF

!.....................  Check if the date is a holiday.  If so, reset the
!                       day based on the holiday arrays settings.
!.....................  NOTE - this approach will not work for region-specific
!                       setting of holidays.
                J = FIND1( TDATE, NHOLIDAY, HOLJDATE )
                IF( J .GT. 0 ) WDAY = HOLALTDY( J )

                MONTH( H,I ) = MON
                DAYOW( H,I ) = WDAY
                DAYOM( H,I ) = MDAY

                CALL NEXTIME( TDATE, TTIME, 10000 )

            END DO

        END DO

! note: Need to create EMWKDAY and update WRDAYMSG also to use it.
!.............  Write message for day of week and date
        CALL WRDAYMSG( JDATE, MESG )

    END IF   ! End new date check

!.........  Initialize day-specific emissions array with average emissions
!           Only need to do this for a new day because the sources with
!           day-specific data might change for each day.
    EMACV = EMAC  ! array

!.........  If day-specific emissions, prepare day-corrections.
!.........  Read day-specific data for each hour, because different time
!           zones may be used for different sources and this approach
!               is much more workable.
    IF( DFLAG ) THEN

        RDFLAG = .TRUE.
!.................  Read source index for this day
        IF ( .NOT. READ3( DNAME, 'INDXD', ALLAYS3,&
        &                  JDATE, JTIME, INDXD      ) ) THEN
            WRITE( MESG,94010 ) 'WARNING: Could not read "INDXD" '//&
            &     'from file "'// TRIM( DNAME )//'", at', JDATE, ':',&
            &      JTIME
            CALL M3MESG( MESG )

            RDFLAG = .FALSE.
            INDXD = 0   ! array
        END IF      !  if read3() failed on dname

!.............  Read start and ending hour(ENDHOUR) for wildfire processing
        IF( FIREFLAG ) THEN
            IF( .NOT. READ3( DNAME, 'BEGHOUR', ALLAYS3,&
            &              JDATE, JTIME, STHOUR      ) ) THEN
                WRITE(MESG,94010)'WARNING: Could not read  '//&
                &   '"BEGHOUR" from file "'// TRIM( DNAME )//'", at',&
                &   JDATE, ':', JTIME
                CALL M3MESG( MESG )
                FIREFLAG = .FALSE.
            END IF      !  if read3() failed on dname

            IF( .NOT. READ3( DNAME, 'ENDHOUR', ALLAYS3,&
            &              JDATE, JTIME, EDHOUR      ) ) THEN
                WRITE(MESG,94010)'WARNING: Could not read  '//&
                &   '"ENDHOUR" from file "'// TRIM( DNAME )//'", at',&
                &   JDATE, ':', JTIME
                CALL M3MESG( MESG )
                FIREFLAG = .FALSE.
            END IF      !  ifread3() failed on dname
        END IF

    END IF          ! if using day-specific emissions

!.........  Set integer hour of day for output time
    HOUR = 1 + MOD( JTIME / 10000 , 24 )

!.........  Determine if this TMAT needs to be updated
    TMATCALC = ( JDATE .NE. LDATE .OR. FIRSTSTP )

!.........  Construct TMAT -- array of effective composite profile coefficients
    IF( TMATCALC ) THEN

        CALL MKTMAT( NSRC, IGRP, NGRP, NGSZ, JDATE, JTIME, TZONE,&
        &             MONTH, DAYOW, DAYOM, PNAME )

    END IF         ! if TMAT is to be calculated

!.........  If hour-specific emissions, profiles, or activity data...
    IF( HFLAG ) THEN

        RHFLAG = .TRUE.
!.............  Read source index for this hour
        IF( .NOT. READ3( HNAME, 'INDXH', ALLAYS3,&
        &                 JDATE, JTIME, INDXH      ) ) THEN

            WRITE( MESG,94010 ) 'WARNING: Could not read "INDXH" '//&
            &     'from file "'// HNAME//'", at', JDATE, ':', JTIME
            CALL M3MESG( MESG )

            RHFLAG = .FALSE.
            INDXH = 0   ! array

        END IF      !  if read3() failed on HNAME

    END IF  !  End if hour-specific data

!.........  Precompute hour/zone correction factor
    HCORR = TZONE + 23

!.........  Loop through the emissions and emission types and compute the
!           hourly emissions, depending on if the input data is a pollutant
!           or an activity
    DO V = 1, NGSZ

        NAMBUF = NAMIN( V )

!.............  Skip blanks that can occur when NGRP > 1
        IF ( NAMBUF .EQ. ' ' ) CYCLE

!.............  Find pollutant/activity in list of all.  Use EAREAD2D b/c
!               EANAM has been update to contain emission types.
!.............  NOTE - this is sloppy because NIPPA has a larger dimension
!               than EAREAD2D for emission types
        PIDX = INDEX1( NAMBUF, NIPPA, EAREAD2D )

!.............  Set units conversion factor for this pollutant/activity
        UFAC = EACNV( PIDX )

!.............  Read hourly emissions/profile/activity if current pollutant or
!               emission type has them
        IF( LHSPOA( V ) .AND. RHFLAG ) THEN

            IF( .NOT. READ3( HNAME, NAMBUF, ALLAYS3,&
            &                 JDATE, JTIME, EMACH     ) ) THEN

                MESG = 'Could not read "' //&
                &       NAMBUF( 1:LEN_TRIM( NAMBUF ) ) //&
                &       '" from file "' // HNAME //'".'
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

            END IF      !  if read3() failed on hourly data

        END IF  ! if this pollutant is hour-specific

!.............  If source-specific profiles are being used, update temporal
!               matrix for the current hour
        IF( LHPROF( V ) ) THEN

            CALL UPDTMAT( NSRC, IGRP, NGSZ, JDATE, TZONE, V,&
            &              HOUR, MONTH, DAYOW, DAYOM )

        END IF

!.............  For all pollutants and activities...

!.............  Apply hourly factors to all sources for current pollutant or
!               activity. Also apply units conversion.
        EMIST( :,V ) = UFAC * EMACV( :,V ) * TMAT( :,V,HOUR )

!.............  If day-specific data are available for current pollutant
        IF( LDSPOA( V ) .AND. RDFLAG ) THEN

!.................  Read day-specific data
            NAMBUF = NAMIN( V )
            IF( .NOT. READ3( DNAME, NAMBUF, ALLAYS3,&
            &                 JDATE, JTIME, EMACD    )) THEN

                MESG = 'Could not read "' //&
                &       NAMBUF( 1:LEN_TRIM( NAMBUF ) ) //&
                &       '" from file "' // DNAME // '".'
                CALL M3EXIT( PROGNAME,JDATE,JTIME,MESG,2 )

            END IF      !  if read3() failed on day-specific data

!.................  Loop through day-specific sources
            DO I = 1, NDYSRC

                S = INDXD( I )                      ! Get source index
                IF( S .EQ. 0 ) CYCLE                ! If no source, skip
                IF( EMACD( I ) .LE. AMISS3 ) THEN   ! No day-specific emis

!.....................  write out warning message(s) of overwriting with daily emissions
                    IF( WARNCNT <= MXWARN ) THEN
                        WARNCNT = WARNCNT + 1
                        CALL FMTCSRC( CSOURC( S ), NCHARS, BUFFER, L2 )
                        WRITE( MESG,94020 ) 'WARNING: Replacing missing daily '//&
                        &    TRIM(NAMBUF)//' emission with annual-based daily '//&
                        &    'emission "', EMIST(S,V), '" for source:'
                        CALL M3MESG( MESG )
                    END IF
                    CYCLE
                END IF

!.....................  Override annual adjusted emissions with day-specific
!                       emissions and hourly profile adjustments

                K   = 1 + MOD( HOUR + HCORR - TZONES( S ), 24 )
                IHR = HRLPROF( S, DAYOW( HOUR, TZONES( S ) ), PIDX )

!.....................  Re-normalizing hourly temporal factors for wildfires only
                IF( FIREFLAG ) THEN

!.........................  Convert local hours to output time zone
                    EPSBEG = 1 + INT( STHOUR( I ) )/10000
                    EPSEND = 1 + INT( EDHOUR( I ) )/10000
                    K1 = 1 + MOD( EPSBEG + HCORR, 24 )
                    K2 = 1 + MOD( EPSEND + HCORR, 24 )

                    TOTHRLFAC = 0.0
                    SUMHRLFAC = 0.0
                    TMPHRLFAC = 0.0

                    TMPHRLFAC( 1:24 ) = HRLFAC( 1:24,IHR )

!.........................  Sum of 24 hourly temporal factors and store org hourly factors
!                           to tmp hourly factors array
                    DO II = 1, 24
                        TOTHRLFAC = TOTHRLFAC + TMPHRLFAC( II )
                    END DO

!.........................  Sum of hourly temporal factors b/n BENHR and ENDHR
                    IF( K1 > K2 ) THEN

                        WRITE( MESG,94010 )'ERROR: Fire end hour:', K2,&
                        &     ' can not be earlier than begin hour:', K1
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                    ELSE   ! when end hour is greater than start hour

                        DO II = K1, K2
                            SUMHRLFAC = SUMHRLFAC + TMPHRLFAC( II )
                        END DO

                        IF( K1 == 1 .AND. K2 == 1 ) THEN
                            TMPHRLFAC( 2:24 ) = 0.0

                        ELSE IF( K1 == 24 .AND. K2 == 24 ) THEN
                            TMPHRLFAC( 1:23 ) = 0.0

                        ELSE
                            TMPHRLFAC( 1   :K1-1 ) = 0.0   ! resetting to zero
                            TMPHRLFAC( K2+1:  24 ) = 0.0   ! resetting to zero

                        END IF

                    END IF

                    NORMFAC = TOTHRLFAC / SUMHRLFAC

!.........................  Re-normalizing hourly temporal factors
                    DO II = 1, 24
                        TMPHRLFAC( II ) = TMPHRLFAC( II ) * NORMFAC
                    END DO

                    EMIST( S,V ) = UFAC * EMACD( I ) * TMPHRLFAC( K )

                ELSE       ! end of computing wildfires hourly emission factors

                    EMIST( S,V ) = UFAC * EMACD( I ) * HRLFAC( K,IHR )

                END IF

!....................  Apply hour-of-day Gentpro temporal profiles to daily inventory
                VIDX = IPOL2D( V,IGRP )
                IF( VIDX < 1 ) CYCLE
                IMET = METPROF( S,VIDX )

                IF( IMET > 0 .AND. HOUR_TPROF == 'DAY' ) THEN

                    IF( ALLOCATED( METVAL ) ) DEALLOCATE( METVAL )
                    ALLOCATE( METVAL( NMETPROF ), STAT = IOS )
                    CALL CHECKMEM( IOS, 'METVAL', PROGNAME )
                    METVAL = 0.0

                    IF( .NOT. READ3( PNAME, 'DAYTOT', 1, JDATE, JTIME, METVAL ) ) THEN
                        MESG = 'Could not read DAYTOT variable from ' // TRIM( PNAME )
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                    END IF
                    TOT = METVAL( IMET )

                    IF( .NOT. READ3( PNAME, 'HRLSRC', 1, JDATE, JTIME, METVAL ) ) THEN
                        MESG = 'Could not read HRLSRC variable from ' // TRIM( PNAME )
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                    END IF

                    EMIST( S,V ) = UFAC * EMACD( I ) * METVAL( IMET ) / TOT

                END IF

            END DO

        END IF  ! if this pollutant is day-specific

!.............  If hourly data are available for current pollutant/activity and
!               the values are emissions (not profiles), then overwrite with
!               this data
        IF( RHFLAG .AND. LHSPOA( V ) .AND. .NOT. LHPROF( V ) ) THEN
            DO I = 1, NHRSRC
                S = INDXH( I )
                IF( S .EQ. 0 ) CYCLE
                IF( EMACH( I ) .LE. AMISS3 ) THEN
                    IF( WARNCNT <= MXWARN ) THEN
                        WARNCNT = WARNCNT + 1
                        CALL FMTCSRC( CSOURC( S ), NCHARS, BUFFER, L2 )
                        WRITE( MESG,94020 ) 'WARNING: Replacing missing hourly '//&
                        &    TRIM(NAMBUF)//' emission with annual-based hourly '//&
                        &    'emission "', EMIST(S,V), '" for source:'&
                        &    //CRLF()//BLANK10//BUFFER(1:L2)
                        CALL M3MESG( MESG )
                    END IF
                    CYCLE
                END IF
                EMIST( S,V ) = EMACH( I )
            END DO
        END IF

!.............  If input data is in an activity (and output an emission type)...
!               E.G. this is for mobile sources
        IF( NAMIN( V ) .NE. NAMOUT( V ) ) THEN

!.................  Loop through sources and apply emission factors to
!                   hourly activity for non-diurnal emissions
!.................  Apply emission factors tp hourly activity data
!.................  Convert to tons (assuming EFs are in grams)
            EMIST( :,V ) = EMIST( :,V ) * EMFAC( :,V )

        END IF              ! End namin != namout

    END DO                  ! End pollutant loop

!.........  Reset previous date and time for future iterations
    LDATE = JDATE
    LTIME = JTIME

!.........  Set controller to turn off first time step setting
    FIRSTSTP = .FALSE.

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I9, :, 1X ) )
94020 FORMAT( 10( A, :, E10.2, :, 1X ) )
94030 FORMAT( 8X, 'at time ', A8 )

END SUBROUTINE GENHEMIS


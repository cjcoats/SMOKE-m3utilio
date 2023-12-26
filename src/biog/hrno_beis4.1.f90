
SUBROUTINE HRNO_BEIS4( JDATE, JTIME, NX, NY, TA, SOILM, SOILT,&
&              WSAT, WRF_WSAT, ISLTYP, RAIN, GROWAGNO, NGROWAGNO,&
&              NONAGNO, PX_VERSION, INITIAL_HOUR, PTYPE,&
&              PULSEDATE, PULSETIME, EMPOL )

!***********************************************************************
!  subroutine body starts at line  150
!
!  DESCRIPTION:
!
!     Uses new NO algorithm NO = Normalized*Tadj*Padj*Fadj*Cadj
!     to estimate NO emissions
!     Information needed to estimate NO emissions
!     Julian Day          (integer)    JDATE
!     Surface Temperature (MCIP field) TA    (K)
!     Rainfall    (MCIP derived field) RAIN  (cm)
!     Soil Moisture       (MCIP field) SOILM (M**3/M**3) (PX_VERSION)
!          (ratio of volume of water per volume of soil)
!     Soil Temperature    (MCIP field) SOILT (K)         (PX_VERSION)
!     Soil Type           (MCIP field) ISLTYP            (PX_VERSION)
!
!     saturation values for soil types (constants)       (PX_VERSION)
!     FOR PX Version, the Temperature adjustment factor accounts for wet and dry soils
!                and  the precipitation adjustment factor accounts for saturated soils
!     FOR the non-PX version, the basic algorithm remains with a temperature adjustment factor (dry soil)
!                     and no adjustment for saturated soils
!
!
!     The following arrays are updated after each call to HRNO
!     PTYPE   type of NO emission pulse
!     PULSEDATE julian date for the beginning of an NO pulse
!     PULSETIME        time for the beginning of an NO pulse
!
!     The calculation are based on the following paper by J.J. Yienger and H. Levy II
!     J.J. Yienger and H. Levy II, Journal of Geophysical Research, vol 100,11447-11464,1995
!
!     The Temperature Adjustment Factor is based on section 4.2 for wet and dry soils with
!       the following modification (PX version):
!       Instead of classifying soils as either 'wet' or 'dry', the wet and dry adjustment is
!       calculated at each grid cell.  A linear interpolation between the wet and dry adjustment
!       factor is made using the relative amount of soil moisture in the top layer (1cm)
!       as the interpolating factor.  The relative amount of soil moisture is determined by
!       taking the MCIP soil moisture field and dividing by the saturation value defined for each
!       soil type in the PX version of MCIP
!       the soil temperature is used in PX version
!
!     The Precipation Adjustment factor is based on section 4.1 with the following modifications.
!       The rainrate is computed from the MCIP directly using a 24 hr daily total.
!       THe types of Pulses as described in YL95 were used to estimate the NO emission
!       rate.
!
!    Also see the following paper for more information:
!    Proceedings of the Air and Waste Management Association/U.S. Environmental Protection
!    Agency EMission Inventory Conference, Raleigh October 26-28, 1999 Raleigh NC
!    by Tom Pierce and Lucille Bender
!
!    REFERENCES
!
!    JACQUEMIN B. AND NOILHAN J. (1990), BOUND.-LAYER METEOROL., 52, 93-134.
!    J.J. Yienger and H. Levy II, Journal of Geophysical Research, vol 100,11447-11464,1995
!    T. Pierce and L. Bender, Examining the Temporal Variability of Ammonia and Nitric Oxide Emissions from Agricultural Processes
!       Proceedings of the Air and Waste Management Association/U.S. Environmental Protection
!        Agency EMission Inventory Conference, Raleigh October 26-28, 1999 Raleigh NC
!
!  PRECONDITIONS REQUIRED:
!     Normalized NO emissions, Surface Temperature, Soil Moisture, Soil type,
!     NO emission pulse type, soil moisture from previous time step, julian date
!     of NO emission pulse start, time of NO emission pulse start,
!     soil type, SOIL TYPES, Land use data
!
!  SUBROUTINES AND FUNCTIONS CALLED (directly or indirectly):
!     PRECIP_ADJ     computes precipitation adjustment factor
!     FERTILIZER_ADJ computes fertlizer adjustment factor
!     VEG_ADJ        computes vegatation adjustment factor
!     GROWSEASON     computes day of growing season
!
!     PULSETYPE      determines type & duration of NO emission pulse from rainrate
!
!  REVISION  HISTORY:
!    10/01 : Prototype by GAP
!    10/03 : modified transition to non growing season for jul-oct of the year
!    08/04 : Converted to SMOKE code style by C. Seppanen
!    03/22 : minor edits and cleanup by J Vukovich; still differs from
!            CMAQ inline BEIS
!    05/23 : adding logic in case where WSAT variable not available from
!            WRF/MCIP
!    Adapted 10/2023 to USE M3UTILIO by Carlie J. Coats, Jr., UNCIE
!
!***********************************************************************
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
!***********************************************************************

    USE M3UTILIO

    IMPLICIT NONE

!.........  INCLUDES
    INCLUDE 'B3V14DIMS3.EXT'     ! biogenic-related constants

!.........  ARGUMENTS and their descriptions
    INTEGER, INTENT(IN   )  :: JDATE   !  current simulation date (YYYYDDD)
    INTEGER, INTENT(IN   )  :: JTIME   !  current simulation time (HHMMSS)
    INTEGER, INTENT(IN   )  :: NX      !  no. columns
    INTEGER, INTENT(IN   )  :: NY      !  no. rows

    REAL,    INTENT(IN   )  ::  TA    ( NX, NY )    !  air temperature (K)
    REAL,    INTENT(IN   )  ::  SOILM ( NX, NY )    !  soil moisture (m3/m3)
    REAL,    INTENT(IN   )  ::  SOILT ( NX, NY )    !  soil temperature (K)
    REAL,    INTENT(INOUT)  ::  WSAT ( NX, NY )     !  WSAT

    REAL,    INTENT(IN   )  ::  ISLTYP( NX, NY )    !  soil type
    REAL,    INTENT(IN   )  ::  RAIN  ( NX, NY )    !  rainfall rate (cm/ 24 hr)
    REAL,    INTENT(IN   )  ::  GROWAGNO  ( NX, NY )    !  norm NO emissions
    REAL,    INTENT(IN   )  ::  NGROWAGNO ( NX, NY )    !  norm NO emissions
    REAL,    INTENT(IN   )  ::  NONAGNO   ( NX, NY )    !  norm NO emissions

    LOGICAL, INTENT(IN   ) :: PX_VERSION         ! true: using PX version of MCIP
    LOGICAL, INTENT(IN   ) :: INITIAL_HOUR       ! true:
    LOGICAL, INTENT(IN   ) :: WRF_WSAT

    INTEGER, INTENT(INOUT) :: PTYPE (NX, NY)     ! 'pulse' type
    INTEGER, INTENT(INOUT) :: PULSEDATE (NX, NY) ! date of pulse start
    INTEGER, INTENT(INOUT) :: PULSETIME (NX, NY) ! date of pulse end

    REAL,    INTENT(  OUT) :: EMPOL ( NX, NY, NSEF )  !  output pol emissions

!.........  Local PARAMETERS
    REAL,    PARAMETER :: CFNODRYFC = ( 1.0 / 3.0 ) * ( 1.0 / 30.0 )
    REAL,    PARAMETER :: EFAC = EXP( -0.103 * 30.0 )

!............MAXSTYPES and SATURATION array only for use in modified
!............ BEIS4 when WSAT not available from WRFV3-MCIP

    INTEGER, PARAMETER :: MAXSTYPES = 11

! Saturation values for 11 soil types from pxpbl.F  (MCIP PX version)
!       PLEIM-XIU LAND-SURFACE AND PBL MODEL (PX-LSM)
! See JACQUEMIN B. AND NOILHAN J. (1990), BOUND.-LAYER METEOROL., 52,
! 93-134.
!   NOTE ONLY WORKS FOR WRFV3/MCIP output

    REAL, PARAMETER :: SATURATION( MAXSTYPES ) = (/&
    &                  0.395, 0.410, 0.435, 0.485,&
    &                  0.451, 0.420, 0.477, 0.476,&
    &                  0.426, 0.482, 0.482         /)

    CHARACTER(16), PARAMETER :: PROGNAME = 'HRNO_BEIS4'   !  program name

!.........  SCRATCH LOCAL VARIABLES and their descriptions:
    REAL            FAC1, FAC2, FAC3, FAC4
    INTEGER         R, C, L      ! counters
    INTEGER         SOILCAT      ! soil category
    INTEGER         LSM_WATER

    REAL            CFNO         ! NO correction factor
    REAL            CFNOGRASS    ! NO correction factor for grasslands
    REAL            TAIR         ! surface temperature
    REAL            TSOI         ! soil temperature
    REAL         :: CFNOWET, CFNODRY, RATIO

    LOGICAL         USE_SOILT    ! true: use soil temp rather than estimate as in BEIS2
    CHARACTER(256)  MESG         ! message buffer

!***********************************************************************
!   begin body of subroutine HRNO

    USE_SOILT = .TRUE.  ! use soil temperature in PX version
    LSM_WATER = 14

    IF ( .NOT. WRF_WSAT ) THEN
        MESG = "Estimating WSAT in HRNO_BEIS4 routine..."
        CALL M3MESG( MESG )
    ENDIF

    IF ( GROWSEASON( JDATE ) .EQ. 0 ) THEN   ! not growing season

!.........  Loop through cells
        DO R = 1, NY
            DO C = 1, NX
                TAIR = TA( C, R )         ! unit in degree K
!..................  Check max and min bounds for temperature
                IF (TAIR .LT. 200.0) THEN
                    WRITE( MESG, 94010 ) 'TAIR=', TAIR,&
                    &    'out of range at (C,R)=', C, R
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF
                TAIR = MIN( TAIR, 303.0 )

                IF ( TAIR .GT. 268.8690 ) THEN
!...................grass (from BEIS2)
                    CFNO = EXP( 0.04686 * TAIR - 14.30579 )
                ELSE
                    CFNO = 0.0
                END IF

                EMPOL( C,R,NSEF-1 ) = CFNO * ( NGROWAGNO( C,R )  &! agriculture
                &               +   NONAGNO( C,R ) )  ! non-ag

            END DO  ! columns
        END DO  ! rows

    ELSE   ! growing season

        DO R = 1, NY
            DO C = 1, NX

                TAIR = TA( C,R )   ! unit in degree K

! Check min bounds for temperature and limit max to 303 deg K
                IF ( TAIR .LT. 200.0 ) THEN
                    WRITE( MESG, 94010 ) 'TAIR=', TAIR,&
                    &     'out of range at (C,R)=', C, R
                    CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
                END IF
                TAIR = MIN( TAIR, 303.0 )

                IF( .NOT. PX_VERSION ) THEN

!............  Calculate NO emissions by going thru temperature cases
!............  Calculate soil temp from air temp

                    TSOI = 0.72 * TAIR + 82.28

!............  TSOI convert to deg C
                    TSOI = MIN( MAX( TSOI, 273.16 ), 303.16 ) - 273.16
                    CFNODRY = CFNODRYFC * TSOI ! see YL 1995 Eqn 9a p11452
                    IF ( TSOI .LE. 10.0 ) THEN ! see YL 1995 Eqn 7b
                        CFNOWET = 0.28 * EFAC * TSOI   ! linear cold case
                    ELSE
                        CFNOWET = EFAC * EXP( 0.103 * TSOI )   !
                    END IF
                    CFNO = 0.5 * ( CFNOWET + CFNODRY )

                    FAC1 = GROWAGNO( C,R ) * CFNO&
                    &     * FERTILIZER_ADJ( JDATE )&
                    &     * VEG_ADJ( JDATE )
!!!!!!!
                    IF ( INITIAL_HOUR ) THEN
                        FAC2 = 1.0
                        PTYPE( C,R ) = 0
                        PULSEDATE( C,R ) = 0
                        PULSETIME( C,R ) = 0
                    ELSE
                        FAC2 = PRECIP_ADJ( JDATE, JTIME, RAIN( C,R ),&
                        &                   PTYPE( C,R ), PULSEDATE( C,R ),&
                        &                   PULSETIME( C,R ) )
                    END IF

                ELSE   ! PX version being used

                    SOILCAT = ISLTYP( C,R )
                    TSOI = 0.0

!............ case where WSAT not available from WRF/MCIP
                    IF ( .NOT. WRF_WSAT ) THEN
                        WSAT( C,R ) =  SATURATION( SOILCAT )
                    ENDIF

                    IF( SOILCAT .NE. LSM_WATER  ) THEN
                        RATIO = SOILM( C,R ) / WSAT( C,R )
                        IF ( USE_SOILT ) THEN
                            TSOI = SOILT( C,R )
                        ELSE
                            TSOI = 0.72 * TAIR + 82.28
                        END IF
!................... Convert to deg C
                        TSOI = MIN( MAX( TSOI, 273.16 ), 303.16 ) - 273.16

                        CFNODRY = CFNODRYFC * TSOI ! see YL 1995 Eqn 9a
                        IF ( TSOI .LE. 10.0 ) THEN ! see YL 1995 Eqn 7b
                            CFNOWET = 0.28 * EFAC * TSOI   ! linear cold
                        ELSE
                            CFNOWET = EFAC * EXP( 0.103 * TSOI )   !
                        END IF
                        CFNO = CFNODRY + RATIO * ( CFNOWET - CFNODRY )
                        FAC1 = GROWAGNO( C,R ) * CFNO&
                        &     * FERTILIZER_ADJ( JDATE )&
                        &     * VEG_ADJ( JDATE )
                    ELSE
                        FAC1 = 0.0
                    END IF   ! if block for ne water

                    IF ( INITIAL_HOUR ) THEN
                        FAC2 = 1.0
                        PTYPE( C,R ) = 0
                        PULSEDATE( C,R ) = 0
                        PULSETIME( C,R ) = 0
                    ELSE
                        FAC2 = PRECIP_ADJ_PX( JDATE, JTIME, RAIN( C,R ),&
                        & SOILM( C,R ), WSAT( C,R ), SOILCAT,&
                        & PTYPE( C,R ), PULSEDATE( C,R ), PULSETIME( C,R ) )
                    END IF

                END IF   ! PX check

                IF ( TAIR .GT. 268.8690 ) THEN
                    CFNOGRASS = EXP( 0.04686 * TAIR - 14.30579 )
                    FAC3 = NGROWAGNO( C,R ) * CFNOGRASS
                    FAC4 = NONAGNO( C,R ) * CFNOGRASS
                ELSE
                    FAC3 = 0.0
                    FAC4 = 0.0
                END IF

                IF( TSOI .LE. 0.0 ) THEN
                    EMPOL( C,R,NSEF-1 ) = 0.0
                ELSE
                    EMPOL( C,R,NSEF-1 ) = MAX(( FAC1*FAC2 ), FAC3) + FAC4
                END IF

            END DO  ! columns
        END DO  ! rows

    ENDIF   ! if block for growing or not growing season

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( A, F10.2, 1X, A, I3, ',', I3 )
94020 FORMAT( A, F10.2, 1X, A, I3, ',', I3, A )

    RETURN



CONTAINS


    REAL FUNCTION PRECIP_ADJ_PX( JDATE, JTIME, RNTOT, SMO,&
    &              WST, SCAT, PTYPE, PULSEDATE, PULSETIME )

!***********************************************************************
!  function body starts at line  386
!
!  DESCRIPTION:
!
!     computes precipitation adjustment factor for estimate of NO emissions
!     uses  julian day, time, soil moisture
!     requires the use of three arrays that are re-used each time step
!     PTYPE, PULSEDATE, PULSETIME
!     These arrays store the type of NO pulse initiated by the rainfall
!     and the starting date and time of the pulse.
!
!  PRECONDITIONS REQUIRED:
!     Soil Moisture current time, Soil Moisture previous time,
!     Soil type, Land Use, PTYPE, PULSEDATE, PULSETIME
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!     precipfact     computes precip adjustment factor from rainrate and time since pulse initiation
!     pulsetype      determines type & duration of NO emission pulse from rainrate
!
!  REVISION  HISTORY:
!    11/01 : Prototype by GAP
!    3/05  : create separate functions for PX vs non-PX versions
!    3/22  : adapted from CMAQ hrno.F logic; adding WSAT
!
!***********************************************************************

        IMPLICIT NONE

!.............  Function arguments
        INTEGER, INTENT(IN) :: JDATE
        INTEGER, INTENT(IN) :: JTIME
        INTEGER, INTENT(IN) :: SCAT

        REAL, INTENT(IN) :: RNTOT
        REAL, INTENT(IN) :: SMO   ! only avilable if PX version
        REAL, INTENT(IN) :: WST

        INTEGER, INTENT(IN OUT) :: PTYPE     ! pulse type
        INTEGER, INTENT(IN OUT) :: PULSEDATE ! date of pulse start
        INTEGER, INTENT(IN OUT) :: PULSETIME ! date of pulse end

!.............  Local parameters
        REAL,    PARAMETER :: SAT_THRES = 0.95
        INTEGER, PARAMETER :: LSM_WATER = 14

!.............  Local variables
        INTEGER PTYPE_TEST

!-----------------------------------------------------------------------------

!.............  Summary of algorithm
!        1. compute rate of change of soil moisture from soil moisture
!        2. estimate rainrate from soilmoisture and soil moisture rate
!        3. compute adjustment using pulsetype, rainrate,ptype, and date/time
!             if stronger NO pulse compared to previous time step , then
!             start a new NO emission pulse
!             otherwise continue present NO pulse
!        4. override adjustment for saturated soils

        IF ( SCAT .NE. LSM_WATER  ) THEN
            IF ( SMO .GE. SAT_THRES * WST ) THEN
                PRECIP_ADJ_PX = 0.0
            ELSE
                PTYPE_TEST = PULSETYPE( RNTOT )
!............... Rainfall class type increase
                IF ( PTYPE_TEST .GT. PTYPE ) THEN
!............... NO emissions pulse generated

                    PULSEDATE = JDATE
                    PULSETIME = JTIME
                    PTYPE = PTYPE_TEST
                END IF

                PRECIP_ADJ_PX = PRECIPFAC( JDATE, JTIME, PULSEDATE,&
                &   PULSETIME, PTYPE )

            END IF

        ELSE   ! if water

            PRECIP_ADJ_PX = 0.0

        END IF

        RETURN

    END FUNCTION PRECIP_ADJ_PX

    REAL FUNCTION PRECIP_ADJ( JDATE, JTIME, RAIN,&
    &                          PTYPE, PULSEDATE, PULSETIME )

!***********************************************************************
!  function body starts at line  386
!
!  DESCRIPTION:
!
!     computes precipitation adjustment factor for estimate of NO emissions
!     uses  julian day, time, soil moisture
!     requires the use of three arrays that are re-used each time step
!     PTYPE, PULSEDATE, PULSETIME
!     These arrays store the type of NO pulse initiated by the rainfall
!     and the starting date and time of the pulse.
!
!  PRECONDITIONS REQUIRED:
!     Soil Moisture current time, Soil Moisture previous time,
!     Soil type, Land Use, PTYPE, PULSEDATE, PULSETIME
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!     precipfact     computes precip adjustment factor from rainrate and time since pulse initiation
!     pulsetype      determines type & duration of NO emission pulse from rainrate
!
!  REVISION  HISTORY:
!    11/01 : Prototype by GAP
!    3/05  : created a non-PX version of this function
!    3/22  : taken from CMAQ hrno.F
!
!***********************************************************************

        IMPLICIT NONE

!.............  Function arguments
        INTEGER, INTENT(IN) :: JDATE
        INTEGER, INTENT(IN) :: JTIME

        REAL, INTENT(IN) :: RAIN

        INTEGER, INTENT(IN OUT) :: PTYPE     ! pulse type
        INTEGER, INTENT(IN OUT) :: PULSEDATE ! date of pulse start
        INTEGER, INTENT(IN OUT) :: PULSETIME ! date of pulse end

!.............  Local variable
        INTEGER PTYPE_TEST

!-----------------------------------------------------------------------------

! Summary of algorithm
!    1. if no rainfall or new rainfall class less than current one,
!       continue existing NO emission pulse
!    2. if new rainfall that increases rainfall class, then create new
!        NO emission pulse using pulsetype, rainrate, ptype, and date/time -
!        if stronger NO pulse compared to previous time step, then start
!        a new NO emission pulse

        PTYPE_TEST = PULSETYPE( RAIN )
!.............Rainfall class type increases
        IF ( PTYPE_TEST .GT. PTYPE ) THEN
!.............NO emissions pulse generated
            PULSEDATE = JDATE
            PULSETIME = JTIME
            PTYPE = PTYPE_TEST
        END IF

        PRECIP_ADJ = PRECIPFAC( JDATE, JTIME, PULSEDATE, PULSETIME,&
        &            PTYPE )

        RETURN

    END FUNCTION PRECIP_ADJ

!-----------------------------------------------------------------------------

!.............  This internal function computes a fertilizer adjustment factor
!               for the given date in yyyyddd format. If it is not growing
!               season, the adjustment factor is 0; otherwise, it ranges from
!               0.0 to 1.0.
    REAL FUNCTION FERTILIZER_ADJ( DATE )

        IMPLICIT NONE

!.............  Function arguments
        INTEGER, INTENT(IN) :: DATE

!.............  Local variables
        INTEGER  GDAY

!-----------------------------------------------------------------------------

        GDAY = GROWSEASON( DATE )

        IF( GDAY == 0 ) THEN
            FERTILIZER_ADJ = 0.
        ELSE IF( GDAY >= 1 .AND. GDAY < 30 ) THEN   ! first month of growing season
            FERTILIZER_ADJ = 1.
        ELSE IF( GDAY >= 30 ) THEN
            FERTILIZER_ADJ = 1. + 30. / 184. - FLOAT( GDAY ) / 184.
        ELSE
            WRITE( MESG,94010 ) 'Invalid date specified; date = ',&
            &                    DATE, 'growing season day = ',&
            &                    GDAY
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010   FORMAT( A, F10.2, 1X, A, I3, ',', I3 )
94020   FORMAT( A, F10.2, 1X, A, I3, ',', I3, A )


        RETURN

    END FUNCTION FERTILIZER_ADJ

!-----------------------------------------------------------------------------

!.............  This internal function computes a vegetation adjustment factor
!               for the given date in yyyyddd format. The adjustment factor
!               ranges from 0.5 to 1.0.
    REAL FUNCTION VEG_ADJ( DATE )

        IMPLICIT NONE

!.............  Function arguments
        INTEGER, INTENT(IN) :: DATE

!.............  Local variables
        INTEGER  GDAY

!-----------------------------------------------------------------------------

        GDAY = GROWSEASON( DATE )

        IF( GDAY <= 30 ) THEN
            VEG_ADJ = 1.
        ELSE IF( GDAY > 30 .AND. GDAY < 60 ) THEN
            VEG_ADJ = 1.5 - ( FLOAT( GDAY ) / 60. )
        ELSE IF( GDAY >= 60 ) THEN
            VEG_ADJ = 0.5
        ELSE
            WRITE( MESG,94010 ) 'Invalid date specified; date = ',&
            &                    DATE, 'growing season day = ',&
            &                    GDAY
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010   FORMAT( A, F10.2, 1X, A, I3, ',', I3 )
94020   FORMAT( A, F10.2, 1X, A, I3, ',', I3, A )

        RETURN

    END FUNCTION VEG_ADJ

!-----------------------------------------------------------------------------

!.............  This internal function computes the day of the growing season
!               corresponding to the given date in yyyyddd format.
    INTEGER FUNCTION GROWSEASON( DATE )

        IMPLICIT NONE

!.............  Function arguments
        INTEGER, INTENT(IN) :: DATE

!.............  Local parameters
        INTEGER, PARAMETER :: GSEASON_START = 0401
        INTEGER, PARAMETER :: GSEASON_END   = 1031

!.............  Local variables
        INTEGER  YEAR, MONTH, DAY
        INTEGER  JDAY, GDAY
        INTEGER  GSTART_MONTH, GSTART_DAY, GSJULIAN_START
        INTEGER  GEND_MONTH, GEND_DAY, GSJULIAN_END

!-----------------------------------------------------------------------------

        YEAR = INT( FLOAT( DATE ) / 1000. )
        JDAY = DATE - YEAR * 1000

        GSTART_MONTH = INT( FLOAT( GSEASON_START ) / 100. )
        GSTART_DAY   = GSEASON_START - GSTART_MONTH * 100
        GSJULIAN_START = JULIAN( YEAR, GSTART_MONTH, GSTART_DAY )

        GEND_MONTH = INT( FLOAT( GSEASON_END ) / 100. )
        GEND_DAY   = GSEASON_END - GEND_MONTH * 100
        GSJULIAN_END = JULIAN( YEAR, GEND_MONTH, GEND_DAY )

        IF( JDAY >= GSJULIAN_START .AND. JDAY <= GSJULIAN_END ) THEN  ! growing season
            GROWSEASON = JDAY - GSJULIAN_START + 1
        ELSE IF( JDAY >= 1 .AND. JDAY <= 366 ) THEN  ! before or after growing season
            GROWSEASON = 0
        ELSE
            WRITE( MESG,94010 ) 'Invalid date specified; date = ',&
            &                    DATE, 'jday = ', JDAY
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010   FORMAT( A, F10.2, 1X, A, I3, ',', I3 )
94020   FORMAT( A, F10.2, 1X, A, I3, ',', I3, A )

        RETURN

    END FUNCTION GROWSEASON

    REAL FUNCTION PRECIPFAC( JDATE, JTIME, PDATE, PTIME, PTYPE )

! Compute a precipitation adjustment factor from a previous 24 hour
! rainfall based on YL 1995.
! The pulse type is an integer ranging from 0 to 3 indicating the type
! of rainfall rate:
! If rainfall < 0.1 cm in last 24 hr, "reset"
! Else if rainfall < 0.5 cm in last 24 hr, and time since last pulse is
! .ge. 2 days,
!    reset; else, precipfact=11.19*...
! Else if rainfall < 1.5 cm in last 24 hr, and time since last pulse is
! .ge. 6 days,
!    reset; else, precipfact=14.68*...
! Else if rainfall >=1.5 cm in last 24 hr, and time since last pulse is
! .ge. 13 days,
!    reset; else, precipfact=18.46*...

        IMPLICIT NONE

! Function arguments:
        INTEGER, INTENT( IN )    :: JDATE, JTIME, PDATE, PTIME
        INTEGER, INTENT( INOUT ) :: PTYPE

! Parameters:
!............... daypersec below = 0.000011574074074
        REAL, PARAMETER :: DAYPERSEC = 1.0 / ( 24.0 * 3600.0 )

! Local variables:
        REAL DAYDIFF, DAYDIF1

!-----------------------------------------------------------------------

        DAYDIFF = FLOAT( SECSDIFF( PDATE, PTIME, JDATE, JTIME ) ) *&
        &              DAYPERSEC
        DAYDIF1 = DAYDIFF + 1.0

        SELECT CASE( PTYPE )
          CASE( 0 )
            PRECIPFAC = 1.0
          CASE( 1 )
            IF ( ( DAYDIFF ) .LT. 2.0 ) THEN
                PRECIPFAC = 11.19 * EXP( -0.805 * DAYDIF1 )
            ELSE
                PTYPE = 0
                PRECIPFAC = 1.0
            END IF
          CASE( 2 )
            IF ( ( DAYDIFF ) .LT. 6.0 ) THEN
                PRECIPFAC = 14.68 * EXP( -0.384 * DAYDIF1 )
            ELSE
                PTYPE = 0
                PRECIPFAC = 1.0
            END IF
          CASE( 3 )
            IF ( ( DAYDIFF ) .LT. 13.0 ) THEN
                PRECIPFAC = 18.46 * EXP( -0.208 * DAYDIF1 )
            ELSE
                PTYPE = 0
                PRECIPFAC = 1.0
            END IF
          CASE DEFAULT
            WRITE( MESG,'( A, I6 )' ) 'Invalid Pulse Type specified ',&
            &                    PTYPE
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END SELECT

        RETURN

    END FUNCTION PRECIPFAC

!-----------------------------------------------------------------------------

!.............  This internal function computes the pulse type from a rainfall
!               rate (see YL 1995).
    INTEGER FUNCTION PULSETYPE( RAINRATE )

        IMPLICIT NONE

!.............  Function arguments
        REAL, INTENT(IN) :: RAINRATE

!-----------------------------------------------------------------------------

        IF( RAINRATE < 0.1 ) THEN
            PULSETYPE = 0
        ELSE IF( RAINRATE < 0.5 ) THEN
            PULSETYPE = 1
        ELSE IF( RAINRATE < 1.5 ) THEN
            PULSETYPE = 2
        ELSE
            PULSETYPE = 3
        END IF

        RETURN

    END FUNCTION PULSETYPE

!-----------------------------------------------------------------------------

END SUBROUTINE HRNO_BEIS4


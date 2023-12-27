
SUBROUTINE HRNO( JDATE, JTIME, NX, NY, TA, SOILM, SOILT,        &
                 ISLTYP, RAIN, GROWAGNO, NGROWAGNO, NONAGNO,    &
                 PX_VERSION, INITIAL_HOUR, PTYPE, PULSEDATE,    &
                 PULSETIME, EMPOL )

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
    !         (ratio of volume of water per volume of soil)
    !     Soil Temperature    (MCIP field) SOILT (K)         (PX_VERSION)
    !     Soil Type           (MCIP field) ISLTYP            (PX_VERSION)
    !
    !     saturation values for soil types (constants)       (PX_VERSION)
    !     FOR PX Version, the Temperature adjustment factor accounts for wet and dry soils
    !         and  the precipitation adjustment factor accounts for saturated soils
    !     FOR the non-PX version, the basic algorithm remains with a temperature adjustment factor (dry soil)
    !         and no adjustment for saturated soils
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
    !     PRECIPFACT     computes precip adjustment factor from rainrate and time since pulse initiation
    !     PULSETYPE      determines type  duration of NO emission pulse from rainrate
    !
    !  REVISION  HISTORY:
    !       10/01 : Prototype by GAP
    !       10/03 : modified transition to non growing season for jul-oct of the year
    !       08/04 : Converted to SMOKE code style by C. Seppanen
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***********************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !         System
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

    !.......  INCLUDES
    INCLUDE 'B3V14DIMS3.h90'         ! biogenic-related constants

    !.......  ARGUMENTS and their descriptions
    INTEGER, INTENT (IN)  :: JDATE       !  current simulation date (YYYYDDD)
    INTEGER, INTENT (IN)  :: JTIME       !  current simulation time (HHMMSS)
    INTEGER, INTENT (IN)  :: NX          !  no. columns
    INTEGER, INTENT (IN)  :: NY          !  no. rows

    REAL, INTENT (IN)  ::  TA    ( NX, NY )        !  air temperature (K)
    REAL, INTENT (IN)  ::  SOILM ( NX, NY )        !  soil moisture (m3/m3)
    REAL, INTENT (IN)  ::  SOILT ( NX, NY )        !  soil temperature (K)
    REAL, INTENT (IN)  ::  ISLTYP( NX, NY )        !  soil type
    REAL, INTENT (IN)  ::  RAIN  ( NX, NY )        !  rainfall rate (cm/ 24 hr)
    REAL, INTENT (IN)  ::  GROWAGNO  ( NX, NY )        !  norm NO emissions
    REAL, INTENT (IN)  ::  NGROWAGNO ( NX, NY )        !  norm NO emissions
    REAL, INTENT (IN)  ::  NONAGNO   ( NX, NY )        !  norm NO emissions

    LOGICAL, INTENT (IN) :: PX_VERSION             ! true: using PX version of MCIP
    LOGICAL, INTENT (IN) :: INITIAL_HOUR           ! true:

    INTEGER, INTENT (IN OUT) :: PTYPE (NX, NY)         ! 'pulse' type
    INTEGER, INTENT (IN OUT) :: PULSEDATE (NX, NY)     ! date of pulse start
    INTEGER, INTENT (IN OUT) :: PULSETIME (NX, NY)     ! date of pulse end


    REAL, INTENT (OUT) :: EMPOL ( NX, NY, NSEF )      !  output pol emissions
    REAL TMP1, TMP2, TMP3, TMP4

    !.......  Local PARAMETERS
    INTEGER, PARAMETER :: MAXSTYPES = 11

    !   Saturation values for 11 soil types from pxpbl.F  (MCIP PX version)
    !   PLEIM-XIU LAND-SURFACE AND PBL MODEL (PX-LSM)
    !   See JACQUEMIN B. AND NOILHAN J. (1990), BOUND.-LAYER METEOROL., 52, 93-134.
    REAL, PARAMETER :: SATURATION( MAXSTYPES ) =    &
                   (/ 0.395, 0.410, 0.435, 0.485,   &
                      0.451, 0.420, 0.477, 0.476,   &
                      0.426, 0.482, 0.482        /)

    CHARACTER(16), PARAMETER :: PROGNAME = 'HRNO'       !  program name

    !.......  SCRATCH LOCAL VARIABLES and their descriptions:
    INTEGER         R, C, L          ! counters
    INTEGER         SOILCAT          ! soil category

    REAL            CFNO             ! NO correction factor
    REAL            CFNOGRASS        ! NO correction factor for grasslands
    REAL            TAIR             ! surface temperature
    REAL            TSOI             ! soil temperature
    REAL            CFNOWET, CFNODRY, RATIO

    LOGICAL         USE_SOILT        ! true: use soil temp rather than estimate as in BEIS2

    CHARACTER(256)  MESG             ! message buffer

    !***********************************************************************
    !   begin body of subroutine HRNO

    USE_SOILT = .TRUE.      ! use soil temperature in PX version

    !.......  Loop through cells
    DO R = 1, NY
    DO C = 1, NX

        TAIR = TA( C, R )                     ! unit in degree K

        !......  Check max and min bounds for temperature
        IF (TAIR < 200.0) THEN
            WRITE( MESG, 94010 ) 'TAIR=', TAIR, 'out of range at (C,R)=', C, R
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
        END IF

        IF (TAIR > 315.0 ) THEN
            WRITE( MESG, 94020 ) 'TAIR=', TAIR, 'out of range at (C,R)=', C, R, &
                ' resetting to 315K'
            CALL M3WARN( PROGNAME, JDATE, JTIME, MESG )
            TAIR = 315.0
        END IF

        IF( .NOT. PX_VERSION ) THEN

            !.......  Calculate NO emissions by going thru temperature cases
            IF( GROWSEASON( JDATE ) == 0) THEN                     ! not growing season
                IF( TAIR > 303.00 ) TAIR = 303.00

                IF ( TAIR > 268.8690 ) THEN
                    CFNO = EXP( 0.04686 * TAIR - 14.30579 )             ! grass (from BEIS2)
                ELSE
                    CFNO = 0.0
                END IF

                EMPOL( C,R,NSEF-1 ) = NGROWAGNO( C,R ) * CFNO + NONAGNO( C,R ) * CFNO               !  non-agriculture

            ELSE               ! growing season

                TSOI = 0.72 * TAIR + 82.28
                IF (TSOI <= 273.16) TSOI = 273.16
                IF (TSOI >= 303.16) TSOI = 303.16

                CFNODRY = (1./3.) * (1./30.) * (TSOI-273.16)      ! see YL 1995 Equa 9a p. 11452
                IF (TSOI <= 283.16) THEN             ! linear cold case
                    CFNOWET = (TSOI-273.16)*EXP(-0.103*30.0)*0.28     ! see YL 1995 Equ 7b
                ELSE                                 ! exponential case
                    CFNOWET = EXP(0.103 * (TSOI-273.16)) * EXP(-0.103 * 30.0)
                END IF
                CFNO = 0.5 * CFNOWET + 0.5 * CFNODRY

                IF ( TAIR > 303.00 ) TAIR = 303.00

                IF ( TAIR > 268.8690 ) THEN
                    CFNOGRASS = EXP( 0.04686 * TAIR - 14.30579 )     ! grass (from BEIS2)
                ELSE
                    CFNOGRASS = 0.0
                END IF




                TMP1 = GROWAGNO( C,R ) * CFNO * FERTILIZER_ADJ( JDATE ) * VEG_ADJ( JDATE )

                IF (INITIAL_HOUR) THEN
                    TMP2 = 1.0
                    PTYPE(C,R) = 0
                    PULSEDATE(C,R) = 0
                    PULSETIME(C,R) = 0
                ELSE
                    TMP2 = PRECIP_ADJ( JDATE, JTIME, RAIN( C,R ),       &
                                       PTYPE( C,R ), PULSEDATE( C,R ),  &
                                       PULSETIME( C,R ) )
                ENDIF
                TMP3 = NGROWAGNO( C,R ) * CFNOGRASS

                TMP4 =  NONAGNO( C,R ) * CFNOGRASS


                EMPOL( C,R,NSEF-1 )=MAX((TMP1*TMP2),TMP3)+TMP4


            END IF

        !.......  If using PX version of MCIP
        ELSE
            IF( GROWSEASON( JDATE ) == 0) THEN             ! not growing season
                IF( TAIR > 303.00 ) TAIR = 303.00

                IF ( TAIR > 268.8690 ) THEN
                    CFNO = EXP( 0.04686 * TAIR - 14.30579 )     ! grass (from BEIS2)
                ELSE
                    CFNO = 0.0
                END IF

                EMPOL( C,R,NSEF-1 ) = NGROWAGNO( C,R ) * CFNO + NONAGNO( C,R ) * CFNO       !  non-agriculture

            ELSE       ! growing season

                IF( USE_SOILT ) THEN
                    TSOI = SOILT( C,R )
                ELSE
                    TSOI = 0.72 * TAIR + 82.28
                END IF

                IF (TSOI <= 273.16) TSOI = 273.16
                IF (TSOI >= 303.16) TSOI = 303.16

                CFNODRY = (1./3.)*(1./30.)*(TSOI-273.16)      ! see YL 1995 Equa 9a p. 11452
                IF (TSOI <= 283.16) THEN             ! linear cold case
                    CFNOWET = (TSOI-273.16)*EXP(-0.103*30.0)*0.28     ! see YL 1995 Equ 7b
                ELSE                                 ! exponential case
                    CFNOWET = EXP(0.103 * (TSOI-273.16)) * EXP(-0.103 * 30.0)
                END IF

                SOILCAT = INT( ISLTYP( C,R ) )
                IF( SOILCAT > 0. .AND. SOILCAT <= MAXSTYPES ) THEN
                    RATIO = SOILM( C,R ) / SATURATION( SOILCAT )
                    CFNO = RATIO*CFNOWET + (1.-RATIO)*CFNODRY
                ELSE
                    CFNO = 0.
                END IF

                IF ( TAIR > 303.00 ) TAIR = 303.00

                IF ( TAIR > 268.8690 ) THEN
                    CFNOGRASS = EXP( 0.04686 * TAIR - 14.30579 )     ! grass (from BEIS2)
                ELSE
                    CFNOGRASS = 0.
                END IF
                TMP1 = GROWAGNO( C,R ) * CFNO * FERTILIZER_ADJ( JDATE ) * VEG_ADJ( JDATE )

                TMP2 = PRECIP_ADJ_PX( JDATE, JTIME, RAIN( C,R ),        &
                                      SOILM( C,R ), ISLTYP( C,R ),      &
                                      PTYPE( C,R ), PULSEDATE( C,R ),   &
                                      PULSETIME( C,R ) )
                TMP3 = NGROWAGNO( C,R ) * CFNOGRASS

                TMP4 =  NONAGNO( C,R ) * CFNOGRASS


                EMPOL( C,R,NSEF-1 )=MAX((TMP1*TMP2),TMP3)+TMP4

            END IF      ! growing season check
        END IF      ! PX version check
    END DO      ! loop over columns
    END DO      ! loop over rows
    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( A, F10.2, 1X, A, I3, ',', I3 )
94020 FORMAT( A, F10.2, 1X, A, I3, ',', I3, A )

    RETURN



CONTAINS


    REAL FUNCTION PRECIP_ADJ_PX( JDATE, JTIME, RAIN, SOILM, ISLTYP,&
                                 PTYPE, PULSEDATE, PULSETIME )

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
        !     pulsetype      determines type  duration of NO emission pulse from rainrate
        !
        !  REVISION  HISTORY:
        !    11/01 : Prototype by GAP
        !    3/05  : create separate functions for PX vs non-PX versions
        !
        !***********************************************************************

        IMPLICIT NONE

        !.......  Function arguments
        INTEGER, INTENT(IN) :: JDATE
        INTEGER, INTENT(IN) :: JTIME

        REAL, INTENT(IN) :: RAIN
        REAL, INTENT(IN) :: SOILM           ! only avilable if PX version
        REAL, INTENT(IN) :: ISLTYP          ! only available if PX version

        INTEGER, INTENT(IN OUT) :: PTYPE             ! pulse type
        INTEGER, INTENT(IN OUT) :: PULSEDATE         ! date of pulse start
        INTEGER, INTENT(IN OUT) :: PULSETIME         ! date of pulse end

        !	     PULSETYPE, PRECIPFACT

        !.......  Local parameters
        REAL, PARAMETER :: SAT_THRES = 0.95
        REAL, PARAMETER :: IRFACT = 2.0

        !.......  Local variables
        INTEGER  SOILCAT             ! soil type category

        !-----------------------------------------------------------------------------

        !.......  Summary of algorithm
        !        1. compute rate of change of soil moisture from soil moisture
        !        2. estimate rainrate from soilmoisture and soil moisture rate
        !        3. compute adjustment using pulsetype, rainrate,ptype, and date/time
        !         if stronger NO pulse compared to previous time step , then
        !         start a new NO emission pulse
        !         otherwise continue present NO pulse
        !        4. override adjustment for saturated soils

        SOILCAT = INT( ISLTYP )

        IF( INITIAL_HOUR ) THEN
            PTYPE = 0
            PULSEDATE = 0
            PULSETIME = 0

            IF( SOILCAT > 0 .AND. SOILCAT <= MAXSTYPES ) THEN
                PRECIP_ADJ_PX = 2.
            ELSE
                PRECIP_ADJ_PX = 1.
            END IF

        ELSE

            IF( PULSETYPE( RAIN ) <= PTYPE ) THEN

                !.......  No new rainfall or new rainfall class less than current one
                PRECIP_ADJ_PX = PRECIPFACT( PTYPE, JDATE, JTIME, PULSEDATE, PULSETIME )

            ELSE

                !.......  Rainfall class type increases (NO emission pulse generated)
                PULSEDATE = JDATE
                PULSETIME = JTIME
                PTYPE = PULSETYPE( RAIN )

                PRECIP_ADJ_PX = PRECIPFACT( PTYPE, JDATE, JTIME, PULSEDATE, PULSETIME )
            END IF


            IF( SOILCAT > 0. .and. SOILCAT <= MAXSTYPES ) THEN
                IF( SOILM >= SAT_THRES * SATURATION( SOILCAT ) ) THEN
                    PRECIP_ADJ_PX = 0.
                END IF
            END IF


        END IF

        RETURN

    END FUNCTION PRECIP_ADJ_PX

    REAL FUNCTION PRECIP_ADJ( JDATE, JTIME, RAIN, PTYPE, PULSEDATE, PULSETIME )

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
        !     pulsetype      determines type  duration of NO emission pulse from rainrate
        !
        !  REVISION  HISTORY:
        !    11/01 : Prototype by GAP
        !    3/05  : created a non-PX version of this function
        !
        !***********************************************************************

        IMPLICIT NONE

        !.......  Function arguments
        INTEGER, INTENT(IN) :: JDATE
        INTEGER, INTENT(IN) :: JTIME

        REAL, INTENT(IN) :: RAIN


        INTEGER, INTENT(IN OUT) :: PTYPE             ! pulse type
        INTEGER, INTENT(IN OUT) :: PULSEDATE         ! date of pulse start
        INTEGER, INTENT(IN OUT) :: PULSETIME         ! date of pulse end

        !-----------------------------------------------------------------------------
        !.......  Summary of algorithm
        !        1. initalize if first hour
        !        2. no rainfall or new rainfall class less than current one, continue existing NO emission pulse
        !        3. if new rainfall that increases rainfall class, then create new NO emission pulse
        !         using pulsetype, rainrate,ptype, and date/time
        !         if stronger NO pulse compared to previous time step , then
        !         start a new NO emission pulse

        IF( INITIAL_HOUR ) THEN
            PTYPE      = 0
            PULSEDATE  = 0
            PULSETIME  = 0
            PRECIP_ADJ = 1.
        ELSE

            IF( PULSETYPE( RAIN ) <= PTYPE ) THEN

                !.......  No new rainfall or new rainfall class less than current one
                PRECIP_ADJ = PRECIPFACT( PTYPE, JDATE, JTIME, PULSEDATE, PULSETIME )

            ELSE

                !.......  Rainfall class type increases (NO emission pulse generated)
                PULSEDATE = JDATE
                PULSETIME = JTIME
                PTYPE     = PULSETYPE( RAIN )

                PRECIP_ADJ = PRECIPFACT( PTYPE, JDATE, JTIME, PULSEDATE, PULSETIME )
            END IF


        END IF

        RETURN

    END FUNCTION PRECIP_ADJ

    !-----------------------------------------------------------------------------
    !.......  This internal function computes a fertilizer adjustment factor
    !         for the given date in yyyyddd format. If it is not growing
    !         season, the adjustment factor is 0; otherwise, it ranges from
    !         0.0 to 1.0.
    REAL FUNCTION FERTILIZER_ADJ( DATE )

        IMPLICIT NONE

        !.......  Function arguments
        INTEGER, INTENT(IN) :: DATE

        !.......  External functions
        !         INTEGER, EXTERNAL :: GROWSEASON

        !.......  Local variables
        INTEGER  GDAY

        !-----------------------------------------------------------------------------

        GDAY = GROWSEASON( DATE )

        IF( GDAY == 0 ) THEN
            FERTILIZER_ADJ = 0.
        ELSE IF( GDAY >= 1 .AND. GDAY < 30 ) THEN           ! first month of growing season
            FERTILIZER_ADJ = 1.
        ELSE IF( GDAY >= 30 ) THEN
            FERTILIZER_ADJ = 1. + 30. / 184. - FLOAT( GDAY ) / 184.
        ELSE
            WRITE( MESG,94010 ) 'Invalid date specified; date = ', DATE,    &
                                'growing season day = ', GDAY
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        !******************  FORMAT  STATEMENTS   ******************************

        !.......   Internal buffering formats...... 94xxx

94010   FORMAT( A, F10.2, 1X, A, I3, ',', I3 )


        RETURN

    END FUNCTION FERTILIZER_ADJ

    !-----------------------------------------------------------------------------
    !.......  This internal function computes a vegetation adjustment factor
    !         for the given date in yyyyddd format. The adjustment factor
    !         ranges from 0.5 to 1.0.
    REAL FUNCTION VEG_ADJ( DATE )

        IMPLICIT NONE

        !.......  Function arguments
        INTEGER, INTENT(IN) :: DATE

        !.......  External functions
        !         INTEGER, EXTERNAL :: GROWSEASON

        !.......  Local variables
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
            WRITE( MESG,94010 ) 'Invalid date specified; date = ', DATE, &
                                'growing season day = ', GDAY
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        !******************  FORMAT  STATEMENTS   ******************************

        !.......   Internal buffering formats...... 94xxx

94010   FORMAT( A, F10.2, 1X, A, I3, ',', I3 )

        RETURN

    END FUNCTION VEG_ADJ

    !-----------------------------------------------------------------------------
    !.......  This internal function computes the day of the growing season
    !         corresponding to the given date in yyyyddd format.
    INTEGER FUNCTION GROWSEASON( DATE )

        IMPLICIT NONE

        !.......  Function arguments
        INTEGER, INTENT(IN) :: DATE

        !.......  Local parameters
        INTEGER, PARAMETER :: GSEASON_START = 0401
        INTEGER, PARAMETER :: GSEASON_END   = 1031

        !.......  Local variables
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

        IF( JDAY >= GSJULIAN_START .AND. JDAY <= GSJULIAN_END ) THEN      ! growing season
            GROWSEASON = JDAY - GSJULIAN_START + 1
        ELSE IF( JDAY >= 1 .AND. JDAY <= 366 ) THEN      ! before or after growing season
            GROWSEASON = 0
        ELSE
            WRITE( MESG,94010 ) 'Invalid date specified; date = ', DATE, 'jday = ', JDAY
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010   FORMAT( A, F10.2, 1X, A, I3, ',', I3 )

        RETURN

    END FUNCTION GROWSEASON

    !-----------------------------------------------------------------------------
    !.......  This internal function computes a precipitation adjustment
    !         factor from YL 1995 based on a rain rate. The pulse type is
    !         and integer ranging from 0 to 3 indicating the type of
    !         rainfall rate.
    REAL FUNCTION PRECIPFACT( IPULSETYPE, JDATE, JTIME, SDATE, STIME )

        IMPLICIT NONE

        !.......  Function arguments
        INTEGER, INTENT(IN OUT) :: IPULSETYPE
        INTEGER, INTENT(IN) :: JDATE
        INTEGER, INTENT(IN) :: JTIME
        INTEGER, INTENT(IN) :: SDATE
        INTEGER, INTENT(IN) :: STIME

        !.......  Local variables
        REAL     DAYDIFF

        !-----------------------------------------------------------------------------

        DAYDIFF = FLOAT( SECSDIFF( SDATE, STIME, JDATE, JTIME ) )/ ( 24.0 * 3600.0 )

        SELECT CASE( IPULSETYPE )
          CASE( 0 )
            PRECIPFACT = 1.
          CASE( 1 )
            IF( DAYDIFF < 2. ) THEN
                PRECIPFACT = 11.19 * EXP(-0.805*( DAYDIFF + 1.0) )
            ELSE
                IPULSETYPE = 0
                PRECIPFACT = 1.
            END IF
          CASE( 2 )
            IF( DAYDIFF < 6. ) THEN
                PRECIPFACT = 14.68 * EXP(-0.384*( DAYDIFF + 1.0) )
            ELSE
                IPULSETYPE = 0
                PRECIPFACT = 1.
            END IF
          CASE DEFAULT
            IF( DAYDIFF < 13. ) THEN
                PRECIPFACT = 18.46 * EXP(-0.208*( DAYDIFF + 1.0) )
            ELSE
                IPULSETYPE = 0
                PRECIPFACT = 1.
            END IF
        END SELECT

        RETURN

    END FUNCTION PRECIPFACT

    !-----------------------------------------------------------------------------
    !.......  This internal function computes the pulse type from a rainfall
    !         rate (see YL 1995).
    INTEGER FUNCTION PULSETYPE( RAINRATE )

        IMPLICIT NONE

        !.......  Function arguments
        REAL, INTENT(IN) :: RAINRATE

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

END SUBROUTINE HRNO

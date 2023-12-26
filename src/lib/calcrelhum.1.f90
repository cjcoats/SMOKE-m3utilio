
REAL FUNCTION CALCRELHUM( TEMPK, PRESSURE, MIXRATIO )

!***********************************************************************
!  function body starts at line 69
!
!  DESCRIPTION:
!       Calculates relative humidity as a percentage based on air temperature,
!       pressure, and water vapor mixing ratio. Uses Lowe's approximation
!       to calculate saturation vapor pressure, then calculates saturation
!       water vapor mixing ratio.
!
!  PRECONDITIONS REQUIRED:
!       Temperature in Kelvin
!       Pressure in pascals
!       Water vapor mixing ratio in kg/kg
!
!  SUBROUTINES AND FUNCTIONS CALLED:  none
!
!  REVISION  HISTORY:
!     12/03: Created by C. Seppanen
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

    IMPLICIT NONE

!.........  Function arguments
    REAL, INTENT(IN) :: TEMPK     ! temperature in Kelvin
    REAL, INTENT(IN) :: PRESSURE  ! pressure in pascals
    REAL, INTENT(IN) :: MIXRATIO  ! mixing ratio in kg/kg

!.........  Local parameters
    REAL, PARAMETER :: A0 = 6984.505294         ! constants for Lowe's
    REAL, PARAMETER :: A1 = -188.9039310        ! approximation of
    REAL, PARAMETER :: A2 = 2.133357675         ! saturation vapor
    REAL, PARAMETER :: A3 = -1.288580973E-2     ! pressure
    REAL, PARAMETER :: A4 = 4.393587233E-5
    REAL, PARAMETER :: A5 = -8.023923082E-8
    REAL, PARAMETER :: A6 = 6.136820929E-11

    REAL, PARAMETER :: WVRATIO = 0.622      ! ratio of MW of water vapor to dry air

!.........  Local variables
    REAL    SVP         ! saturation vapor pressure
    REAL    SMR         ! saturation mixing ratio
    REAL    RELHUM      ! relative humidity

    CHARACTER(16) :: PROGNAME = 'CALCRELHUM' ! program name

!***************************************************************
!   begin body of function CALCRELHUM

!.........  Calculate saturation vapor pressure; uses Lowe's (1977)
!           polynomial approximation to the Clausius Clayperon equation
    SVP = A0 +&
    &        TEMPK * ( A1 +&
    &            TEMPK * ( A2 +&
    &                TEMPK * ( A3 +&
    &                    TEMPK * ( A4 +&
    &                        TEMPK * ( A5 +&
    &                            TEMPK * A6 )))))

!.........  Convert saturation vapor pressure from millibars to pascals
    SVP = SVP * 100.0

!.........  Calculate saturation mixing ratio
    SMR = WVRATIO * ( SVP / ( PRESSURE - SVP ))

!.........  Calculate relative humidity %
    RELHUM = ( MIXRATIO / SMR ) * 100

!.........  Make sure relative humidity is not more than 100%
!           or less than 0%
    RELHUM = MIN( RELHUM, 100.0 )
    RELHUM = MAX( RELHUM, 0.0 )

    CALCRELHUM = RELHUM

    RETURN

END FUNCTION CALCRELHUM


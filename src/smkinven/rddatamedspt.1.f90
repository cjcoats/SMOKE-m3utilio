
SUBROUTINE RDDATAMEDSPT( LINE, READDATA, READPOL, NPOLPERLN,&
&                        IYEAR, CORS, BLID, DESC, HT, DM, TK,&
&                        FL, VL, SIC, LAT, LON, HDRFLAG )

!***********************************************************************
!  subroutine body starts at line 156
!
!  DESCRIPTION:
!      This subroutine processes a line from an MEDS format point-source inventory
!      file and returns the inventory data values.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!      Created on 10/2013 by B.H. Baek
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
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
!.........  This module contains the information about the source category
    USE MODINFO, ONLY: TMPNAM

    USE MODSOURC, ONLY: NMEDGRD, CMEDGRD

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   SUBROUTINE ARGUMENTS
    CHARACTER(*),       INTENT  (IN) :: LINE                  ! input line
    CHARACTER(*),       INTENT (OUT) ::&
    &                          READDATA( NPOLPERLN,NPTPPOL3 )  ! array of data values
    CHARACTER(NAMLEN3), INTENT (OUT) :: READPOL( NPOLPERLN )  ! array of pollutant names
    INTEGER,            INTENT(INOUT):: NPOLPERLN             ! no. pollutants per line
    INTEGER,            INTENT (OUT) :: IYEAR                 ! inventory year
    CHARACTER(ORSLEN3), INTENT (OUT) :: CORS                  ! DOE plant ID
    CHARACTER(BLRLEN3), INTENT (OUT) :: BLID                  ! boiler ID
    CHARACTER(40),      INTENT (OUT) :: DESC                  ! plant description
    CHARACTER(16),      INTENT (OUT) :: HT                    ! stack height
    CHARACTER(16),      INTENT (OUT) :: DM                    ! stack diameter
    CHARACTER(16),      INTENT (OUT) :: TK                    ! exit temperature
    CHARACTER(16),      INTENT (OUT) :: FL                    ! flow rate
    CHARACTER(9),       INTENT (OUT) :: VL                    ! exit velocity
    CHARACTER(SICLEN3), INTENT (OUT) :: SIC                   ! SIC
    CHARACTER(16),      INTENT (OUT) :: LAT                   ! stack latitude
    CHARACTER(16),      INTENT (OUT) :: LON                   ! stack longitude
    LOGICAL,            INTENT (OUT) :: HDRFLAG               ! true: line is a header line

!...........   Local parameters, indpendent
    INTEGER, PARAMETER :: MXPOLFIL = 53  ! maximum pollutants in file

!...........   Local parameter arrays...

!...........   Other local variables
    INTEGER         I,J,N     ! counters and indices

    INTEGER         ROW, COL ! tmp grid row and col index
    INTEGER         INY      !  inventory year
    INTEGER         IOS      !  i/o status

    CHARACTER(300)      MESG                 !  message buffer
    CHARACTER(CHRLEN3)  ROWCOL           !  Row/Col & GAI lookup variables

    CHARACTER(16) :: PROGNAME = 'RDDATAMEDSPT' ! Program name

!***********************************************************************
!   begin body of subroutine RDDATAMEDSPT

    HDRFLAG = .FALSE.
    NPOLPERLN = 6

    READPOL( 1 ) = 'CO'
    READPOL( 2 ) = 'NOX'
    READPOL( 3 ) = 'SOX'
    READPOL( 4 ) = 'TOG'
    READPOL( 5 ) = 'PM'
    READPOL( 6 ) = 'NH3'

!.........  set year
    INY = STR2INT( LINE( 59:60 ) )
    IYEAR = 2000 + INY

!.........  Read source data
    CORS = ''  ! DOE plant ID
    BLID = ''  ! boiler ID
    DESC = ''  ! plant description

    HT   = ADJUSTL( LINE( 74:78 ) ) ! stack height
    DM   = '0.1'         ! stack diameter
    TK   = '273.15'      ! exit temperature
    FL   = '0.1'         ! flow rate
    VL   = '0.1'         ! exit velocity
    SIC  = ADJUSTL( LINE( 23:36 ) ) ! SIC

    COL  = STR2INT( LINE( 37:39 ) )
    ROW  = STR2INT( LINE( 40:42 ) )
    WRITE( ROWCOL,'( 2I3.3 )' ) COL, ROW
    N = INDEX1( ROWCOL, NMEDGRD, CMEDGRD( :,1 ) )

    LON  = ADJUSTL( CMEDGRD( N,2 ) )  ! grid cell longitude
    LAT  = ADJUSTL( CMEDGRD( N,3 ) )  ! grid cell latitude

!.........  Set all MEDS annual/avg to zero since all are daily/hourly inv
    READDATA = '0.0'

    RETURN

END SUBROUTINE RDDATAMEDSPT

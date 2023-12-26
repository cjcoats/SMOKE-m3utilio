
INTEGER FUNCTION CVTVEHTYPE( SMKVEH )

!***********************************************************************
!  subroutine body starts at line 53
!
!  DESCRIPTION:
!       Converts inventory vehicle type to MOBILE6 type
!
!  PRECONDITIONS REQUIRED: none
!
!  SUBROUTINES AND FUNCTIONS CALLED: none
!
!  REVISION  HISTORY:
!     10/01: Created by C. Seppanen
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
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

!.........  Function arguments
    INTEGER, INTENT (IN) :: SMKVEH   ! vehicle type in SMOKE code

!...........   LOCAL VARIABLES and their descriptions:
    LOGICAL :: EFLAG      = .FALSE.   ! true: error found

    CHARACTER(300)          MESG      !  message buffer

    CHARACTER(16) :: PROGNAME = 'CVTVEHTYPE' ! program name

!***********************************************************************
!   begin body of function CVTVEHTYPE

    SELECT CASE( SMKVEH )
      CASE( 100 )              ! light duty gasoline vehicles
        CVTVEHTYPE = 1
      CASE( 102 )              ! light duty gasoline trucks 1
        CVTVEHTYPE = 2
      CASE( 104 )              ! light duty gasoline trucks 2
        CVTVEHTYPE = 3
      CASE( 107 )              ! heavy duty gasoline vehicles
        CVTVEHTYPE = 4
      CASE( 3000 )             ! light duty diesel vehicles
        CVTVEHTYPE = 5
      CASE( 3006 )             ! light duty diesel trucks
        CVTVEHTYPE = 6
      CASE( 3007 )             ! heavy duty diesel vehicles
        CVTVEHTYPE = 7
      CASE( 108 )              ! motorcycles
        CVTVEHTYPE = 8
      CASE DEFAULT
        EFLAG = .TRUE.

        WRITE( MESG, '(A,I8,A)' ) 'ERROR: Vehicle type ', SMKVEH,&
        &                     'is not recognized.'
        CALL M3MESG( MESG )
    END SELECT

    IF( EFLAG ) THEN
        MESG = 'Problem converting vehicle type'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    RETURN

END FUNCTION CVTVEHTYPE

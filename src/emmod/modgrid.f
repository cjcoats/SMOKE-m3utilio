
        MODULE MODGRID

!***********************************************************************
!  Module body starts at line
!
!  DESCRIPTION:
!     This module contains the public variables and allocatable arrays
!     used for gridding.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION HISTORY:
!     Created 10/2000 by M. Houyoux
!
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!***************************************************************************
!
! Project Title: EDSS Tools Library
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

        IMPLICIT NONE

!.........  Horizontal grid information
        CHARACTER(NAMLEN3), PUBLIC, SAVE :: GRDNM = ' '  ! grid name
        CHARACTER(NAMLEN3), PUBLIC, SAVE :: COORD = ' '  ! coord system name
        INTEGER, PUBLIC, SAVE :: GDTYP = -1     ! i/o api grid type code
        REAL(8), PUBLIC, SAVE :: P_ALP = 0.D0   ! projection alpha
        REAL(8), PUBLIC, SAVE :: P_BET = 0.D0   ! projection beta
        REAL(8), PUBLIC, SAVE :: P_GAM = 0.D0   ! projection gamma
        REAL(8), PUBLIC, SAVE :: XCENT = 0.D0   ! x-center of projection
        REAL(8), PUBLIC, SAVE :: YCENT = 0.D0   ! y-center of projection
        REAL(8), PUBLIC, SAVE :: XORIG = 0.D0   ! x-origin of grid
        REAL(8), PUBLIC, SAVE :: YORIG = 0.D0   ! y-origin of grid
        REAL(8), PUBLIC, SAVE :: XCELL = 0.D0   ! x-dim of cells
        REAL(8), PUBLIC, SAVE :: YCELL = 0.D0   ! y-dim of cells
        INTEGER, PUBLIC, SAVE :: NCOLS = 0      ! number of columns in grid
        INTEGER, PUBLIC, SAVE :: NROWS = 0      ! number of rows in grid
        INTEGER, PUBLIC, SAVE :: NGRID = 0      ! number of cells in grid
        INTEGER, PUBLIC, SAVE :: XDIFF = 0      ! subgrid fewer cols (nxsub = nx - xdiff)
        INTEGER, PUBLIC, SAVE :: YDIFF = 0      ! subgrid fewer cols (nysub = ny - ydiff)
        INTEGER, PUBLIC, SAVE :: XOFF  = 0      ! subgrid offset (x-sub = x - xoff)
        INTEGER, PUBLIC, SAVE :: YOFF  = 0      ! subgrid offset
        INTEGER, PUBLIC, SAVE :: XOFF_A= 0      ! tmp subgrid offset (x-sub = x - xoff)
        INTEGER, PUBLIC, SAVE :: YOFF_A= 0      ! tmp subgrid offset
        LOGICAL, PUBLIC, SAVE :: OFFLAG = .FALSE. ! true: subgrid offset has been set

!.........  Vertical structure information
        INTEGER, PUBLIC, SAVE :: NLAYS =  1      ! number of layers
        INTEGER, PUBLIC, SAVE :: VGTYP = -1      ! type of vertical coordinates
        REAL   , PUBLIC, SAVE :: VGTOP = 0.0     ! model-top, for sigma coord types
        REAL   , ALLOCATABLE, PUBLIC, SAVE :: VGLVS( : ) ! vertical coordinate values

        END MODULE MODGRID

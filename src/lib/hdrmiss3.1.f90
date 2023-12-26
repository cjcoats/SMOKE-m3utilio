
SUBROUTINE HDRMISS3

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!     This subroutine initializes the I/O API header to missing and default
!     values appropriate for emission files.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 3/99 by M Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!***********************************************************************
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

!***********************************************************************
!   begin body of subroutine HDRMISS3

    FTYPE3D = GRDDED3
    SDATE3D = 0
    STIME3D = 0
    TSTEP3D = 0
    NCOLS3D = 1
    NROWS3D = 1
    NLAYS3D = 1
    NTHIK3D = 1
    NVARS3D = 0
    GDTYP3D = IMISS3
    P_ALP3D = BADVAL3
    P_BET3D = BADVAL3
    P_GAM3D = BADVAL3
    XCENT3D = BADVAL3
    YCENT3D = BADVAL3
    XORIG3D = BADVAL3
    YORIG3D = BADVAL3
    XCELL3D = BADVAL3
    YCELL3D = BADVAL3
    VGTYP3D = IMISS3
    VGTOP3D = BADVAL3
    VGLVS3D = 0.      ! array
    GDNAM3D = ' '

    VNAME3D = ' '     ! array
    VTYPE3D = 0       ! array
    UNITS3D = ' '     ! array
    VDESC3D = ' '     ! array
    FDESC3D = ' '     ! array

    RETURN

END SUBROUTINE HDRMISS3

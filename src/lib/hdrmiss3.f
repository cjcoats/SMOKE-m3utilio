
        SUBROUTINE HDRMISS3

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C     This subroutine initializes the I/O API header to missing and default
C     values appropriate for emission files.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 3/99 by M Houyoux
C       Version 11/2023 by CJC:  USE M3UTILIO and related changes
C
C***********************************************************************
C
C Project Title: EDSS Tools Library
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$
C
C****************************************************************************

        USE M3UTILIO

        IMPLICIT NONE

C***********************************************************************
C   begin body of subroutine HDRMISS3

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

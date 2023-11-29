
        SUBROUTINE CONVRTXY( NSRC, CTYPE, GDNAM, P_ALP, P_BET,
     &                       P_GAM, XCENT, YCENT, XVALS, YVALS )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine converts coordinates from lat-lon to the grid
C      defined by the subroutine arguments.
C
C  PRECONDITIONS REQUIRED:
C
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Added POLGRD3 as a supported coord sys type - E. Giroux CNRC 03/2004
!
!       Version 10/2016 by C. Coats:  USE M3UTILIO
C
C       Version 11/2023 by C. Coats:  USE MODGCTP's XY2XY
C
C****************************************************************************/
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
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
C***************************************************************************
        USE M3UTILIO
        USE MODGCTP

        IMPLICIT NONE

C.........  SUBROUTINE ARGUMENTS
        INTEGER,      INTENT(IN   ) :: NSRC          !  actual source count
        INTEGER,      INTENT(IN   ) :: CTYPE         !  coord sys type
        CHARACTER(*), INTENT(IN   ) :: GDNAM         !  grid name
        REAL(8),      INTENT(IN   ) :: P_ALP         !  first, second, third map
        REAL(8),      INTENT(IN   ) :: P_BET         !  projection descriptive
        REAL(8),      INTENT(IN   ) :: P_GAM         !  parameters
        REAL(8),      INTENT(IN   ) :: XCENT         !  lon for coord-system X=0
        REAL(8),      INTENT(IN   ) :: YCENT         !  lat for coord-system Y=0
        REAL(8),      INTENT(INOUT) :: XVALS( NSRC ) !  x location (input lon)
        REAL(8),      INTENT(INOUT) :: YVALS( NSRC ) !  y location (input lat)

C...........   Other local variables

        REAL(8)     XOUTS( NSRC ), YOUTS( NSRC )

C***********************************************************************
C   begin body of subroutine CONVRTXY

        CALL XY2XY( LATGRD3, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,
     &              CTYPE,   P_ALP, P_BET, P_GAM, XCENT, YCENT,
     &              NSRC, XVALS, YVALS, XOUTS, YOUTS )

        XVALS = XOUTS     !  since CONVRTXY does transformation in-place
        YVALS = YOUTS
        
        RETURN

        END SUBROUTINE CONVRTXY


        SUBROUTINE WRDAYMSG( JDATE, MESG )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      Writes a text message to stdout and log file about which day is
C      being processed
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 9/99 by M. Houyoux
C       Version 11/2023 by CJC:  USE M3UTILIO and related changes
C
C**************************************************************************
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

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: JDATE    ! Julian date
        CHARACTER(*), INTENT(OUT) :: MESG     ! message buffer

C...........   Local variables
        INTEGER         DAY          !  day of week number

C***********************************************************************
C   begin body of subroutine WRDAYMSG

        DAY = WKDAY( JDATE )

        MESG= 'Processing '// TRIM( DAYS( DAY ) )// ' '// MMDDYY( JDATE )
        CALL M3MSG2( MESG )

        RETURN

        END SUBROUTINE WRDAYMSG


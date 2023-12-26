
SUBROUTINE WRPTREF( NPSRC, IDIU, IWEK, IMON )

!***********************************************************************
!  subroutine body starts at line 78
!
!  DESCRIPTION:
!      This subroutine writes a temporal x-ref file in SMOKE format. It is
!      only used for the EMS-95 input files.
!
!  PRECONDITIONS REQUIRED:
!      Input arrays populated
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!      Subroutines: I/O API subroutine
!
!  REVISION  HISTORY:
!      Created 10/98 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!****************************************************************************/
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

!.........  MODULES for public variables
!.........  This module contains the information about the source category
    USE MODINFO, ONLY: CRL

    IMPLICIT NONE

!...........   INCLUDES

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!.........  SUBROUTINE ARGUMENTS
    INTEGER, INTENT (IN) :: NPSRC          !  actual source count
    INTEGER, INTENT (IN) :: IDIU( NPSRC )  !  source FIPS (county) ID
    INTEGER, INTENT (IN) :: IWEK( NPSRC )  !  source SCC
    INTEGER, INTENT (IN) :: IMON( NPSRC )  !  source SIC

!...........   Other local variables
    INTEGER         S

    INTEGER         FDEV             !  output file unit number

    CHARACTER(300)  MESG             !  message buffer

    CHARACTER(16) :: PROGNAME = 'WRPTREF' !  program name

!***********************************************************************
!   begin body of program WRPTREF

    MESG = 'Enter the name of the TEMPORAL X-REF output file'
    FDEV = PROMPTFFILE( MESG, .FALSE., .TRUE.,&
    &                    CRL // 'TREF_ALT', PROGNAME )

    MESG = 'Writing out TEMPORAL CROSS-REFERENCE file...'
    CALL M3MSG2( MESG )

    DO S = 1, NPSRC

        WRITE( FDEV, 93040, ERR=6001 ) 1, IWEK( S ), IDIU( S )

    ENDDO

    RETURN

6001 MESG = 'ERROR writing temporal x-ref file'
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93040 FORMAT( I5, ',' ,I5, ',' ,I5 )

END SUBROUTINE WRPTREF

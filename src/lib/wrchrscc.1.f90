
SUBROUTINE WRCHRSCC( CSCC )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine writes character string SCC codes
!
!  PRECONDITIONS REQUIRED:
!      Index sorted in order of increasing ISCC value and ISCC populated
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!      Subroutines: I/O API subroutines
!
!  REVISION  HISTORY:
!      Created 10/98 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!************************************************************************
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
!.........  This module contains the lists of unique source characteristics
    USE MODLISTS, ONLY: NINVSCC, INVSCC

!.........  This module contains the information about the source category
    USE MODINFO, ONLY: CRL, NSRC

    IMPLICIT NONE

!...........   INCLUDES

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!.........  SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: CSCC( NSRC )  !  unsorted SCCs

!...........   Other local variables
    INTEGER       J, S, JSTAT
    INTEGER       FDEV              !  output file unit number

    CHARACTER(300)     MESG, TEXT   !  message buffer
    CHARACTER(SCCLEN3) TSCC, LSCC   !  current and previous 10-digit SCC

    CHARACTER(16), PARAMETER :: PNAME = 'WRCHRSCC' !  program name

!***********************************************************************
!   begin body of subroutine WRCHRSCC

!.............  Get unique-SCC output file
    FDEV = PROMPTFFILE(&
    &        'Enter the name of the ACTUAL SCC output file',&
    &        .FALSE., .TRUE., CRL // 'SCC', PNAME )

!.........  Write out SCCs list (this subroutine expect the data structure
!           that is being provided.

    CALL M3MSG2( 'Writing out ACTUAL SCC file...' )

!.........  Call the routine that generates the unique SCC list
    CALL GENUSLST

!.........  Write unique SCCs list
    DO J = 1, NINVSCC

        WRITE( FDEV, '(A)', ERR=6001, IOSTAT=JSTAT, IOMSG=TEXT ) INVSCC( J )

    END DO

    RETURN

6001 WRITE( MESG, '(A,I10)' ) 'Problem writing SCC file.  IOSTAT=', JSTAT
    CALL PERROR( TEXT )
    CALL M3EXIT( PNAME, 0, 0, MESG, 2 )

END SUBROUTINE WRCHRSCC


SUBROUTINE RDSMAT( FNAME, VDESC, SMAT )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine reads a speciation matrix for any source category.
!      The variable description is provided as the name for reading, and
!      this routine finds which variable name to read and reads in the data.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created ??/???? by ???
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
!.........  This module is required by the FileSetAPI
    USE MODFILESET

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
    INCLUDE 'SETDECL.EXT'   !  FileSetAPI function declarations

!.........  SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: FNAME     ! speciation matrix file name
    CHARACTER(*), INTENT (IN) :: VDESC     ! variable description to read
    REAL        , INTENT(OUT) :: SMAT( * ) ! coeffs for sources

!.........  Other local variables
    INTEGER         J, L, L1, L2, V       !  counters and indices

    CHARACTER(MXDLEN3) DBUF    !  variable description buffer
    CHARACTER(NAMLEN3) VBUF    !  variable name buffer
    CHARACTER(300)     MESG    !  message buffer

    CHARACTER(16) :: PROGNAME = 'RDSMAT' ! program name

!***********************************************************************
!   begin body of subroutine RDSMAT

!.........  Retrieve file header
    IF ( .NOT. DESCSET( FNAME, ALLFILES ) ) THEN
        MESG = 'Could not get description of file ' // FNAME
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Make sure variable name has proper spacing that is the same as
!           file header
    L1 = INDEX( VDESC, SPJOIN )
    L2 = LEN_TRIM( VDESC )
    DBUF = VDESC( 1:L1-1 )
    DBUF = DBUF( 1:NAMLEN3 ) // SPJOIN // VDESC( L1+1:L2 )

!.........  Find variable description in list of descriptions
    J = INDEX1( DBUF, NVARSET, VDESCSET )

    IF( J .LE. 0 ) THEN
        MESG = 'INTERNAL ERROR: Speciation variable description "'//&
        &       DBUF // '" not found in speciation matrix.'
        CALL M3MSG2( MESG )
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
    END IF

!.........  Read variable and print nice error message if cannot
    VBUF = VNAMESET( J )
    IF ( .NOT. READSET( FNAME, VBUF, 1, ALLFILES,&
    &                    0, 0, SMAT ) ) THEN

        L  = LEN_TRIM( FNAME )
        L1 = LEN_TRIM( VBUF )
        MESG = 'Could not read speciation matrix from file "' //&
        &       FNAME( 1:L ) // '"' // CRLF() // BLANK10 //&
        &       'for variable "' // VBUF( 1:L1 ) //&
        &       '" with description "' // VDESC( 1:L2 ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF    !  if readset() failed for speciation matrix

    RETURN

END SUBROUTINE RDSMAT

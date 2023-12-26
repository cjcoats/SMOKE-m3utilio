
SUBROUTINE RDGMAT( FNAME, NGRID, NMAT1, NMAT2, NX, IX, CX )

    !***************************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine reads a gridding matrix for any source category
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version ??/???? by ???
    !       Version 11/2023 by CJC:  USE M3UTILIO, ".f90" source format, and
    !       related changes
    !***************************************************************************
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

    IMPLICIT NONE

    !.......  SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: FNAME           ! gridding matrix name
    INTEGER     , INTENT (IN) :: NGRID           ! number of grid cells
    INTEGER     , INTENT (IN) :: NMAT1           ! dim 1 for matrix
    INTEGER     , INTENT (IN) :: NMAT2           ! dim 2 for matrix
    INTEGER     , INTENT(OUT) :: NX( NGRID )     ! number of sources per cell
    INTEGER                   :: IX( NMAT1 )     ! list of sources per cell
    REAL                      :: CX( NMAT2 )     ! coefficients for sources

    !.......  Other local variables
    INTEGER         C           !  tmp cell number
    INTEGER         NSUM        !  count of gridding matrix size

    CHARACTER(300)  MESG        !  message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'RDGMAT'     ! program name

    !***********************************************************************
    !   begin body of subroutine RDGMAT

    !.......  Read matrix
    IF ( .NOT. READ3( FNAME, 'ALL', 1, 0, 0, NX ) ) THEN

        MESG = 'Could not read gridding matrix from file "' // TRIM( FNAME ) // '".'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF          !  if read3() failed for gridding matrix

    !.......  Check to make sure that data are consistent with header
    NSUM = 0
    DO C = 1, NGRID
        NSUM = NSUM + NX( C )
    ENDDO

    IF( NSUM .GT. NMAT1 ) THEN

        MESG = 'ERROR: Gridding matrix dimension is inconsistent '//    &
               'with records count    !' // CRLF() // '          ' //   &
               'Delete gridding matrix and recreate it.'

        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    RETURN

END SUBROUTINE RDGMAT

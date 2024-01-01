
SUBROUTINE RDUMAT( FNAME, NSRC, NMAT1, NMAT2, NU, IU, CU )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine reads an ungridding matrix for any source category
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version ??/???? by ???
    !       Version 11/2023 by CJC:  USE M3UTILIO, ".f90" source format, and 
    !       related changes
    !****************************************************************************/
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !         System
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
    CHARACTER(*), INTENT (IN) :: FNAME           ! ungridding matrix name
    INTEGER     , INTENT (IN) :: NSRC            ! number of sources
    INTEGER     , INTENT (IN) :: NMAT1           ! dim 1 for matrix
    INTEGER     , INTENT (IN) :: NMAT2           ! dim 2 for matrix
    INTEGER     , INTENT(OUT) :: NU( NSRC  )     ! number of cells per source
    INTEGER                   :: IU( NMAT1 )     ! list of cells per source
    REAL                      :: CU( NMAT2 )     ! coefficients for cells

    !.......  Other local variables
    CHARACTER(300)  MESG        !  message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'RDUMAT'     ! program name

    !***********************************************************************
    !   begin body of subroutine RDUMAT

    MESG = 'Reading ungridding matrix...'
    CALL M3MSG2( MESG )

    IF ( .NOT. READ3( FNAME, 'ALL', 1, 0, 0, NU ) ) THEN

        MESG = 'Could not read ungridding matrix from file "' // TRIM( FNAME ) // '".'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF          !  if read3() failed for ungridding matrix

    RETURN

END SUBROUTINE RDUMAT

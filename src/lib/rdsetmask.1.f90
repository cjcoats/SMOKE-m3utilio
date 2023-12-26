
SUBROUTINE RDSETMASK( FNAME, JDATE, JTIME, NDIM, NDIM2, NVLIST,&
&                      VNAMES, VINDX, OUTVAR )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine reads in certain variables from a list of variables,
!      depending on the index, which indicates if those variables are expected
!      to be present in a file set.  The index is a pointer which says which
!      part of the output array to read the variable into. The variables are
!      assumed to be reals.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Version ??/???? by ???
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!****************************************************************************/
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
!***************************************************************************
    USE M3UTILIO

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'SETDECL.EXT'   !  FileSetAPI function declarations

!.........  SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: FNAME     ! i/o api file name
    INTEGER     , INTENT (IN) :: JDATE     ! Julian date (YYYYDDD)
    INTEGER     , INTENT (IN) :: JTIME     ! time (HHMMSS)
    INTEGER     , INTENT (IN) :: NDIM      ! dimension for output array
    INTEGER     , INTENT (IN) :: NDIM2     ! 2nd dimension for output array
    INTEGER     , INTENT (IN) :: NVLIST    ! variable description to read
    CHARACTER(*), INTENT (IN) :: VNAMES( NVLIST ) ! variable names
    INTEGER     , INTENT (IN) :: VINDX ( NVLIST ) ! var index to OUTVAR
    REAL        , INTENT(OUT) :: OUTVAR( NDIM,NDIM2 ) ! coeffs for sources

!.........  Other local variables
    INTEGER         J, L, L1, V       !  counters and indices

    CHARACTER(NAMLEN3) VBUF    !  variable name buffer
    CHARACTER(300)     MESG    !  message buffer

    CHARACTER(16) :: PROGNAME = 'RDSETMASK' ! program name

!***********************************************************************
!   begin body of subroutine RDSETMASK

!.........  Loop through variables that are possibilities for reading
    DO V = 1, NVLIST

!.............  Check variable index to see if this variable is to be read
        J = VINDX( V )
        IF( J .EQ. 0 ) CYCLE  ! to go bottom of loop

!.............  Bounds check
        IF( J .GT. NDIM2 ) THEN

            MESG = 'INTERNAL ERROR: Prevented overflow for read '//&
            &       'of variable "' // TRIM( VBUF ) // '"'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        END IF

!.............  Read variable and print nice error message if cannot
        VBUF = VNAMES( V )
        IF ( .NOT. READSET( FNAME, VBUF, 1, ALLFILES, JDATE,&
        &                    JTIME, OUTVAR( 1,J ) ) ) THEN

            L  = LEN_TRIM( FNAME )
            L1 = LEN_TRIM( VBUF )
            MESG = 'Could not read variable "' // VBUF( 1:L1 ) //&
            &       '" from file "' // FNAME( 1:L ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF    !  if readset() failed for file

    END DO        !  End loop on possible variables

    RETURN

END SUBROUTINE RDSETMASK

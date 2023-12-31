
        SUBROUTINE RD3MASK( FNAME, JDATE, JTIME, NDIM, NDIM2, NVLIST,
     &                      VNAMES, VINDX, OUTVAR )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine reads in certain variables from a list of variables,
C      depending on the index, which indicates if those variables are expected
C      to be present in a file.  The index is a pointer which says which
C      part of the output array to read the variable into. The variables are
C      assumed to be reals.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Version ??/???? by ???
C       Version 11/2023 by CJC:  USE M3UTILIO and related changes
C****************************************************************************/
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
C***************************************************************************
        USE M3UTILIO

        IMPLICIT NONE

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: FNAME     ! i/o api file name
        INTEGER     , INTENT (IN) :: JDATE     ! Julian date (YYYYDDD)
        INTEGER     , INTENT (IN) :: JTIME     ! time (HHMMSS)
        INTEGER     , INTENT (IN) :: NDIM      ! dimension for output array
        INTEGER     , INTENT (IN) :: NDIM2     ! 2nd dimension for output array
        INTEGER     , INTENT (IN) :: NVLIST    ! variable description to read
        CHARACTER(*), INTENT (IN) :: VNAMES( NVLIST ) ! variable names
        INTEGER     , INTENT (IN) :: VINDX ( NVLIST ) ! var index to OUTVAR
        REAL        , INTENT(OUT) :: OUTVAR( NDIM,NDIM2 ) ! coeffs for sources

C.........  Other local variables
        INTEGER         J, L, L1, V       !  counters and indices

        CHARACTER(NAMLEN3) VBUF    !  variable name buffer
        CHARACTER(300)     MESG    !  message buffer

        CHARACTER(16) :: PROGNAME = 'RD3MASK' ! program name

C***********************************************************************
C   begin body of subroutine RD3MASK

C.........  Loop through variables that are possibilities for reading
        DO V = 1, NVLIST

C.............  Check variable index to see if this variable is to be read
            J = VINDX( V )
            IF( J .EQ. 0 ) CYCLE  ! to go bottom of loop

            VBUF = VNAMES( V )

C.............  Bounds check
            IF( J .GT. NDIM2 ) THEN

                MESG = 'INTERNAL ERROR: Prevented overflow for read '//
     &                 'of variable "' // TRIM( VBUF ) // '"'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

            END IF

C.............  Read variable and print nice error message if cannot
            IF ( .NOT. READ3(
     &           FNAME, VBUF, 1, JDATE, JTIME, OUTVAR( 1,J ) ) ) THEN

                MESG = 'Could not read variable "' // TRIM( VBUF ) //
     &                 '" from file "' // TRIM( FNAME ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF    !  if read3() failed for file

        END DO        !  End loop on possible variables

        RETURN

        END SUBROUTINE RD3MASK

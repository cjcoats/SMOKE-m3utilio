
SUBROUTINE GENGGMAT( GNAME, FDEV, MXSCEL, NASRC, NMATX, &
                     NX, IX, CX, NCOEF, CMAX, CMIN )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created  2/0209 by B. Baek
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !************************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !                System
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
    !************************************************************************
    USE M3UTILIO

    !.......   MODULES for public variables
    !.......   This module is the source inventory arrays
    USE MODSOURC, ONLY: XLOCA, YLOCA, CIFIP, CELLID, CSOURC

    !.......  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: NGRID, NCOLS, NROWS

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: NCHARS, NSRC

    IMPLICIT NONE

    !.......   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: GNAME             ! gridding mtx logical name
    INTEGER     , INTENT (IN) :: FDEV              ! surg codes report file
    INTEGER     , INTENT (IN) :: MXSCEL            ! max sources per cell
    INTEGER     , INTENT (IN) :: NASRC             ! no. mobile sources
    INTEGER     , INTENT (IN) :: NMATX             ! no. source-cell intersects
    INTEGER     , INTENT(OUT) :: NX  ( NGRID )     ! no. srcs per cell
    INTEGER     , INTENT(OUT) :: IX  ( NMATX )     ! src IDs
    REAL        , INTENT(OUT) :: CX  ( NMATX )     ! gridding coefficients
    INTEGER     , INTENT(OUT) :: NCOEF             ! no. of gridding coeffs
    INTEGER     , INTENT(OUT) :: CMAX              ! max no. of sources per cell
    INTEGER     , INTENT(OUT) :: CMIN              ! min no. of sources per cell

    !.......   EXTERNAL FUNCTIONS
    LOGICAL, EXTERNAL :: DSCM3GRD

    !.......   Local allocatable arrays...

    !.......   Scratch Gridding Matrix (subscripted by source-within-cell, cell)
    INTEGER     IS ( MXSCEL, NGRID )     ! source IDs for each cell
    REAL        CS ( MXSCEL, NGRID )     ! factors

    !.......   Temporary array for flagging sources that are outside the
    !              domain and for flagging sources with all zero surrogate fractions
    LOGICAL     INDOMAIN( NASRC )     ! true: src is in the domain

    !.......   Temporary arrays for storing surrogate codes to use
    INTEGER     SURGID1( NSRC )     ! primary surrogate code
    INTEGER     SURGID2( NSRC )     ! secondary surrogate code

    !.......   Other local variables
    INTEGER         C, F, I, J, K, N, S     !  indices and counters.

    INTEGER         COL         ! tmp column
    INTEGER         ID1,ID2     ! tmp primary and secondary surg codes
    INTEGER         IOS         ! i/o status
    INTEGER         ISIDX       ! tmp surrogate ID code index
    INTEGER         JMAX        ! counter for storing correct max dimensions
    INTEGER         L2          ! string length
    INTEGER         NCEL        ! tmp number of cells
    INTEGER         NNOSRG      ! no. of cy/st/co codes with no surrogates
    INTEGER         ROW         ! tmp row

    REAL            FRAC        ! tmp surrogate fraction

    LOGICAL      :: EFLAG = .FALSE.      !  true: error detected
    LOGICAL      :: LFLAG = .FALSE.      !  true: location data available
    LOGICAL      :: XYSET = .FALSE.     ! true: X/Y available for src

    CHARACTER(FIPLEN3) CFIP       !  tmp country/state/county code
    CHARACTER(FIPLEN3) LFIP       !  cy/st/co code from previous iteration
    CHARACTER(16)   COORUNIT      !  coordinate system projection units
    CHARACTER(80)   GDESC         !  grid description
    CHARACTER(256)  MESG          !  message buffer

    CHARACTER(SRCLEN3)    CSRC      ! tmp source chars string

    CHARACTER(16), PARAMETER :: PROGNAME = 'GENGGMAT'     ! program name

    !***********************************************************************
    !   begin body of subroutine GENGGMAT

    !.......  Initialize number of sources per cell counter
    NX = 0       ! Array

    SURGID1 = 0       ! array
    SURGID2 = 0       ! array

    !.......  Set flag to indicate that XLOCA/YLOCA are available
    LFLAG = ALLOCATED( XLOCA )

    !.......   Compute gridding matrix:
    !.......       First case:   explicit link (ILINK > 0
    !.......       Second case:  some LNKDEF entry applies
    !.......       Third case:   FIP+roadtype cross-reference match
    !.......       Fourth case:  state+roadtype cross-reference match
    !.......       fifth case:   roadtype cross-reference match
    !.......       sixth case:   fallback default

    MESG = 'Computing gridding matrix and statistics...'
    CALL M3MSG2( MESG )

    LFIP  = ' '
    NNOSRG   = 0
    JMAX  = -1

    DO S = 1, NSRC

        CFIP = CIFIP ( S )
        C    = CELLID( S )
        CSRC = CSOURC( S )

        !.......  Initialize sources as being in the domain
        INDOMAIN( S ) = .TRUE.

        !.......  Special case for source already assigned a grid cell
        IF( C .GT. 0 ) THEN

            J = NX( C )

            !.......  Check that the maximum number of sources per cell is ok
            IF( J .LT. MXSCEL ) THEN
                J = J + 1
                IS( 1,C ) = S
                CS( 1,C ) = 1.0

            ELSE
                IF( J+1 .GT. JMAX ) JMAX = J+1
            END IF

            NX( C ) = J

            CYCLE               ! To head of loop over sources
        END IF

    END DO            !  end loop on sources S, computing gridding matrix.

    !.......  Abort if overflow occurred
    IF ( JMAX .GT. MXSCEL ) THEN

        WRITE( MESG,94010 )                                     &
         'INTERNAL ERROR: Gridding matrix not written.' //      &
         CRLF() // BLANK10 // 'Arrays would have overflowed.'// &
         CRLF() // BLANK10 // 'Current maximum sources per cell (MXSCEL) =', MXSCEL, '.' // &
         CRLF() // BLANK10 // 'Actual  maximum sources per cell          =', JMAX  , '.'
        CALL M3MSG2( MESG )
        CALL M3EXIT( PROGNAME, 0, 0, 'Failure', 2 )

    END IF

    !.......  Compress matrix into I/O representation from scratch
    !.......  representation and compute statistics.

    K    = 0
    CMAX = 0
    CMIN = 99999999
    DO C = 1, NGRID      ! Loop through cells

        J = NX( C )

        IF (      J .GT. CMAX ) THEN
            CMAX = J
        ELSE IF ( J .GT. 0 .AND. J .LT. CMIN ) THEN
            CMIN = J
        END IF

        DO N = 1, J      ! Loop through sources in this cell
            K = K + 1
            IF ( K .LE. NMATX ) THEN
                S       = IS( N,C )
                IX( K ) = S
                CX( K ) = CS( N,C )
            END IF
        END DO

    END DO        !  end of loop on cells C for this FIP

    NCOEF = K

    !.......  Write gridding matrix
    MESG = 'Writing out GRIDDING MATRIX file...'
    CALL M3MSG2( MESG )

    IF( .NOT. WRITE3( GNAME, 'ALL', 0, 0, NX ) ) THEN
        MESG = 'Error writing GRIDDING MATRIX file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Write output surrogates codes
    DO S = 1, NSRC
        WRITE( FDEV,93360 ) S, SURGID1( S ), SURGID2( S )
    END DO

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Formatted file I/O formats...... 93xxx

93360 FORMAT( I8, 1X, I8, 1X, I8 )

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I9, :, 1X ) )

END SUBROUTINE GENGGMAT


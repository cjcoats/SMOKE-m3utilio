
SUBROUTINE GENPGMAT( FNAME, NPSRC, NGRID, XLOCA, YLOCA, VFLAG,&
                     NX, IX, NCOEF, CMAX, CMIN )

    !***********************************************************************
    !  subroutine body starts at line 102
    !
    !  DESCRIPTION:
    !      This subroutine creates the point source gridding matrix and
    !      writes it to a I/O API file.  This subroutine is created to be
    !      consisent with those for area and mobile sources, which are needed
    !      to ensure that the sparse matrix is stored contiguously in memory.
    !
    !  PRECONDITIONS REQUIRED:
    !      File must be opened and its logical name input through the FNAME
    !      argument.  Memory for NX and IX must be allocated prior to the
    !      subroutine call to ensure that it is contiguous. The x and y coordinates
    !      must already be converted to the output grid.
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created by M. Houyoux 1/99
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
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

    IMPLICIT NONE

    !.......   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: FNAME     ! matrix output inventory logical name
    INTEGER     , INTENT (IN) :: NPSRC              ! Actual source count
    INTEGER     , INTENT (IN) :: NGRID              ! Actual grid cell count
    REAL(8)     , INTENT (IN) :: XLOCA( NPSRC )     ! X-coordinate in proper coordinate system
    REAL(8)     , INTENT (IN) :: YLOCA( NPSRC )     ! Y-coordinate in proper coordinate system
    LOGICAL     , INTENT (IN) :: VFLAG              ! true: using variable grid
    INTEGER     , INTENT(OUT) :: NX( NGRID )        ! Number of sources per cell
    INTEGER     , INTENT(OUT) :: IX( NPSRC )        ! Source number, w/ order based on NX
    INTEGER     , INTENT(OUT) :: NCOEF              ! Number of coefficients
    INTEGER     , INTENT(OUT) :: CMAX               ! Max no. of srcs per cell
    INTEGER     , INTENT(OUT) :: CMIN               ! Min no. of srcs per cell

    !.......   EXTERNAL FUNCTIONS and their descriptions:
    LOGICAL, EXTERNAL :: DSCM3GRD

    !.......   Scratch Gridding Matrix

    INTEGER     INDX( NPSRC )
    INTEGER     GN  ( NPSRC )
    INTEGER     SN  ( NPSRC )

    !.......   Other local variables

    INTEGER         C, I, J, K, K2, R     !  indices and counters.

    INTEGER         IOS       ! i/o status
    INTEGER         KSAV      ! saved value of K
    INTEGER         K2SAV     ! saved value of K2
    INTEGER         NEXCLD    ! number of sources not in grid.  Will
        !    appear as blanks in SN( ) and GN( )

    LOGICAL         GFLAG         !  generate output gridding matrix if true

    CHARACTER(300)  MESG          !  message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'GENPGMAT'     ! program name

    !***********************************************************************
    !   begin body of subroutine GENPGMAT

    !.......  Initialize temporary storage table. NOTE - uses I/O API
    !           definitions of grid from include file
    IF( VFLAG ) THEN
        CALL GENPTVCEL( NPSRC, NGRID, XLOCA, YLOCA, NEXCLD, NX, INDX, GN, SN )
    ELSE
        CALL GENPTCEL(  NPSRC, NGRID, XLOCA, YLOCA, NEXCLD, NX, INDX, GN, SN )
    END IF

    !.......  Initialize statistics
    CMAX = NX( 1 )
    CMIN = CMAX

    !.......  Sort temporary gridding matrix
    CALL SORTI2( NPSRC, INDX, GN, SN )

    GFLAG = ( FNAME( 1:5 ) .NE. 'NONE ' )

    !.......   Compute gridding matrix and/or compute statistics

    K     = 0
    KSAV  = 0
    K2SAV = 0
    IF ( GFLAG ) THEN              ! write gridding matrix to file

        CALL M3MSG2( 'Computing gridding matrix and statistics...' )

        !.......   Compress matrix into I/O representation from scratch representation.
        !.......   Compute statistics

        DO R = 1, NGRID

            J = NX( R )

            IF( J .GT. CMAX ) THEN
                CMAX = J
            ELSEIF( J .LT. CMIN ) THEN
                CMIN = J
            ENDIF

            IF( J .NE. 0 ) THEN

                DO C = 1, J
                    K  = K + 1
                    K2 = K + NEXCLD

                    IF( K .LE. NPSRC .AND. K2 .LE. NPSRC ) IX( K ) = SN( INDX( K2 ) )

                    IF( K  .GT. KSAV )  KSAV  = K
                    IF( K2 .GT. K2SAV ) K2SAV = K2

                ENDDO

            ENDIF

        ENDDO        !  end of loop on cells K for this FIP

        !.......  Give error(s) if memory allocation exceeded
        IF( K .GT. NPSRC ) THEN
            WRITE( MESG,94010 )         &
                    'INTERNAL ERROR: Number of gridding coefficients K=', K,    &
                    'exceeds NPSRC=', NPSRC
            CALL M3MSG2( MESG )
        ENDIF

        IF( K2 .GT. NPSRC ) THEN
            WRITE( MESG,94010 )         &
                    'INTERNAL ERROR: Number of gridding coefficients K2=', K2,  &
                    'exceeds NPSRC=', NPSRC
            CALL M3MSG2( MESG )
        ENDIF

        !.......  Stop program if memory exceedance errors were written
        IF( K .GT. NPSRC .OR. K2 .GT. NPSRC )  CALL M3EXIT( PROGNAME, 0, 0, 'Overflow', 2 )

        CALL M3MSG2( 'Writing out GRIDDING MATRIX file...' )

        IF( .NOT. WRITE3( FNAME, 'ALL', 0, 0, NX ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, 'Error writing GRIDDING MATRIX file.', 2 )
        ENDIF

    ELSE                !....  not GFLAG:  report matrix stats only

        CALL M3MSG2( 'Computing gridding statistics...' )

        DO R = 1, NGRID      ! Must start loop at 1 to get correct K value

            J = NX( R )
            K = K + J

            IF( J .GT. CMAX ) THEN
                CMAX = J
            ELSEIF ( J .LT. CMIN ) THEN
                CMIN = J
            ENDIF

        ENDDO        !  end of loop on cells K for this FIP

    ENDIF      !  if gflag, or not

    NCOEF = K

    WRITE( MESG,94010 ) 'NOTE: Number of sources excluded from grid was', NEXCLD

    CALL M3MSG2( MESG )

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I9, :, 1X ) )

END SUBROUTINE GENPGMAT


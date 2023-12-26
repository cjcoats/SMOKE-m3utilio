
SUBROUTINE RDRMAT( FNAME, NSREAC, SPECNUM, IDX, REPEM,&
&                   PRJFC, MKTPN, RFAC )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine reads a reactivity matrix for any source category.
!      Certain variable names are expected to be in the files, so if the
!      reactivity matrix writer is changed, then this routine needs to be
!      changed as well.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Version ??/???? by ???
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
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

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!.........  SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: FNAME     ! speciation matrix file name
    INTEGER     , INTENT (IN) :: NSREAC    ! no. reactivity sources
    INTEGER     , INTENT (IN) :: SPECNUM   ! no. source-spec species (all)
    INTEGER     , INTENT(OUT) :: IDX  ( NSREAC ) ! source index
    REAL        , INTENT(OUT) :: REPEM( NSREAC ) ! replacement emissions
    REAL        , INTENT(OUT) :: PRJFC( NSREAC ) ! projection factors
    REAL        , INTENT(OUT) :: MKTPN( NSREAC ) ! market penetration
    REAL        , INTENT(OUT) :: RFAC ( NSREAC, SPECNUM ) ! spc coeffs

!...........   Error message strings
    CHARACTER(23), PARAMETER :: PART1 = 'Error reading variable '
    CHARACTER(23), PARAMETER :: PART3 = ' from REACTIVITY MATRIX'

!.........  Other local variables
    INTEGER            N, V       !  counters and indices

    CHARACTER(300)     MESG    ! message buffer
    CHARACTER(NAMLEN3) INVAR    ! tmp inventory pollutant name

    CHARACTER(16) :: PROGNAME = 'RDRMAT' ! program name

!***********************************************************************
!   begin body of subroutine RDRMAT

!.........  Retrieve file header
    IF ( .NOT. DESC3( FNAME ) ) THEN
        MESG = 'Could not get description of file ' // FNAME
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Loop through variables and read them
    V = 0
    DO N = 1, NVARS3D

        INVAR = VNAME3D( N )

        MESG = PART1 // INVAR( 1:LEN_TRIM( INVAR ) ) // PART3

!.............  Section for reading non-species variables
        SELECT CASE( INVAR )

          CASE( 'SRCID' )

            IF( .NOT. READ3( FNAME,INVAR,ALLAYS3,0,0,IDX ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ENDIF
            CYCLE

          CASE( 'REPEMIS' )

            IF( .NOT. READ3( FNAME,INVAR,ALLAYS3,0,0,REPEM ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ENDIF
            CYCLE

          CASE( 'PRJFAC' )

            IF( .NOT. READ3( FNAME,INVAR,ALLAYS3,0,0,PRJFC ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ENDIF
            CYCLE

          CASE( 'MKTPEN' )

            IF( .NOT. READ3( FNAME,INVAR,ALLAYS3,0,0,MKTPN ) ) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ENDIF
            CYCLE

        END SELECT

!.............  Section for reading species variables

        V = V + 1   ! Increment count for species-specific variables

        IF( V .LE. SPECNUM ) THEN

            IF( .NOT. READ3(FNAME,INVAR,ALLAYS3,0,0,RFAC(1,V))) THEN
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ENDIF

        END IF

    END DO  ! End loop on variables indicated for reading

    IF( V .GT. SPECNUM ) THEN

        WRITE( MESG,94010 )&
        &       'INTERNAL ERROR: Dimension mismatch. Memory ' //&
        &       'needed for reactivity factors was', V, CRLF() //&
        &       BLANK10 // 'but', SPECNUM, 'was allocated.'
        CALL M3MSG2( MESG )
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

    END IF

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I10, :, 1X ) )

END SUBROUTINE RDRMAT

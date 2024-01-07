
SUBROUTINE OPENCTMP( PKTTYP, IDEV )

    !***********************************************************************
    !  subroutine body starts at line 75
    !
    !  DESCRIPTION:
    !      This subroutine opens the temporary files that will contain the
    !      indices to the control packet data tables
    !
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version ??/???? by ???
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***************************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !       System
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

    !.....   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.....   SUBROUTINE ARGUMENTS:

    CHARACTER(*), INTENT (IN) :: PKTTYP       ! packet type
    INTEGER     , INTENT(OUT) :: IDEV         ! logical file number


    !.....   Other local variables

    INTEGER          IOS                       ! i/o status
    INTEGER, SAVE :: LP                        ! length of tmp file path

    LOGICAL, SAVE :: FIRSTIME = .TRUE.      ! true: first time subroutine called

    CHARACTER(200), SAVE :: PATHNM             ! path name for tmp file
    CHARACTER(220)          FILENM                    ! file name
    CHARACTER(256)          MESG                      ! message buffer
    CHARACTER(16), PARAMETER :: PROGNAME = 'OPENCTMP'     ! program name

    !***********************************************************************
    !   Begin body of subroutine OPENCTMP

    IF( FIRSTIME ) THEN

        FIRSTIME = .FALSE.

        MESG = 'Path where temporary control files will be written'
        CALL ENVSTR( 'SMK_TMPDIR', MESG, '.', PATHNM, IOS )
        LP = LEN_TRIM( PATHNM )

        IF( IOS .NE. 0 ) THEN
            MESG = 'WARNING: Large temporary files being placed '// &
                   'executable directory because ' // CRLF() //     &
                   BLANK10 // 'environment variable SMK_TMPPATH '// &
                   'is not set properly'
            CALL M3MSG2( MESG )
        END IF

    END IF

    SELECT CASE( PKTTYP )

      CASE ( 'CTG' )

        FILENM = PATHNM( 1:LP ) // '/cntlmat_tmp_ctg'
        IDEV = JUNIT()
        OPEN( IDEV, ERR=1006, FILE=FILENM )

      CASE ( 'CONTROL', 'EMS_CONTROL' )

        FILENM = PATHNM( 1:LP ) // '/cntlmat_tmp_ctl'
        IDEV = GETEFILE( FILENM, .FALSE., .TRUE., PROGNAME )

      CASE ( 'ALLOWABLE' )

        FILENM = PATHNM( 1:LP ) // '/cntlmat_tmp_alw'
        IDEV = GETEFILE( FILENM, .FALSE., .TRUE., PROGNAME )

      CASE ( 'PROJECTION' )

        FILENM = PATHNM( 1:LP ) // '/cntlmat_tmp_proj'
        IDEV = GETEFILE( FILENM, .FALSE., .TRUE., PROGNAME )

      CASE ( 'MACT' )

        FILENM = PATHNM( 1:LP ) // '/cntlmat_tmp_mact'
        IDEV = GETEFILE( FILENM, .FALSE., .TRUE., PROGNAME )

    END SELECT

    IF ( IDEV .LT. 0 ) THEN
        MESG = 'Could not open '//FILENM
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    RETURN

    !.....  Open file error
1006 MESG = 'Could not open temporary ASCII file for control info.'
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

END SUBROUTINE OPENCTMP

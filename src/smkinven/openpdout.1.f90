
SUBROUTINE OPENPDOUT( NPDSRC, NVAR, TZONE, SDATE, STIME, TSTEP,&
&                      FILFMT, TYPNAM, PFLAG, EAIDX,  SPSTAT,&
&                      FNAME, RDEV )

!***********************************************************************
!  subroutine body starts at line 96
!
!  DESCRIPTION:
!      This subroutine opens the day- or hour-specific output files
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!      Subroutines: I/O API subroutines
!
!  REVISION  HISTORY:
!      Created 12/99 by M. Houyoux
!
!       Version 10/2016 by C. Coats:  USE M3UTILIO
!
!****************************************************************************
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
!.........  This module contains the information about the source category
    USE MODINFO, ONLY: CATDESC, CRL, BYEAR, NIPOL, NIACT,&
    &                   EANAM

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters

!...........   EXTERNAL FUNCTIONS and their descriptions
    CHARACTER(16), EXTERNAL :: VERCHAR

!...........   SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN) :: NPDSRC    ! no. period-specific sources
    INTEGER     , INTENT (IN) :: NVAR      ! no. output variables
    INTEGER     , INTENT (IN) :: TZONE     ! time zone of date/time
    INTEGER     , INTENT (IN) :: SDATE     ! Julian start date
    INTEGER     , INTENT (IN) :: STIME     ! start time HHMMSS
    INTEGER     , INTENT (IN) :: TSTEP     ! time step HHMMSS
    INTEGER     , INTENT (IN) :: FILFMT    ! format of period-specific data
    CHARACTER(*), INTENT (IN) :: TYPNAM    ! 'day' or 'hour'
    LOGICAL     , INTENT (IN) :: PFLAG     ! true: creating profiles
    INTEGER     , INTENT (IN) :: EAIDX( NVAR ) ! pol/act index
    INTEGER     , INTENT (IN) :: SPSTAT( MXSPDAT ) ! true: special val exists
    CHARACTER(*), INTENT(OUT) :: FNAME     ! logical file name
    INTEGER     , INTENT(IN OUT) :: RDEV   ! report unit number

!...........   LOCAL PARAMETERS
    CHARACTER(50), PARAMETER ::&
    &CVSW = '$Name SMOKEv5.0_Jun2023$' ! CVS release tag

!...........   Other local variables

    INTEGER       J, K, V   ! counter and indices

    CHARACTER(5)       CTZONE      ! string of time zone
    CHARACTER(5)       SCRNAM      ! upcase(typnam)
    CHARACTER(NAMLEN3) VARNAM      ! name for integer index
    CHARACTER(NAMLEN3) VBUF        ! tmp buffer for variable names
    CHARACTER(300)     MESG        ! message buffer

    CHARACTER(16) :: PROGNAME = 'OPENPDOUT' ! program name

!***********************************************************************
!   begin body of subroutine OPENPDOUT

!.........  Write time zone to character string
    WRITE( CTZONE,94000 ) TZONE

!.........  Determine integer index variable name
    IF( TYPNAM .EQ. 'day' ) THEN
        VARNAM = 'INDXD'
    ELSE IF( TYPNAM .EQ. 'hour' ) THEN
        VARNAM = 'INDXH'
    ELSE
        MESG = 'INTERNAL ERROR: Type name "'// TYPNAM //&
        &       '" not recognized'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Set header up with missing
    CALL HDRMISS3

!.........  Set of header with actual settings
    SDATE3D = SDATE
    STIME3D = STIME
    TSTEP3D = TSTEP
    NROWS3D = NPDSRC     !  number of rows = # of period sources

!.........  Define variable for source index
    VNAME3D( 1 ) = VARNAM
    VTYPE3D( 1 ) = M3INT
    UNITS3D( 1 ) = 'n/a'
    VDESC3D( 1 ) = 'source IDs'

!.........  Define hour-specific emissions, if any

    J = 1
    DO V = 1, NVAR
        VBUF = EANAM( EAIDX( V ) )
        J = J + 1
        VNAME3D( J ) = VBUF
        VTYPE3D( J ) = M3REAL

!............. If outputs are profiles instead of data values
        IF( PFLAG ) THEN
            UNITS3D( J ) = 'n/a'
            VDESC3D( J ) = TYPNAM // '-specific ' //&
            &                   TRIM( VBUF ) // ' diurnal profile'

!............. If outputs are data values...
        ELSE
            UNITS3D( J ) = 'ton/' // TYPNAM
            VDESC3D( J ) = TYPNAM // '-specific ' //&
            &                   TRIM( VBUF ) // ' data'
        END IF

    END DO

!.........  Define hour-specific special data values, if any
    K = 0
    DO V = 1, MXSPDAT
        IF( SPSTAT( V ) > 0 ) THEN
            K = K + 1
            J = J + 1
            VNAME3D( J ) = SPDATNAM( V )
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = SPDATUNT( V )

            VDESC3D( J ) = TYPNAM // '-specific ' //&
            &               TRIM( SPDATDSC( V ) ) // ' data'

        END IF
    END DO
    NVARS3D = J

!.........  Define general description info
    FDESC3D( 1 ) = CATDESC // TRIM( TYPNAM ) //&
    &              '-specific source inventory'
    FDESC3D( 2 ) = '/FROM/ ' // PROGNAME
    FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )

    IF( NIPOL .GT. 0 ) THEN
        WRITE( FDESC3D( 4 ),94010 ) '/POLLUTANTS/', NIPOL
    END IF

    IF( NIACT .GT. 0 ) THEN
        WRITE( FDESC3D( 5 ),94010 ) '/ACTIVITIES/', NIACT
    END IF

    IF( K .GT. 0 ) THEN
        WRITE( FDESC3D( 6 ),94010 ) '/SPECIAL DATA/', K
    END IF

    WRITE( FDESC3D( 7 ),94010 ) '/BASE YEAR/ '    , BYEAR
    FDESC3D( 8 ) = '/TZONE/ ' // CTZONE

!.........  Set up default file name and prompting message
    SCRNAM = TYPNAM
    CALL UPCASE( SCRNAM )
    IF( PFLAG ) THEN
        MESG = 'Enter logical name for ' // TRIM( SCRNAM ) //&
        &        '-SPECIFIC PROFILES output file'
        FNAME = CRL // TYPNAM // 'PRO'

    ELSE
        MESG = 'Enter logical name for ' // TRIM( SCRNAM ) //&
        &        '-SPECIFIC output file'
        FNAME = CRL // TYPNAM
    END IF

!.........  Prompt for output file
    CALL UPCASE( FNAME )
    FNAME = PROMPTMFILE( MESG, FSUNKN3, FNAME, PROGNAME )

!.........  If format is CEM format, prompt for report output file name
    IF ( FILFMT .EQ. CEMFMT .AND. RDEV .LE. 0 ) THEN
        MESG = 'Enter logical name for the CEM MATCHING REPORT'
        RDEV = PROMPTFFILE( MESG, .FALSE., .TRUE.,&
        &                    'REPINVEN', PROGNAME )

    END IF

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94000 FORMAT( I2.2 )

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE OPENPDOUT


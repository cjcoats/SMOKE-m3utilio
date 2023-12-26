
SUBROUTINE OPENGRWOUT( ENAME, FYEAR, NAME1, SFLAG,&
&                       OFLAG, ODEV, DDEV, VDEV, RDEV,&
&                       ONAME, VARPATH )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine opens files needed by the Grwinven program.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 2/2000 by M. Houyoux
!       Updated 6/8/2005 by M. Houyoux for ORL format
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
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

!.........  MODULES for public variables
!.........  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, NIPOL, NIACT

    IMPLICIT NONE

!...........   INCLUDES

    INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
    INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

!...........   SUBROUTINE ARGUMENTS
    CHARACTER(16), INTENT (IN) :: ENAME ! emis input inven logical name
    INTEGER      , INTENT (IN) :: FYEAR ! future year or 0 for no projection
    CHARACTER(80), INTENT (IN) :: NAME1 ! physical name part for i/o api
    LOGICAL      , INTENT (IN) :: SFLAG ! true: output SMOKE file
    LOGICAL      , INTENT (IN) :: OFLAG ! true: output ORL file
    INTEGER      , INTENT(OUT) :: ODEV  ! output map inventory file
    INTEGER      , INTENT(OUT) :: DDEV  ! IDA output emissions file number
    INTEGER      , INTENT(OUT) :: VDEV  ! IDA output activity file number
    INTEGER      , INTENT(OUT) :: RDEV  ! ORL output emissions file number
    CHARACTER(16), INTENT(OUT) :: ONAME ! output logical main i/o api
    CHARACTER(PHYLEN3), INTENT( OUT ) :: VARPATH ! path for pol/act output files

!...........   EXTERNAL FUNCTIONS and their descriptions

    CHARACTER(MXDLEN3), EXTERNAL :: GETCFDSC
    INTEGER           , EXTERNAL :: GETIFDSC
    CHARACTER(16)     , EXTERNAL :: VERCHAR

!...........   LOCAL PARAMETERS
    CHARACTER(16), PARAMETER :: PROGNAME = 'OPENGRWOUT' ! program name
    CHARACTER(50), PARAMETER ::&
    &CVSW = '$Name SMOKEv5.0_Jun2023$' ! CVS release tag

!...........   Other local variables

    INTEGER         BYEAR     ! tmp base year
    INTEGER         CATIDX    ! index for source category
    INTEGER         IOS       ! i/o status
    INTEGER         L, N         ! counters and indices

    CHARACTER(16)   ANAME     ! tmp dummy buffer
    CHARACTER(16)   INAME     ! tmp output IDA file name
    CHARACTER(16)   NNAME     ! tmp output map logical file name
    CHARACTER(30)   EVNAME    ! tmp environment variable name
    CHARACTER(256)  MESG      ! message buffer
    CHARACTER(PHYLEN3) APHYS   ! ASCII physical file name
    CHARACTER(PHYLEN3) EPHYS   ! I/O API physical file name
    CHARACTER(PHYLEN3) PATH    ! path name
    CHARACTER(MXDLEN3) IFDESC2, IFDESC3 ! fields 2 & 3 from inven FDESC

!***********************************************************************
!   begin body of subroutine OPENGRWOUT

!.........  Initialize outputs
    ONAME = 'NONE'
    DDEV  = 0
    VDEV  = 0
    RDEV  = 0

!.........  Get inventory file names given source category
    CALL GETINAME( CATEGORY, NNAME, ANAME )

!.........  Open map-formatted output inventory file without prompting
    ODEV = 0
    IF( SFLAG ) THEN
        NNAME = TRIM( NNAME ) // '_O'
        ODEV = GETEFILE( NNAME, .FALSE., .TRUE., PROGNAME )
        IF ( ODEV .LT. 0 ) THEN     !  failure to open

            MESG = 'Could not open INVENTORY MAP file:' // CRLF() //&
            &        BLANK10 // TRIM( NNAME ) // '.'
            CALL M3MSG2( MESG )

            MESG = 'Ending program.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF      !  if getefile() failed

!.............  Evaluate physical file name of inventory map
        MESG = 'Inventory map file name'
        CALL ENVSTR( NNAME, MESG, BLANK16, APHYS, IOS )

        IF( IOS .NE. 0 ) THEN
            MESG = 'Unable to evaluate environment variable "' //&
            &       TRIM( NNAME ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!.............  Determine path of inventory map
        L = LEN_TRIM( APHYS )
        DO N = L, 1, -1

            IF( APHYS( N:N ) .EQ. '/' .OR.&
            &    APHYS( N:N ) .EQ. '\'      ) THEN
                PATH = APHYS( 1:N )
                EXIT
            END IF

        END DO

!.............  Re-read header for input inventory file to initialize
!               header for output inventory file
        IF( .NOT. DESCSET( ENAME,ALLFILES ) ) THEN
            MESG = 'Could not read description for "' //&
            &       TRIM( ENAME ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!.............  Store non-category-specific header information
        IFDESC2 = GETCFDSC( FDESC3D, '/FROM/', .TRUE. )
        IFDESC3 = GETCFDSC( FDESC3D, '/VERSION/', .TRUE. )
        BYEAR   = GETIFDSC( FDESC3D, '/BASE YEAR/', .TRUE. )

        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )

!.............  If projection year is non-zero and we haven't simply
!               applied a factor to the inventory for the same year,
!               store future year in FDESC3D
        IF( FYEAR .NE. 0 .AND. FYEAR .NE. BYEAR ) THEN
            WRITE( FDESC3D( 14 ),94010 ) '/PROJECTED YEAR/ ', FYEAR
        END IF

        FDESC3D( 15 ) = '/INVEN FROM/ ' // IFDESC2
        FDESC3D( 16 ) = '/INVEN VERSION/ ' // IFDESC3

!.............  Build output file physical name
        ONAME = 'IOAPI_INV'
        EPHYS = TRIM( PATH ) // TRIM( NAME1 ) // '.ncf'

!.............  Set output logical file name
        IF( .NOT. SETENVVAR( ONAME, EPHYS ) ) THEN
            MESG = 'Could not set logical file name for ' //&
            &       'file:'// CRLF() // BLANK10 // TRIM( EPHYS )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

!............  Open I/O API file
        ELSE IF( .NOT. OPENSET( ONAME, FSNEW3, PROGNAME ) ) THEN
            MESG = 'Could not open I/O API inventory file '//&
            &       'for file name:' // CRLF() // BLANK10 //&
            &       TRIM( EPHYS )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

    END IF

!.........  Provide variable path
    VARPATH = TRIM( PATH ) // TRIM( NAME1 ) // '_dat'

!.........  Get index for source category to use for output file names
    CATIDX = INDEX1( CATEGORY, NCAT, CATLIST )

!.........  Prompt for emissions ORL file
    IF( OFLAG ) THEN
        MESG  = 'Enter logical name for the ORL EMISSIONS ' //&
        &        'output file'
        INAME = ANAMLIST( CATIDX ) // '_O'

        RDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., INAME, PROGNAME )
    END IF

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE OPENGRWOUT


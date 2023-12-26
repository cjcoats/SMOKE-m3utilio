
SUBROUTINE OPENRMAT( ENAME, RPOL, SFLAG, LFLAG,&
&                     BYEARIN, PYEAR, NSREAC, NMSPC, SPECIES,&
&                     SDEV, SNAME, LNAME, SVNAMES, LVNAMES )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      Open the mass-based and mole-based reactivity matrices
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 3/99 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
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

!.........  MODULES for public variables
!.........  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, CATLEN, CRL, NSRC

!.........This module is required by the FileSetAPI
    USE MODFILESET

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
    INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

!...........   EXTERNAL FUNCTIONS and their descriptions:
    CHARACTER(IODLEN3), EXTERNAL :: GETCFDSC
    INTEGER           , EXTERNAL :: GETIFDSC
    CHARACTER(16)     , EXTERNAL :: VERCHAR

!...........   LOCAL PARAMETERS
    INTEGER      , PARAMETER :: NBASVAR = 4
    CHARACTER(16), PARAMETER:: PROGNAME = 'OPENRMAT' ! program name
    CHARACTER(50), PARAMETER :: CVSW = '$Name SMOKEv5.0_Jun2023$' ! CVS release tag

!.........  SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: ENAME      ! emissions inven logical name
    CHARACTER(*), INTENT (IN) :: RPOL       ! pollutant for matrices
    LOGICAL     , INTENT (IN) :: SFLAG      ! true: open mass-based file
    LOGICAL     , INTENT (IN) :: LFLAG      ! true: open mole-based file
    INTEGER     , INTENT (IN) :: BYEARIN    ! base year of proj factors
    INTEGER     , INTENT (IN) :: PYEAR      ! projected year of proj factors
    INTEGER     , INTENT (IN) :: NSREAC     ! number of reactvty sources
    INTEGER     , INTENT (IN) :: NMSPC      ! number of reactivity species
    CHARACTER(*), INTENT (IN) :: SPECIES( NMSPC ) ! model spec nams
    INTEGER     , INTENT(OUT) :: SDEV       ! unit number of supplement file
    CHARACTER(*), INTENT(OUT) :: SNAME      ! mass-based spec file name
    CHARACTER(*), INTENT(OUT) :: LNAME      ! mole-based spec file name
    CHARACTER(*), INTENT(OUT) :: SVNAMES( NMSPC )   ! species mass out vars
    CHARACTER(*), INTENT(OUT) :: LVNAMES( NMSPC )   ! species mole out vars

!.........  Other local variables
    INTEGER          I, J           !  counters and indices
    INTEGER          IOS            !  i/o status

    CHARACTER(NAMLEN3) NAMBUF   ! file name buffer
    CHARACTER(300)     MESG     ! message buffer

    CHARACTER(IOVLEN3) CPOL     ! pollutant name buffer
    CHARACTER(IODLEN3) IFDESC2, IFDESC3 ! fields 2 & 3 from inven FDESC

!***********************************************************************
!   begin body of subroutine OPENRMAT

!.........  Get header information from inventory file

    IF ( .NOT. DESCSET( ENAME,-1 ) ) THEN
        MESG = 'Could not get description of file "'&
        &       // TRIM( ENAME ) // '".'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IFDESC2 = GETCFDSC( FDESC3D, '/FROM/', .TRUE. )
    IFDESC3 = GETCFDSC( FDESC3D, '/VERSION/', .TRUE. )

!.........  Initialize variable names for I/O API files
    SVNAMES = ' '  ! Array
    LVNAMES = ' '  ! Array

!.........  Initialize I/O API output file headers
    CALL HDRMISS3

!.........  Set I/O API header parms that need values
    NVARSET = MIN( NBASVAR + NMSPC, MXVARS3 )
    NROWS3D = NSREAC
    NTHIK3D = NSRC

    FDESC3D( 1 ) = CATEGORY( 1:CATLEN ) // ' reactivity matrix'
    FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
    FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )

    WRITE( FDESC3D( 5 ), '(A,I4)' ) '/BASE YEAR/ ', BYEARIN
    WRITE( FDESC3D( 6 ), '(A,I4)' ) '/PROJECTED YEAR/ ', PYEAR
    WRITE( FDESC3D( 7 ), '(A,I4)' ) '/SPECIES VARS/ ', NMSPC
    WRITE( FDESC3D( 8 ), '(A,I4)' ) '/CTYPE/ ', CTYPREAC

    FDESC3D( 11 ) = '/INVEN FROM/ ' // IFDESC2
    FDESC3D( 12 ) = '/INVEN VERSION/ ' // IFDESC3

    IF( ALLOCATED( VTYPESET ) )&
    &    DEALLOCATE( VTYPESET, VNAMESET, VUNITSET, VDESCSET )
    ALLOCATE( VNAMESET( NVARSET ),&
    &          VTYPESET( NVARSET ),&
    &          VUNITSET( NVARSET ),&
    &          VDESCSET( NVARSET ), STAT=IOS )
    CALL CHECKMEM( IOS, 'VNAMESET...VDESCSET', PROGNAME )

!.........  Also deallocate the number of variables per file so
!           that this will be set automatically by openset
    DEALLOCATE( VARS_PER_FILE )

!.........  Set up non-speciation variables
    J = 1
    VNAMESET( J ) = 'SRCID'
    VTYPESET( J ) = M3INT
    VUNITSET( J ) = 'n/a'
    VDESCSET( J ) = 'Inventory source ID (position)'
    J = J + 1

    VNAMESET( J )= 'REPEMIS'
    VTYPESET( J )= M3REAL
    VUNITSET( J )= 'tons/day'
    VDESCSET( J )= 'Reactivity base-year emissions'
    J = J + 1

    VNAMESET( J )= 'PRJFAC'
    VTYPESET( J )= M3REAL
    VUNITSET( J )= 'n/a'
    VDESCSET( J )= 'Reactivity projection factor'

    J = J + 1
    VNAMESET( J ) = 'MKTPEN'
    VTYPESET( J ) = M3REAL
    VUNITSET( J ) = 'fraction'
    VDESCSET( J ) = 'Reactivity control market penetration'

!.........  Make sure program has not been modified improperly to cause a
!           hard-to-detect error

    IF( J .NE. NBASVAR ) THEN
        MESG = 'INTERNAL ERROR: Number of variables NBASVAR is ' //&
        &       'inconsistent with other subroutine modifications'
        CALL M3MSG2( MESG )
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
    END IF

!.........  Set speciation variables for I/O API header
    CPOL = ADJUSTL( RPOL )
    DO I = 1, NMSPC

        J = J + 1

!.............  Set SVNAMES and LVNAMES here instead of VNAME3D, so that it
!               will be easier to change these to have different name roots
!               if needed in the future
!.............  Check total number of output variables with I/O API max.
        IF( J .LE. MXVARS3 ) THEN
            VTYPESET( J ) = M3REAL
            VDESCSET( J ) = CPOL // SPJOIN // SPECIES( I )
        ENDIF

        WRITE( SVNAMES( I ), '(A4,I3.3)' ) 'SVAR', I
        WRITE( LVNAMES( I ), '(A4,I3.3)' ) 'SVAR', I

    END DO

!.........  Error if number of variables is passed maximum because we can't
!           store the names of the variables.
!.........  DO NOT end program here because it will be ended when the write
!           attempt is made for these extra variables.
    IF( J .GT. MXVARS3 ) THEN

        WRITE( MESG, 94010 )&
        &       'Maximum I/O API variables exceeded:' //&
        &       CRLF() // BLANK10 // 'Max: ', MXVARS3, 'Actual:', J
        CALL M3MSG2( MESG )

    ENDIF

    MESG = 'Enter logical name(s) for ...'
    CALL M3MSG2( MESG )

!.........  Set up variables specifically for mass-based file, and open it
    IF( SFLAG ) THEN

        FDESC3D( 4 ) = '/SMATTYPE/ ' // ' Mass'

        I = 0
        DO J = NBASVAR+1, NVARSET

            I = I + 1
            VNAMESET( J ) =  SVNAMES( I )
            VUNITSET( J ) = 'gm/ton'

        END DO

        MESG = 'I/O API MASS-BASED REACTIVITY MATRIX for ' // RPOL

        NAMBUF = PROMPTSET( MESG, FSUNKN3, CRL // 'RMAT_S',&
        &                    PROGNAME )
        SNAME = NAMBUF

    ENDIF

!.........  Set up variables specifically for mole-based file, and open it
    IF( LFLAG ) THEN

        FDESC3D( 4 ) = '/SMATTYPE/ ' // ' Mole'

        I = 0
        DO J = NBASVAR+1, NVARSET

            I = I + 1
            VNAMESET( J ) =  LVNAMES( I )
            VUNITSET( J ) = 'mole/ton'

        END DO

        MESG = 'I/O API MOLE-BASED REACTIVITY MATRIX for ' // RPOL

        NAMBUF = PROMPTSET( MESG, FSUNKN3, CRL // 'RMAT_L',&
        &                    PROGNAME )
        LNAME = NAMBUF

    ENDIF

!.........  Open the supplementary file (for SCCs and SPROFs)

    MESG = 'ASCII REACTIVITY MATRIX SUPPLEMENT file for ' // RPOL

    SDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., CRL // 'RSUP',&
    &                    PROGNAME )

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE OPENRMAT


SUBROUTINE OPENPMAT( ENAME, BYEARIN, PYEAR, PNAME )

    !***********************************************************************
    !  subroutine body starts at line 87
    !
    !  DESCRIPTION:
    !      Open the projection matrix.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version ??/???? by ???
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***********************************************************************
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

    !.......  MODULES for public variables
    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, CATLEN, CRL, NSRC

    !.......  This module contains the control packet data and control matrices
    USE MODCNTRL, ONLY: POLSFLAG, NVPROJ, PNAMPROJ

    !.......This module is required by the FileSetAPI
    USE MODFILESET

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters
    INCLUDE 'SETDECL.h90'       !  FileSetAPI variables and functions

    !.......   EXTERNAL FUNCTIONS and their descriptions:
    CHARACTER(IODLEN3), EXTERNAL :: GETCFDSC
    INTEGER           , EXTERNAL :: GETIFDSC
    CHARACTER(16)     , EXTERNAL :: VERCHAR

    !.......   LOCAL PARAMETERS

    CHARACTER(16), PARAMETER :: PROGNAME = 'OPENPMAT'     ! program name
    CHARACTER(50), PARAMETER ::     CVSW = '$Name SMOKEv5.0_Jun2023$'     ! CVS release tag

    !.......  SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: ENAME          ! emissions inven logical name
    INTEGER     , INTENT (IN) :: BYEARIN        ! base year of proj factors
    INTEGER     , INTENT (IN) :: PYEAR          ! projected year of proj factors
    CHARACTER(*), INTENT(OUT) :: PNAME          ! projection file name

    !.......  Other local variables
    INTEGER          I, J               !  counters and indices
    INTEGER          IOS                !  i/o status

    CHARACTER(NAMLEN3) NAMBUF       ! file name buffer
    CHARACTER(300)     MESG         ! message buffer

    CHARACTER(IODLEN3) IFDESC2, IFDESC3     ! fields 2  3 from inven FDESC

    !***********************************************************************
    !   begin body of subroutine OPENPMAT

    !.......  Get header information from inventory file

    IF ( .NOT. DESCSET( ENAME,-1 ) ) THEN
        MESG = 'Could not get description of file "' // TRIM( ENAME ) // '".'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IFDESC2 = GETCFDSC( FDESC3D, '/FROM/', .TRUE. )
    IFDESC3 = GETCFDSC( FDESC3D, '/VERSION/', .TRUE. )

    !.......  Initialize I/O API output file headers
    CALL HDRMISS3

    !.......  Set I/O API header parms that need values
    NROWS3D = NSRC

    FDESC3D( 1 ) = CATEGORY( 1:CATLEN ) // ' projection matrix'
    FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
    FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )

    WRITE( FDESC3D( 4 ), '(A,I4)' ) '/CTYPE/ ', CTYPPROJ
    WRITE( FDESC3D( 5 ), '(A,I4)' ) '/BASE YEAR/ ', BYEARIN
    WRITE( FDESC3D( 6 ), '(A,I4)' ) '/PROJECTED YEAR/ ', PYEAR

    FDESC3D( 11 ) = '/INVEN FROM/ ' // IFDESC2
    FDESC3D( 12 ) = '/INVEN VERSION/ ' // IFDESC3

    IF( ALLOCATED( VTYPESET ) ) DEALLOCATE( VTYPESET, VNAMESET, VUNITSET, VDESCSET )

    !.......  Set number of variables in output file based on whether pollutant-specific
    !           assignments are being used
    IF( POLSFLAG ) THEN
        NVARSET = NVPROJ
    ELSE
        NVARSET = 1
    END IF

    ALLOCATE( VNAMESET( NVARSET ),      &
              VTYPESET( NVARSET ),      &
              VUNITSET( NVARSET ),      &
              VDESCSET( NVARSET ), STAT=IOS )
    CALL CHECKMEM( IOS, 'VNAMESET...VDESCSET', PROGNAME )

    !.......  Also deallocate the number of variables per file so
    !           that this will be set automatically by openset
    DEALLOCATE( VARS_PER_FILE )

    !.......  If pollutant-specific assignments, then set up projection matrix
    !           with one variable for each pollutant being projected
    IF( POLSFLAG ) THEN
        DO J = 1, NVPROJ
            VNAMESET( J )= PNAMPROJ( J )      ! Lowercase used to permit inv data named "PFAC"
            VTYPESET( J )= M3REAL
            VUNITSET( J )= 'n/a'
            VDESCSET( J )= 'Projection factor for ' // PNAMPROJ( J )
        END DO

    !.......  If no pollutant-specific assignments, then set up projection
    !           matrix to reflect that all pollutants affected the
    ELSE
        J = 1
        VNAMESET( J )= 'pfac'      ! Lowercase used to permit inv data named "PFAC"
        VTYPESET( J )= M3REAL
        VUNITSET( J )= 'n/a'
        VDESCSET( J )= 'Projection factor'
    END IF

    MESG = 'Enter logical name for projection matrix...'
    CALL M3MSG2( MESG )

    !.......  Open projection matrix.
    !.......  Using NAMBUF is needed for HP to ensure string length consistencies

    MESG = 'I/O API PROJECTION MATRIX'

    NAMBUF = PROMPTSET( MESG, FSUNKN3, CRL // 'PMAT',&
                        PROGNAME )
    PNAME = NAMBUF

    RETURN

END SUBROUTINE OPENPMAT

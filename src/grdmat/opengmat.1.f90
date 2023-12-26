
SUBROUTINE OPENGMAT( NMATX, NMATXU, UFLAG, INVPROG, INVVERS,&
&                     VFLAG, NFLAG, GNAME, UNAME, FDEV )

!***********************************************************************
!  subroutine body starts at line 95
!
!  DESCRIPTION:
!      This subroutine sets up the header and variables for the gridding matrix
!      and ungridding matrix for mobile sources.  It prompts for the file
!      names and opens the files.
!
!  PRECONDITIONS REQUIRED:
!      Correct number of sources and matrix dimensions are set.  Grid parameters
!      are defined and passed the FDESC3D include file.
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 5/99 by M. Houyoux
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
    USE MODINFO, ONLY: CATEGORY, CRL, CATDESC, NSRC

    IMPLICIT NONE

!...........   INCLUDES

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN) :: NMATX   ! no. of source-cell intersections
    INTEGER     , INTENT (IN) :: NMATXU  ! no. of county-cell intersections (for all sources)
    LOGICAL     , INTENT (IN) :: UFLAG   ! true: open ungridding matrix
    CHARACTER(*), INTENT (IN) :: INVPROG ! inventory program
    CHARACTER(*), INTENT (IN) :: INVVERS ! inventory program version
    LOGICAL     , INTENT (IN) :: VFLAG   ! true: using variable grid
    LOGICAL     , INTENT (IN) :: NFLAG   ! true: using raw netcdf gridded inv file
    CHARACTER(*), INTENT(OUT) :: GNAME   ! gridding matrix logical name
    CHARACTER(*), INTENT(OUT) :: UNAME   ! ungridding matrix logical name
    INTEGER     , INTENT(OUT) :: FDEV    ! report file

!...........   EXTERNAL FUNCTIONS and their descriptions
    LOGICAL      , EXTERNAL :: DSCM3GRD
    CHARACTER(16), EXTERNAL :: VERCHAR

!...........   LOCAL PARAMETERS
    CHARACTER(16), PARAMETER :: PROGNAME = 'OPENGMAT' ! program name
    CHARACTER(50), PARAMETER :: CVSW =&
    &     '$Name SMOKEv5.0_Jun2023$' ! CVS release tag

!...........   Other local variables

    INTEGER       I, J, L1, L2, V     ! counter and indices
    INTEGER       NCOLS     ! number of columns
    INTEGER       NIOVARS   ! number of I/O API file non-emis variables
    INTEGER       NPOLMAX   ! max no of pollutants, based on I/O API
    INTEGER       NNPVAR    ! no. non-pollutant inventory variables
    INTEGER       NROWS     ! number of rows

    CHARACTER(16)  COORD3D   ! coordinate system name
    CHARACTER(16)  COORUN3D  ! coordinate system projection units
    CHARACTER(80)  GDESC     ! grid description
    CHARACTER(16)  GDNAME    ! grid name
    CHARACTER(256) MESG      ! message buffer

    CHARACTER(NAMLEN3) NAMBUF   ! file name buffer

!***********************************************************************
!   begin body of subroutine OPENGMAT

!.........  Initialize header by getting the Models-3 grid information file
    IF( .NOT. DSCM3GRD( GDNAME, GDESC, COORD3D, GDTYP3D, COORUN3D,&
    &                    P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,&
    &                    XORIG3D, YORIG3D, XCELL3D, YCELL3D,&
    &                    NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

        MESG = 'Could not get Models-3 grid description.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Set up I/O API output file header for gridding matrix
    NCOLS   = NCOLS3D
    NROWS   = NROWS3D
    FTYPE3D = SMATRX3
    NCOLS3D = NMATX
    NROWS3D = NCOLS * NROWS
    NTHIK3D = NSRC

    IF( CATEGORY .EQ. 'POINT' ) THEN
        NVARS3D = 0
    ELSE
        NVARS3D = 1
        VNAME3D( 1 ) = CRL // 'GRDMAT'
        UNITS3D( 1 ) = 'n/a'
        VDESC3D( 1 ) = CATDESC // ' source gridding coefficients'
        VTYPE3D( 1 ) = M3REAL
    END IF

    FDESC3D( 1  ) = CATDESC // 'source gridding matrix'
    FDESC3D( 2  ) = '/FROM/ ' // PROGNAME
    FDESC3D( 3  ) = '/VERSION/ ' // VERCHAR( CVSW )
    FDESC3D( 4  ) = '/GDESC/ ' // GDESC
    WRITE( FDESC3D(5), 94010 ) '/NCOLS3D/ ', NCOLS
    WRITE( FDESC3D(6), 94010 ) '/NROWS3D/ ', NROWS

    FDESC3D( 11 ) = '/INVEN FROM/ ' // INVPROG
    FDESC3D( 12 ) = '/INVEN VERSION/ ' // INVVERS

    IF( VFLAG ) THEN
        FDESC3D( 13 ) = '/VARIABLE GRID/ ' // GDNAME
    END IF

!.........  Get name of gridding matrix
    IF( .NOT. NFLAG ) THEN
        NAMBUF = PROMPTMFILE(&
        &     'Enter logical name for GRIDDING MATRIX output file',&
        &     FSUNKN3, CRL // 'GMAT', PROGNAME )
        GNAME = NAMBUF
    ELSE
        GNAME = CRL // 'GMAT'
    END IF

!.........  Set up I/O API output file header for ungridding matrix
!.........  Leave everything the same as for the gridding matrix, but a couple
!           of header items

    UNAME  = 'NONE'
    IF( UFLAG ) THEN

        NCOLS3D = NMATXU
        NROWS3D = NSRC
        NTHIK3D = NCOLS * NROWS
        FDESC3D( 1 ) = CATDESC // ' source ungridding matrix'

        VNAME3D( 1 ) = CRL // 'UGRDMAT'
        UNITS3D( 1 ) = 'n/a'
        VDESC3D( 1 ) = CATDESC // ' source ungridding coefficients'
        VTYPE3D( 1 ) = M3REAL

        NAMBUF = PROMPTMFILE(&
        &     'Enter logical name for UNGRIDDING MATRIX output file',&
        &     FSUNKN3, CRL // 'UMAT', PROGNAME )
        UNAME = NAMBUF

    END IF

!.........  Open report file
    IF( CATEGORY .EQ. 'AREA'   .OR.&
    &    CATEGORY .EQ. 'MOBILE'      ) THEN

        MESG = 'Enter logical name for the GRIDDING SUPPLEMENTAL '//&
        &       'file'
        FDEV = PROMPTFFILE( MESG, .FALSE., .TRUE.,&
        &                    CRL // 'GSUP', PROGNAME )

    END IF

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE OPENGMAT


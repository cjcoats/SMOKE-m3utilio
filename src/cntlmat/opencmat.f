
        SUBROUTINE OPENCMAT( ENAME, MATTYP, MNAME )

C***********************************************************************
C  subroutine body starts at line 88
C
C  DESCRIPTION:
C      This subroutine opens the file in which the control matrix
C      will be written.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Version ??/???? by ???
C       Version 11/2023 by CJC:  USE M3UTILIO and related changes
C
C***************************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
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

C.........  MODULES for public variables
C.........  This module contains the control packet data and control matrices
        USE MODCNTRL, ONLY: NVCMULT, PNAMMULT

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, CATLEN, CRL, NSRC

C.........This module is required by the FileSetAPI
        USE MODFILESET

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

C.........  SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: ENAME      ! emissions inven logical name
        CHARACTER(*), INTENT (IN) :: MATTYP     ! matrix type
        CHARACTER(*), INTENT(OUT) :: MNAME      ! controls file name

C...........   EXTERNAL FUNCTIONS and their descriptions:
        CHARACTER(IODLEN3), EXTERNAL :: GETCFDSC
        INTEGER           , EXTERNAL :: GETIFDSC
        CHARACTER(16)     , EXTERNAL :: VERCHAR

C...........   LOCAL PARAMETERS
        CHARACTER(50), PARAMETER :: CVSW =
     &      '$Name SMOKEv5.0_Jun2023$' ! CVS release tag
        CHARACTER(16), PARAMETER :: PROGNAME = 'OPENCMAT' ! program name

C.........  Other local variables
        INTEGER          J, N           !  counters and indices
        INTEGER          IOS            !  i/o status
        INTEGER					 NVARS          !  number of output variables

        CHARACTER(NAMLEN3) NAMBUF   ! file name buffer
        CHARACTER(300)     MESG     ! message buffer

        CHARACTER(IODLEN3) IFDESC2, IFDESC3 ! fields 2 & 3 from inven FDESC

C***********************************************************************
C   begin body of subroutine OPENCMAT

C.........  Get header information from inventory file

        IF ( .NOT. DESCSET( ENAME,-1 ) ) THEN
            MESG = 'Could not get description of file "'
     &             // TRIM( ENAME ) // '".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IFDESC2 = GETCFDSC( FDESC3D, '/FROM/', .TRUE. )
        IFDESC3 = GETCFDSC( FDESC3D, '/VERSION/', .TRUE. )

C.........  Initialize I/O API output file headers
        CALL HDRMISS3

C.........  Set I/O API header parms that need values
        NROWS3D = NSRC
        NVARSET = NVCMULT

        FDESC3D( 1 ) = CATEGORY( 1:CATLEN ) // ' control matrix'
        FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
        WRITE( FDESC3D( 4 ), '(A,I4)' ) '/CTYPE/ ', CTYPMULT

        FDESC3D( 11 ) = '/INVEN FROM/ ' // IFDESC2
        FDESC3D( 12 ) = '/INVEN VERSION/ ' // IFDESC3

C.........  Deallocate, then allocate, output arrays
C.........  Also deallocate the number of variables per file so
C           that this will be set automatically by openset
        IF( ALLOCATED( VTYPESET ) )
     &      DEALLOCATE( VNAMESET, VTYPESET, VUNITSET, VDESCSET  )
        DEALLOCATE( VARS_PER_FILE )

        NVARS = ( 3 * NVARSET ) + NVARSET

        ALLOCATE( VNAMESET( NVARS ),
     &            VTYPESET( NVARS ),
     &            VUNITSET( NVARS ),
     &            VDESCSET( NVARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VNAMESET...VDESCSET', PROGNAME )

C.........  Set up non-speciation variables
        N = 0

        DO J = 1,NVARSET

					 N = N + 1
           VNAMESET( N )= PNAMMULT( J )
           VTYPESET( N )= M3REAL
           VUNITSET( N )= 'fraction'
           VDESCSET( N )= 'Multiplicative control factor for pollutant '
     &                   // TRIM( PNAMMULT( J ) )

           N = N + 1
           VNAMESET( N )= 'CE_'//PNAMMULT( J )
           VTYPESET( N )= M3REAL
           VUNITSET( N )= 'fraction'
           VDESCSET( N )= 'Control efficiency for pollutant '
     &                   // TRIM( PNAMMULT( J ) )

           N = N + 1
           VNAMESET( N )= 'RE_'//PNAMMULT( J )
           VTYPESET( N )= M3REAL
           VUNITSET( N )= 'fraction'
           VDESCSET( N )= 'Rule effectivness for pollutant '
     &                   // TRIM( PNAMMULT( J ) )

	         N = N + 1
           VNAMESET( N )= 'RP_'//PNAMMULT( J )
           VTYPESET( N )= M3REAL
           VUNITSET( N )= 'fraction'
           VDESCSET( N )= 'Rule penetration for pollutant '
     &                   // TRIM( PNAMMULT( J ) )

        END DO

        NVARSET = N

        MESG = 'Enter logical name for control matrix...'
        CALL M3MSG2( MESG )

C.........  Open control matrix.
C.........  Using NAMBUF is needed for HP to ensure string length consistencies
        MESG = 'I/O API ' // MATTYP // ' CONTROL MATRIX'

        NAMBUF = PROMPTSET( MESG, FSUNKN3, CRL // 'CMAT', PROGNAME )
        MNAME = NAMBUF

        RETURN

        END SUBROUTINE OPENCMAT


SUBROUTINE OPENSMAT( ENAME, SFLAG, LFLAG, NOPOL, MXSPEC,        &
                     MXTAG, EALLOUT, EAIDX, SPCNAMES, MOLUNITS, &
                     SDEV, SNAME, LNAME, SVNAMES, LVNAMES )

    !***********************************************************************
    !  subroutine body starts at line 106
    !
    !  DESCRIPTION:
    !    Open the mass-based and mole-based speciation matrices
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !    Created 2/99 by M. Houyoux 
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90", and
    !       related changes
    !***************************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !    System
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
    USE MODINFO, ONLY: CATEGORY, CATLEN, CRL, NIPPA

    !.......  This module contains the tagging arrays
    USE MODTAG, ONLY: TAGNUM, TAGNAME

    !.......  This module is required by the FileSetAPI
    USE MODFILESET

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters
    INCLUDE 'SETDECL.h90'       !  FileSetAPI function declarations

    !.......  SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: ENAME          ! emissions inven logical name
    LOGICAL     , INTENT (IN) :: SFLAG          ! true: open mass-based file
    LOGICAL     , INTENT (IN) :: LFLAG          ! true: open mole-based file
    INTEGER     , INTENT (IN) :: NOPOL          ! no. output pollutants
    INTEGER     , INTENT (IN) :: MXSPEC         ! max no. of spec per pol
    INTEGER     , INTENT (IN) :: MXTAG          ! max no. of tags per spec/pol
    CHARACTER(*), INTENT (IN) :: EALLOUT ( NIPPA )     ! output pol/emistypes
    INTEGER     , INTENT (IN) :: EAIDX   ( NIPPA )     ! index to SPCNAMES
    CHARACTER(*), INTENT (IN) :: SPCNAMES( MXSPEC, NOPOL )     ! model spec nams
    CHARACTER(*), INTENT (IN) :: MOLUNITS( MXSPEC, NOPOL )     ! mole-based unts
    INTEGER     , INTENT(OUT) :: SDEV                ! suplmt file unit no.
    CHARACTER(*), INTENT(OUT) :: SNAME               ! mass-based spec file name
    CHARACTER(*), INTENT(OUT) :: LNAME               ! mole-based spec file name
    CHARACTER(*), INTENT(OUT) :: SVNAMES( 0:MXTAG, MXSPEC, NIPPA )       ! mass out vars
    CHARACTER(*), INTENT(OUT) :: LVNAMES( 0:MXTAG, MXSPEC, NIPPA )       ! mole out vars

    !.......   EXTERNAL FUNCTIONS and their descriptions:

    CHARACTER(80), EXTERNAL :: GETCFDSC
    CHARACTER(16), EXTERNAL :: VERCHAR

    !.......   LOCAL PARAMETERS
    CHARACTER(16), PARAMETER :: PROGNAME = 'OPENSMAT'     ! program name
    CHARACTER(50), PARAMETER :: CVSW     = '$Name SMOKEv5.0_Jun2023  $'      ! CVS revision tag

    !.......  Count of species per inventory pollutant/emission type
    INTEGER    NSPEC( NIPPA )

    !.......  Other local variables
    INTEGER          I, J, K, V, T         !  counters and indices
    INTEGER          IOS                !  I/O status

    INTEGER          FMTLEN       ! length of non-blank CTMP
    INTEGER          ICNT         ! cntr for the total number of output vars
    INTEGER          NCNT         ! cntr for number of species per inv pol

    CHARACTER(12)    CTMP         ! character buffer for variable count
    CHARACTER(56)    NAMFMT       ! format for name of variables
    CHARACTER(300)   MESG         ! message buffer

    CHARACTER(NAMLEN3) NAMBUF         ! file name buffer
    CHARACTER(MXDLEN3) IFDESC2, IFDESC3     ! fields 2  3 from inven FDESC
    CHARACTER(SPNLEN3) PCODE          ! current speciation profile code
    CHARACTER(SPNLEN3) PREVCODE       ! previous speciation profile code

    !***********************************************************************
    !   begin body of subroutine OPENSMAT

    !.......  Initialize variable names
    SVNAMES = ' '      ! Array
    LVNAMES = ' '      ! Array

    !.......  Create the output variable names and count them.  Handle the special
    !    case where the total number of variables exceeds the I/O API max.
    !.......  The names of the output variables have been set up so that it will
    !    be easy to make the mass-based and the mole-based ones different.
    ICNT = 0
    DO K = 1, NIPPA

        V = EAIDX( K )

        NCNT = 0
        DO J = 1, MXSPEC

            !.......  End inner loop if species is blank
            IF( SPCNAMES( J,V ) .EQ. ' ' ) EXIT
            NCNT = NCNT + 1

            !.......  Loop through tags for this pol/spec
            DO T = 0, TAGNUM( J,V )

                !.......  Count total number of output variables
                ICNT = ICNT + 1

                !.......  Create custom format statement for building
                !    variable names. This is needed when number
                !    of variables exceeds 999, since the original
                !    format statement was I3.3 for ICNT. This actually
                !    happened for some tagging cases at EPA.
                WRITE( CTMP, '(I12)' ) ICNT
                CTMP = ADJUSTL( CTMP )
                FMTLEN = MAX( LEN( TRIM( CTMP ) ), 3 )      ! Max with 3 to replicate previous version's behavior
                WRITE( NAMFMT, '(A,I2.2,A,I2.2,A)' )    &
                       '(A4,I', FMTLEN, '.', FMTLEN, ')'

                WRITE( SVNAMES( T,J,K ), NAMFMT ) 'SVAR', ICNT
                WRITE( LVNAMES( T,J,K ), NAMFMT ) 'SVAR', ICNT

            END DO      ! end loop on tags
        END DO          ! end loop on species

        NSPEC( K ) = NCNT

    END DO

    !.......  Set up file header(s) for opening I/O API output(s). Base this on
    !    inventory header...

    !.......  Get header information from inventory file
    IF ( .NOT. DESCSET( ENAME,-1 ) ) THEN
        MESG = 'Could not get description of file "'    &
               // TRIM( ENAME ) // '".'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IFDESC2 = GETCFDSC( FDESC3D, '/FROM/', .TRUE. )
    IFDESC3 = GETCFDSC( FDESC3D, '/VERSION/', .TRUE. )

    FDESC3D = ' '       ! array

    FDESC3D( 1 ) = CATEGORY( 1:CATLEN ) // ' speciation matrix'
    FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
    FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )

    FDESC3D( 11 ) = '/INVEN FROM/ ' // IFDESC2
    FDESC3D( 12 ) = '/INVEN VERSION/ ' // IFDESC3

    !.......  Resset output number of variables
    NVARSET = ICNT

    !.......  Deallocate, then allocate, output arrays
    DEALLOCATE( VNAMESET, VTYPESET, VUNITSET, VDESCSET )
    ALLOCATE( VNAMESET( NVARSET ),    &
              VTYPESET( NVARSET ),    &
              VUNITSET( NVARSET ),    &
              VDESCSET( NVARSET ), STAT=IOS )
    CALL CHECKMEM( IOS, 'VNAMESET...VDESCSET', PROGNAME )

    !.......  Also deallocate the number of variables per file so
    !    that this will be set automatically by openset
    DEALLOCATE( VARS_PER_FILE )

    !.......  Set up variable descriptions that will be used to indicate the
    !    inventory pollutant and model species names

    VTYPESET = 0        ! array initialization
    VNAMESET = ' '      ! array initialization
    VUNITSET = ' '      ! array initialization
    VDESCSET = ' '      ! array initialization

    I = 0
    DO K = 1, NIPPA

        V = EAIDX( K )

        DO J = 1, NSPEC( K )
        DO T = 0, TAGNUM( J,V )
            I = I + 1
            VDESCSET( I ) = TRIM( EALLOUT( K )// SPJOIN//    &
                            TRIM( SPCNAMES( J,V ) ) // TAGNAME( T,J,V ) )
            VTYPESET( I ) = M3REAL
        END DO
        END DO
    END DO

    !.......  Set up variables specifically for mass-based file, and open it
    IF( SFLAG ) THEN

        FDESC3D( 4 ) = '/SMATTYPE/ ' // MASSSTR

        I = 0
        DO K = 1, NIPPA

            V = EAIDX( K )

            DO J = 1, NSPEC( K )
            DO T = 0, TAGNUM( J,V )
                I = I + 1
                VNAMESET( I ) = SVNAMES( T,J,K )
                VUNITSET( I ) = SMASUNIT
            END DO
            END DO
        END DO

        !.......  Open with NAMBUF for HP
        NAMBUF = PROMPTSET( 'Enter logical name for MASS-BASED SPECIATION MATRIX',    &
                            FSUNKN3, CRL // 'SMAT_S', PROGNAME )
        SNAME = NAMBUF

    END IF

    !.......  Set up variables specifically for mole-based file, and open it
    IF( LFLAG ) THEN

        FDESC3D( 4 ) = '/SMATTYPE/ ' // MOLESTR

        I = 0
        DO K = 1, NIPPA

            V = EAIDX( K )

            DO J = 1, NSPEC( K )
            DO T = 0, TAGNUM( J,V )
                I = I + 1
                VNAMESET( I ) = LVNAMES( T,J,K )
                VUNITSET( I ) = MOLUNITS( J,V )
            END DO
            END DO
        END DO

        !.......  Open with NAMBUF for HP
        NAMBUF = PROMPTSET( 'Enter logical name for MOLE-BASED SPECIATION MATRIX',    &
                            FSUNKN3, CRL // 'SMAT_L', PROGNAME )
        LNAME = NAMBUF

    END IF

    !.......  Open supplemental speciation file
    MESG = 'Enter logical name for the SPECIATION SUPPLEMENTAL '//    &
           'file'
    SDEV = PROMPTFFILE( MESG, .FALSE., .TRUE., CRL // 'SSUP', PROGNAME )


    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE OPENSMAT

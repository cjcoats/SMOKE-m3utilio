
SUBROUTINE OPENREPIN( ENAME, ANAME, CUNAME, GNAME, LNAME,           &
                      PRNAME, SLNAME, SSNAME, TNAME, RDEV,          &
                      SDEV, GDEV, PDEV, TDEV, EDEV, YDEV, NDEV,     &
                      NIDEV, NPDEV, ADEV, NMDEV, NNDEV, NODEV )

    !***********************************************************************
    !  subroutine OPENREPIN body starts at line
    !
    !  DESCRIPTION:
    !      The purpose of this subroutine is to open all of the necessary
    !      files for the Smkreport routine and set the episode information
    !      for the calling program.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 7/2000 by M. Houyoux
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
    !****************************************************************************
    USE M3UTILIO

    !.......  MODULES for public variables
    !.......  This module contains Smkreport-specific settings
    USE MODREPRT, ONLY: AFLAG, TSTEP, NSTEPS, DATAMISS, MXINDAT,    &
                        TSFLAG, TFLAG, GFLAG, GSFLAG, SLFLAG,       &
                        PSFLAG, SSFLAG, PRFLAG, PRRPTFLG, NMATX,    &
                        PRBYR, PRPYR, PYEAR, CHKPFX, CUFLAG,        &
                        LFLAG, EMLAYS, VFLAG, YFLAG, NFLAG,         &
                        ASCREC, ASCDATA, STIME, SDATE, ETIME,       &
                        EDATE, TZONE, NIFLAG, NMFLAG, NNFLAG,       &
                        NOFLAG, SDFLAG

    !.......  This module contains the temporal profile tables
    USE MODTMPRL, ONLY: NTPDAT, TPNAME, TPUNIT, TPDESC

    !.......  This module contains report arrays for each output bin
    USE MODREPBN, ONLY: NSVARS

    !.......  This module contains the control packet data and control matrices
    USE MODCNTRL, ONLY: NVPROJ, PNAMPROJ, NVCMULT, PNAMMULT

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, CRL, CATDESC, NIPPA, NCHARS, JSCC, &
                       JSTACK, PLTIDX, MXCHRS, NSRC, INVPIDX, BYEAR,&
                       EANAM, EAUNIT, SC_BEGP, SC_ENDP

    !.......  This module is required for the FileSetAPI
    USE MODFILESET

    IMPLICIT NONE

    !.......  INCLUDES:

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters
    INCLUDE 'SETDECL.h90'       !  FileSetAPI function declarations

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT(OUT) :: ENAME      ! name for I/O API inven input
    CHARACTER(*), INTENT(OUT) :: ANAME      ! name for ASCII inven input
    CHARACTER(*), INTENT(OUT) :: CUNAME     ! multiplicative control matrix name
    CHARACTER(*), INTENT(OUT) :: GNAME      ! gridding matrix name
    CHARACTER(*), INTENT(OUT) :: LNAME      ! layer fractions file name
    CHARACTER(*), INTENT(OUT) :: PRNAME     ! projection matrix name
    CHARACTER(*), INTENT(OUT) :: SLNAME     ! speciation matrix name
    CHARACTER(*), INTENT(OUT) :: SSNAME     ! speciation matrix name
    CHARACTER(*), INTENT(OUT) :: TNAME      ! hourly emissions file
    INTEGER     , INTENT(OUT) :: RDEV(3)    ! control report files
    INTEGER     , INTENT(OUT) :: SDEV       ! unit no.: ASCII inven file
    INTEGER     , INTENT(OUT) :: GDEV       ! gridding supplemental file
    INTEGER     , INTENT(OUT) :: PDEV       ! speciation supplemental file
    INTEGER     , INTENT(OUT) :: TDEV       ! temporal supplemental file
    INTEGER     , INTENT(OUT) :: EDEV       ! unit no.: elevated ID file (PELV)
    INTEGER     , INTENT(OUT) :: YDEV       ! unit no.: cy/st/co file
    INTEGER     , INTENT(OUT) :: NDEV       ! unit no.: SCC descriptions
    INTEGER     , INTENT(OUT) :: NIDEV      ! unit no.: SIC descriptions
    INTEGER     , INTENT(OUT) :: NPDEV      ! unit no.: GSPRO descriptions
    INTEGER     , INTENT(OUT) :: NMDEV      ! unit no.: MACT descriptions
    INTEGER     , INTENT(OUT) :: NNDEV      ! unit no.: NAICS descriptions
    INTEGER     , INTENT(OUT) :: NODEV      ! unit no.: ORIS descriptions
    INTEGER     , INTENT(OUT) :: ADEV       ! unit no.: ASCII elevated file

    !....... EXTERNAL FUNCTIONS and their descriptions:

    CHARACTER(MXDLEN3),EXTERNAL :: GETCFDSC
    INTEGER           ,EXTERNAL :: GETIFDSC
    LOGICAL           ,EXTERNAL :: USEEXPGEO

    !.......  Temporary array for speciation variable names
    CHARACTER(MXDLEN3), ALLOCATABLE :: SLVNAMS( : )

    !.......  Local units and logical file names
    INTEGER         IDEV          ! tmp unit number if ENAME is map file
    INTEGER      :: MDEV = 0      ! unit no. emission processes file

    CHARACTER(16)   INAME         ! tmp name for inven file of unknown fmt

    !.......  Other local variables

    INTEGER         I, J, L, L1, L2, N, V           ! counters and indices

    INTEGER         IOS               ! tmp I/O status
    INTEGER         ISD               ! start time of ASCII elevated file
    INTEGER         IED               ! end time of ASCII elevated file
    INTEGER         TSTEP_T           ! unused time step from environment

    LOGICAL      :: EFLAG = .FALSE.      ! true: error found
    LOGICAL      :: TIMEFLAG = .FALSE.      ! true: time info already init

    CHARACTER(16)   NAMBUF           ! tmp file name buffer
    CHARACTER(16)   UNITS            ! units of ASCII elevated file
    CHARACTER(256)  MESG             ! message buffer
    CHARACTER(300)  LINE             ! tmp line buffer

    CHARACTER(NAMLEN3) GRDENV          ! gridded output units from envrmt

    CHARACTER(16) :: PROGNAME = 'OPENREPIN'     ! program name

    !***********************************************************************
    !   begin body of subroutine OPENREPIN

    IF( .NOT. AFLAG ) THEN
        !.......  Get inventory file names given source category
        CALL GETINAME( CATEGORY, ENAME, ANAME )

        !.......  Prompt for and open inventory file
        INAME = ENAME
        MESG = 'Enter logical name for the MAP INVENTORY file'
        IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INAME, PROGNAME )

        !.......  Open and read map file
        CALL RDINVMAP( INAME, IDEV, ENAME, ANAME, SDEV )

        !.......  Store source-category-specific header information,
        !           including the inventory pollutants in the file (if any).  Note that
        !           the I/O API header info is passed by include file and the
        !           results are stored in module MODINFO.
        CALL GETSINFO( ENAME )

        !.......  Store non-category-specific header information
        TSTEP  = 000000
        NSTEPS = 1

        !.......  Reset the maximum input data if any reports did not select
        !           specific data values.  MXINDAT might get larger than needed.
        IF( DATAMISS ) THEN
            MXINDAT = MAX( NIPPA, MXINDAT )
        END IF

        !.......  Determine the year and projection status of the inventory
        !           CALL CHECK_INVYEAR( ENAME, APRJFLAG, FDESC3D )

    ELSE
        NCHARS = 3
        JSCC = 0
        JSTACK = 3

        ALLOCATE( SC_BEGP( NCHARS ),    &
              SC_ENDP( NCHARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SC_BEGP,SC_ENDP', PROGNAME )

        PLTIDX = 2
        MXCHRS = MXPTCHR3

        DO I = 1, NCHARS
            SC_BEGP( I ) = PTBEGL3( I )
            SC_ENDP( I ) = PTENDL3( I )
        END DO

    END IF

    !.......  For temporal inputs, prompt for hourly file
    IF( TFLAG ) THEN

        MESG = 'Enter logical name for the HOURLY  EMISSIONS file'
        TNAME = PROMPTSET( MESG, FSREAD3, CRL//'TMP', PROGNAME )

        !.......  Set parameters and pollutants from hourly file
        CALL RETRIEVE_SET_HEADER( TNAME )
        CALL CHKSRCNO( CATDESC, TNAME, NROWS3D, NSRC, EFLAG )
        CALL UPDATE_TIME_INFO( TNAME )

        !.......  Determine average day emissions status from hourly file
        INVPIDX = GETIFDSC( FDESC3D, '/AVERAGE DAY/', .FALSE. )
        IF( INVPIDX .EQ. 1 ) THEN
            MESG = 'NOTE: Average day emissions in hourly emissions file'
            CALL M3MSG2( MESG )
        END IF

        !.......  Store variable number, names, and units from the hourly
        !               emissions fileTPNAME
        NTPDAT = NVARSET
        ALLOCATE( TPNAME( NTPDAT ),     &
                  TPUNIT( NTPDAT ),     &
                  TPDESC( NTPDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TPNAME...TPDESC', PROGNAME )

        TPNAME = VNAMESET( 1:NTPDAT )          ! array
        TPUNIT = VUNITSET( 1:NTPDAT )          ! array
        TPDESC = VDESCSET( 1:NTPDAT )          ! array

        !.......  Determine the year and projection status of the hourly
        !           CALL CHECK_INVYEAR( TNAME, PRJFLAG, FDESC3D )

    END IF

    IF( TSFLAG ) THEN

        MESG = 'Enter logical name for the TEMPORAL SUPPLEMENTAL file'
        TDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., CRL//'TSUP', PROGNAME )

    END IF

    !.......  Open gridding matrix and compare number of sources
    IF( GFLAG ) THEN

        GNAME = PROMPTMFILE( 'Enter logical name for the GRIDDING MATRIX',&
                             FSREAD3, CRL//'GMAT', PROGNAME )

        CALL RETRIEVE_IOAPI_HEADER( GNAME )
        CALL CHKSRCNO( CATDESC, GNAME, NTHIK3D, NSRC, EFLAG )

        !.......  Initialize grid description
        CALL CHKGRID( CATDESC, 'GMAT', 0, EFLAG )

        !.......  Store gridding matrix size
        NMATX = NCOLS3D

    END IF

    IF( GSFLAG ) THEN

        MESG = 'Enter logical name for the GRIDDING SUPPLEMENTAL '//&
               'file'
        GDEV = PROMPTFFILE( MESG, .TRUE., .TRUE.,&
                            CRL//'GSUP', PROGNAME )

    END IF

    !.......  Open mole speciation matrix, compare number of sources, store
    !           speciation variable descriptions, and store mass or moles.
    IF( SLFLAG ) THEN

        SLNAME = PROMPTSET( 'Enter logical name for the MOLE SPECIATION MATRIX',&
                            FSREAD3, CRL//'SMAT_L', PROGNAME )

        CALL RETRIEVE_SET_HEADER( SLNAME )
        CALL CHKSRCNO( CATDESC, SLNAME, NROWS3D, NSRC, EFLAG )

        NSVARS = NVARSET

        ALLOCATE( SLVNAMS( NSVARS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SLVNAMS', PROGNAME )

        SLVNAMS = VDESCSET      ! array

    END IF      ! end of mole speciation open

    IF( PSFLAG ) THEN

        MESG = 'Enter logical name for the SPECIATION SUPPLEMENTAL file'
        PDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., CRL//'SSUP', PROGNAME )

    END IF

    !.......  Open mass speciation matrix, compare number of sources, store
    !           speciation variable descriptions, and store mass or moles.
    IF( SSFLAG ) THEN

        SSNAME = PROMPTSET( 'Enter logical name for the MASS SPECIATION MATRIX',&
                           FSREAD3, CRL//'SMAT_S', PROGNAME )

        CALL RETRIEVE_SET_HEADER( SSNAME )
        CALL CHKSRCNO( CATDESC, SSNAME, NROWS3D, NSRC, EFLAG )

        !.......  Compare matrix header with mole-based, if available
        IF( SLFLAG ) THEN

            !.......  Check the number of variables
            IF( NSVARS .NE. NVARSET ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Inconsistent number '//     &
                  'of speciation variables.'// CRLF()// BLANK10//       &
                  'Mole file:', NSVARS,'; Mass file:', NVARSET
                CALL M3MSG2( MESG )

            END IF

            !.......  Check the dscriptions of variables
            DO V = 1, NSVARS

                IF( SLVNAMS( V ) .NE. VDESCSET( V ) ) THEN

                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Inconsistent '//    &
                      'variable descriptions in speciation '//      &
                      'matrices for variable', V, CRLF() //         &
                      BLANK10 //'Mole file: ', SLVNAMS( V ) //      &
                      CRLF() // BLANK10 //'Mass file: ', VDESCSET(V)
                    CALL M3MSG2( MESG )

                END IF

            END DO

        !.......  Otherwise, set number of speciation variables
        ELSE

            NSVARS  = NVARSET

        END IF

    END IF      ! end of mass speciation open

    !.......  Open projection matrix, compare number of sources,
    !               and store projection variable names.
    IF( PRFLAG ) THEN

        MESG = 'Enter logical name for the PROJECTION MATRIX'
        PRNAME = PROMPTSET( MESG, FSREAD3, CRL//'PMAT', PROGNAME )

        CALL RETRIEVE_SET_HEADER( PRNAME )
        CALL CHKSRCNO( CATDESC, PRNAME, NROWS3D, NSRC, EFLAG )
        NVPROJ = NVARS3D

        !.......  Set allocation size depending on whether report is also
        !             read in.
        I = NVPROJ
        IF( PRRPTFLG ) I = 2 * I

        !.......  Allocate memory for variables (including possible
        !             other variables for report check) and store names
        ALLOCATE( PNAMPROJ( I ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PNAMPROJ', PROGNAME )
        PNAMPROJ = ' '      ! array

        CALL STORE_VNAMES( 1, 1, NVPROJ, PNAMPROJ )

        IF( NVPROJ .NE. 1 ) THEN
            MESG = 'INTERNAL ERROR: Smkreport is not set up to ' // &
                   'support more than 1 variable in ' //            &
                   CRLF() // BLANK10 // 'the projection matrix.'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        PRBYR = GETIFDSC( FDESC3D, '/BASE YEAR/', .TRUE. )
        PRPYR = GETIFDSC( FDESC3D, '/PROJECTED YEAR/', .TRUE. )

        IF( PYEAR .LE. 0 ) THEN
            IF ( PRBYR .NE. BYEAR ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Inventory base year',    BYEAR,    &
                      'is not consistent with projection base year', PRBYR
                CALL M3MSG2( MESG )

            ELSE
                WRITE( MESG,94010 ) 'NOTE: Base year', BYEAR,           &
                       'is consistent between the inventory and '//     &
                       CRLF() // BLANK10 // 'projection matrix.'
                CALL M3MSG2( MESG )
            END IF

        ELSE

            IF ( PRBYR .NE. PYEAR ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Inventory projected'//      &
                       ' data year', PYEAR, 'is not consistent '//      &
                       'with projection base year', PRBYR
                CALL M3MSG2( MESG )

            ELSE
                WRITE( MESG,94010 ) 'WARNING: Inventory projected'//    &
                       ' year', PYEAR, 'is being projected to',         &
                       PRPYR, 'by projection matrix.'
                CALL M3MSG2( MESG )

            END IF
        END IF
    END IF      ! end of multiplicative control open

    !.......  Open projections report
    IF( PRRPTFLG ) THEN

        MESG = 'Enter logical name for input PROJECTION REPORT from Cntlmat'
        RDEV(1) = PROMPTFFILE( MESG, .TRUE., .TRUE., CRL // 'PROJRPT', PROGNAME )

        DO V = NVPROJ+1, 2*NVPROJ
            J = V - NVPROJ
            PNAMPROJ( V ) = CHKPFX // PNAMPROJ( J )
        END DO

    END IF

    !.......  Open multiplicative control matrix, compare number of sources,
    !               and store control variable names.
    IF( CUFLAG ) THEN

        MESG = 'Enter logical name for the ' //&
               'MULTIPLICATIVE CONTROL MATRIX'
        CUNAME = PROMPTSET( MESG, FSREAD3, CRL//'CMAT', PROGNAME )

        CALL RETRIEVE_SET_HEADER( CUNAME )
        CALL CHKSRCNO( CATDESC, CUNAME, NROWS3D, NSRC, EFLAG )
        NVCMULT = NVARS3D
        ALLOCATE( PNAMMULT( NVCMULT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PNAMMULT', PROGNAME )
        CALL STORE_VNAMES( 1, 1, NVCMULT, PNAMMULT )

    END IF      ! end of multiplicative control open

    !.......  Open additive control matrix, compare number of sources,
    !               and store control variable names.
    !        IF( CAFLAG ) THEN
    !            MESG= 'INTERNAL ERROR: Area additive controls not ' //
    !                 'yet implemented in ' // PROGNAME
    !            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    !        END IF      ! end of additive control open

    !.......  Open reactivity control matrix, compare number of sources, and
    !               store control variable descriptions, and store mass or moles.
    !        IF( CRFLAG ) THEN
    !            RNAME = PROMPTMFILE(
    !                'Enter logical name for the REACTIVITY MATRIX',
    !                FSREAD3, CRL//'RMAT', PROGNAME )

    !            CALL RETRIEVE_IOAPI_HEADER( RNAME )
    !            CALL CHKSRCNO( CATDESC, RNAME, NTHIK3D, NSRC, EFLAG )
    !            NRMATV = NVARS3D
    !            NSREAC = NROWS3D
    !            ALLOCATE( RVDESC( NRMATV ), STAT=IOS )
    !            CALL CHECKMEM( IOS, 'RVDESC', PROGNAME )
    !            CALL STORE_VDESCS( 1, 1, NRMATV, RVDESC )

    !.......  Retrieve the number of speciation factors
    !            RNMSPC = GETIFDSC( FDESC3D, '/SPECIES VARS/', .TRUE. )

    !.......  Check the year and projection year of the matrix
    !           CALL CHECK_INVYEAR( ARNAME, APRJFLAG, FDESC3D )

    !        END IF      ! end of reactivity control open

    !.......  Open layer fractions file, compare number of sources, check
    !           met information, and store the vertical coordinates info
    IF( LFLAG ) THEN

        MESG= 'Enter logical name for the POINT LAYER FRACTIONS MATRIX'
        LNAME = PROMPTMFILE( MESG, FSREAD3, 'PLAY', PROGNAME )

        CALL RETRIEVE_IOAPI_HEADER( LNAME )
        CALL CHKSRCNO( CATDESC, LNAME, NROWS3D, NSRC, EFLAG )
        CALL UPDATE_TIME_INFO( LNAME )
        EMLAYS = NLAYS3D

    END IF      ! End of layer fractions open

    !.......  Open elevated/low-level
    IF( VFLAG ) THEN

    !.......  Open elevated/plume-in-grid file
        MESG = 'Enter logical name for the ELEVATED/PING file'
        EDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'PELV', PROGNAME      )

    END IF

    !.......  Get country, state, and county names, if needed
    IF( YFLAG .AND. .NOT. USEEXPGEO() ) THEN

        MESG = 'Enter logical name for COUNTRY, STATE, AND ' //&
               'COUNTY file'
        YDEV = PROMPTFFILE( MESG,.TRUE.,.TRUE.,'COSTCY',PROGNAME )

    END IF

    !.......  Get SCC descriptions, if needed
    IF( NFLAG ) THEN

        MESG = 'Enter logical name for SCC DESCRIPTIONS'
        NDEV = PROMPTFFILE( MESG,.TRUE.,.TRUE.,'SCCDESC',PROGNAME )

    END IF

    !.......  Get SIC descriptions, if needed
    IF( NIFLAG ) THEN

        MESG = 'Enter logical name for SIC DESCRIPTIONS'
        NIDEV = PROMPTFFILE( MESG,.TRUE.,.TRUE.,'SICDESC',PROGNAME )

    END IF

    !.......  Get GSPRO descriptions, if needed
    IF( SDFLAG ) THEN

        MESG = 'Enter logical name for GSPRO DESCRIPTIONS'
        NPDEV = PROMPTFFILE( MESG,.TRUE.,.TRUE.,'GSPRODESC',PROGNAME )

    END IF

    !.......  Get MACT descriptions, if needed
    IF( NMFLAG ) THEN

        MESG = 'Enter logical name for MACT DESCRIPTIONS'
        NMDEV = PROMPTFFILE( MESG,.TRUE.,.TRUE.,'MACTDESC',PROGNAME )

    END IF

    !.......  Get NAICS descriptions, if needed
    IF( NNFLAG ) THEN

        MESG = 'Enter logical name for NAICS DESCRIPTIONS'
        NNDEV = PROMPTFFILE(&
                    MESG,.TRUE.,.TRUE.,'NAICSDESC',PROGNAME )

    END IF

    !.......  Get ORIS descriptions, if needed
    IF( NOFLAG ) THEN

        MESG = 'Enter logical name for ORIS DESCRIPTIONS'
        NODEV = PROMPTFFILE( MESG,.TRUE.,.TRUE.,'ORISDESC',PROGNAME )

    END IF

    !.......  Open ASCII elevation file output by SMKMERGE, if needed
    IF( AFLAG ) THEN

        CALL ENVSTR( 'MRG_GRDOUT_UNIT', ' ', ' ', GRDENV, IOS)

        IF( GRDENV( 1:1 ) .EQ. 'm' ) THEN

            ADEV = PROMPTFFILE( 'Enter name for ASCII ELEVATED SOURCES file',   &
                                .TRUE., .TRUE., 'ELEVTS_L', PROGNAME )

        ELSE

            ADEV = PROMPTFFILE( 'Enter name for ASCII ELEVATED SOURCES file',   &
                                .TRUE., .TRUE., 'ELEVTS_S', PROGNAME )

        END IF

        !.......  Read ASCII elevated file
        MESG = 'Reading ASCII elevated file...'
        CALL M3MSG2( MESG )

        ASCREC = 0

        !......  Read in units
        ASCREC = ASCREC + 1
        READ( ADEV, '(10X,A)') UNITS

        !......  Skip header lines
        DO I = 1, 2
            ASCREC = ASCREC + 1
            READ( ADEV, '(A)' ) LINE
        END DO

        !.......  Get number of species and point stacks from file
        ASCREC = ASCREC + 1
        READ( ADEV, '(I10,10X,I10)' ) ASCDATA, NSRC

        NIPPA = ASCDATA
        NSVARS = NIPPA
        MXINDAT = MAX( NIPPA, MXINDAT )

        ALLOCATE( EANAM( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EANAM', PROGNAME )
        ALLOCATE( EAUNIT( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EAUNIT', PROGNAME )
        EANAM = ''
        EAUNIT = TRIM( UNITS )

        DO I = 1, 4
            ASCREC = ASCREC + 1
            READ( ADEV, '(A)' ) LINE
        END DO

        !......  Read in list of species
        DO I = 1, ASCDATA
            ASCREC = ASCREC + 1
            READ( ADEV, '(A)' ) LINE
            EANAM( I ) = TRIM( LINE )
        END DO

        !......  Read in start and end dates and times
        ASCREC = ASCREC + 1
        READ( ADEV, '(4I10)' ) ISD, STIME, IED, ETIME
        SDATE = ISD + 1900000
        EDATE = IED + 1900000
        BYEAR = INT( SDATE/1000 )
        TSTEP = 10000
        NSTEPS = 1 + SECSDIFF( SDATE, STIME, EDATE, ETIME )/3600

    END IF

    !.......  If there were any errors inputing files or while comparing
    !           with one another, then abort
    IF( EFLAG ) THEN

        MESG = 'Problems opening input files. See ERROR(S) above.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    !.......  If we are using temporalized emissions, then update date/time and
    !           duration using environment variable settings, then prompt.
    IF( TFLAG ) THEN

        !.......  Write explanation
        MESG = 'For time-based reports, enter the starting date, '//    &
               'starting time, and output' // CRLF() // BLANK10 //      &
               'duration.  Defaults are set based on the input ' //     &
               'files.'
        CALL M3MSG2( MESG )

        !.......  Subselect dates and times
        CALL GETM3EPI( TZONE, SDATE, STIME, TSTEP_T, NSTEPS )
        EDATE = SDATE
        ETIME = STIME
        CALL NEXTIME( EDATE, ETIME, ( NSTEPS-1 ) * TSTEP )

    END IF       !  if have temporalized inputs and outputs

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats.......94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

    !******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

    !.......  This internal subprogram tries to retrieve the I/O API header
    !               and aborts if it was not successful
    SUBROUTINE RETRIEVE_IOAPI_HEADER( FILNAM )

        !.......  Subprogram arguments
        CHARACTER(*), INTENT (IN) :: FILNAM

        !----------------------------------------------------------------------

        IF ( .NOT. DESC3( FILNAM ) ) THEN

            MESG = 'Could not get description of file "' //&
                   TRIM( FILNAM ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ENDIF

    END SUBROUTINE RETRIEVE_IOAPI_HEADER

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    !.......  This subprogram tries to retrieve the description for a file
    !               set and aborts if it was not successful
    SUBROUTINE RETRIEVE_SET_HEADER( FILNAM )

        INCLUDE 'SETDECL.h90'       !  FileSetAPI function declarations

        !.......  Subprogram arguments
        CHARACTER(*) FILNAM

        !----------------------------------------------------------------------

        IF ( .NOT. DESCSET( FILNAM, ALLFILES ) ) THEN

            MESG = 'Could not get description of file set "' //&
                   TRIM( FILNAM ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ENDIF

    END SUBROUTINE RETRIEVE_SET_HEADER

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------

    !.......  This subprogram updates the time (episode) information
    !               and compares to the existing information, if it has been
    !               previously set.
    SUBROUTINE UPDATE_TIME_INFO( FILNAM )

        !.......  Subprogram arguments
        CHARACTER(*) FILNAM

        !.......  Local variables
        INTEGER ISECS           ! number of seconds different between dates/times
        INTEGER ED              ! tmp ending date
        INTEGER ET              ! tmp ending time
        INTEGER LOCZONE         ! tmp time zone

        INTEGER, SAVE :: EDATE = 0            ! Ending date
        INTEGER, SAVE :: ETIME = 0            ! Ending time

        LOGICAL, SAVE :: IFLAG = .FALSE.          ! true: episode settings init
        LOGICAL, SAVE :: ZFLAG = .FALSE.          ! true: time zone init

        CHARACTER(300)   MESG                     ! message buffer

        !----------------------------------------------------------------------

        !.......  If time information has already been initialized...
        IF( IFLAG ) THEN
            ISECS = SECSDIFF( SDATE, STIME, SDATE3D, STIME3D )

            IF( ISECS .GT. 0 ) THEN          ! SDATE3D/STIME3D are later
                SDATE = SDATE3D
                STIME = STIME3D
            END IF

            ED = SDATE3D
            ET = STIME3D
            CALL NEXTIME( ED, ET, ( MXREC3D-1 ) * TSTEP3D )

            ISECS = SECSDIFF( EDATE, ETIME, ED, ET )

            IF( ISECS .LT. 0 ) THEN          ! ED/ET are earlier
                EDATE = ED
                ETIME = ET
            END IF

            NSTEPS = 1+ SECSDIFF( SDATE, STIME, EDATE, ETIME )/ 3600

            IF( TFLAG .AND. NSTEPS .LE. 0 ) THEN
                MESG = 'Because of file ' // FILNAM //  &
                       ', dates and times do not overlap at all        !'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

            !.......  Check time step
            IF( TSTEP3D .NE. TSTEP ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Time step is not one hour in ' // FILNAM // ' file        !'
                CALL M3MSG2( MESG )
            END IF

        !.......  If time information needs to be initialized...
        ELSE
            SDATE  = SDATE3D
            STIME  = STIME3D
            NSTEPS = MXREC3D
            TSTEP  = TSTEP3D

            EDATE  = SDATE
            ETIME  = STIME
            CALL NEXTIME( EDATE, ETIME, ( NSTEPS-1 ) * TSTEP )

            IFLAG = .TRUE.

        END IF

        !.......  Make sure that time step is one hour
        IF( TSTEP .NE. 10000 ) THEN

            EFLAG = .TRUE.
            MESG = 'ERROR: Time step is not one hour in ' // FILNAM // ' file        !'
            CALL M3MSG2( MESG )

        END IF

        !.......  Use layers to screen for non-layer-fractions files
        !.......  If not layer-fractions file, retrieve and compare time zone
        IF( NLAYS3D .LE. 1 ) THEN

            LOCZONE = GETIFDSC( FDESC3D, '/TZONE/', .TRUE. )

            IF( ZFLAG .AND. LOCZONE .NE. TZONE ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )                                 &
                   'Time zone ', LOCZONE, 'in ' // FILNAM //        &
                   ' hourly emissions file is not consistent ' //   &
                   'with initialized value of', TZONE
                CALL M3MSG2( MESG )

            ELSE IF( .NOT. ZFLAG ) THEN
                ZFLAG = .TRUE.
                TZONE = LOCZONE

                MESG = 'NOTE: Time zone initialized using ' //      &
                       FILNAM // ' hourly emissions file.'

                CALL M3MSG2( MESG )
            END IF
        END IF

        !------------------  FORMAT  STATEMENTS   -----------------------------

        !.......   Internal buffering formats.......94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

    END SUBROUTINE UPDATE_TIME_INFO

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    !.......  This subprogram initializes and checks the inventory year
    !               of the emissions and the projection status
    SUBROUTINE CHECK_INVYEAR( FNAME, PRJFLAG, IODESC )

        !.......  Subprogram arguments
        CHARACTER(*), INTENT (IN)     :: FNAME
        LOGICAL     , INTENT (IN OUT) :: PRJFLAG
        CHARACTER(*), INTENT (IN)     :: IODESC( * )

        !.......  Local variables
        INTEGER           YY              ! tmp year

        INTEGER , SAVE :: BYEAR           ! base year
        INTEGER , SAVE :: FLEN            ! name length of savnam

        LOGICAL           STRICT          ! flag for strict checks or not

        CHARACTER(20)     BUFFER          ! program name buffer
        CHARACTER(NAMLEN3), SAVE :: SAVNAM          ! name of file used to init

        !----------------------------------------------------------------------

        STRICT = .TRUE.

        !.......  First determine whether to abort when projected year does not
        !               match.  This is used for reactivity matrices, which will
        !               always have a projection year, even if the inventory isn't
        !               projected.
        IF( .NOT. PRJFLAG ) THEN
            BUFFER = GETCFDSC( FDESC3D, '/FROM/', .FALSE. )
            IF( BUFFER .EQ. 'OPENRMAT' ) STRICT = .FALSE.
        END IF

        !.......  If time information has already been initialized...
        IF( TIMEFLAG ) THEN

            YY = GETIFDSC( IODESC, '/PROJECTED YEAR/', .FALSE. )
            IF( YY .LE. 0 ) THEN

                YY = GETIFDSC( IODESC, '/BASE YEAR/', .FALSE. )
                IF( YY .NE. BYEAR ) THEN
                    WRITE( MESG,94010 )                                 &
                          'Base year of ' // FNAME // ' file:', YY,     &
                          CRLF() // BLANK10 //                          &
                          ', does not equal emissions year of ' //      &
                          SAVNAM( 1:FLEN ) // ' file:', BYEAR
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

            ELSE IF ( STRICT            .AND.&
                      YY     .GT. 0     .AND.&
                      YY     .NE. PYEAR      ) THEN

                WRITE( MESG,94010 )                                     &
                      'Projected year of ' // FNAME // ' file:', YY,    &
                      CRLF() // BLANK10 //                              &
                      ', does not equal emissions year of ' //          &
                      SAVNAM( 1:FLEN ) // ' file:', PYEAR
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

        !.......  If year information needs to be initialized...
        ELSE

            BYEAR = GETIFDSC( IODESC, '/BASE YEAR/', .FALSE. )
            PYEAR = GETIFDSC( IODESC, '/PROJECTED YEAR/', .FALSE. )

            IF( YY .GT. 0 ) THEN
                PRJFLAG = .TRUE.
            ELSE
                PYEAR = BYEAR
            END IF

            SAVNAM = FNAME
            FLEN   = LEN_TRIM( SAVNAM )
            TIMEFLAG  = .TRUE.

        END IF

        !------------------  FORMAT  STATEMENTS   -----------------------------

        !.......   Internal buffering formats.......94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

    END SUBROUTINE CHECK_INVYEAR

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    !.......  This subprogram stores I/O API NetCDF variable names into
    !               a local array based on indices in subprogram call.
    SUBROUTINE STORE_VNAMES( ISTART, INCRMT, NNAM, NAMES )

        !.......  Subprogram arguments
        INTEGER      ISTART                ! starting position in VNAMES of names
        INTEGER      INCRMT                ! increment of VNAMES for names
        INTEGER      NNAM                  ! number of names
        CHARACTER(*) NAMES( NNAM )         ! stored variable names

        !.......  Local variables
        INTEGER  I, J

        !----------------------------------------------------------------------

        J = ISTART
        DO I = 1, NNAM

            NAMES( I ) = VNAMESET( J )
            J = J + INCRMT

        END DO

    END SUBROUTINE STORE_VNAMES

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    !.......  This subprogram stores I/O API NetCDF variable descriptions into
    !               a local array based on indices in subprogram call.
    SUBROUTINE STORE_VDESCS( ISTART, INCRMT, NDESC, DESCS )

        !.......  Subprogram arguments
        INTEGER      ISTART                ! starting position in VDESCS of names
        INTEGER      INCRMT                ! increment of VDESCS for names
        INTEGER      NDESC                 ! number of descriptions
        CHARACTER(*) DESCS( NDESC )        ! stored variable descriptions

        !.......  Local variables
        INTEGER  I, J, L

        !----------------------------------------------------------------------

        DESCS = ' '

        J = ISTART
        DO I = 1, NDESC

            L = LEN_TRIM( VDESCSET( J ) )
            DESCS( I ) = VDESCSET( J )( 1:L )
            J = J + INCRMT

        END DO

    END SUBROUTINE STORE_VDESCS

END SUBROUTINE OPENREPIN

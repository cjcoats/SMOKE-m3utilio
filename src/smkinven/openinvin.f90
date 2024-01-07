
SUBROUTINE OPENINVIN( CATEGORY, ADEV, DDEV, HDEV, RDEV, SDEV,   &
                      PDEV, ZDEV, CDEV, ODEV, UDEV,             &
                      YDEV, ENAME, INNAME, IDNAME, IHNAME )

    !***********************************************************************
    !  subroutine body starts at line 119
    !
    !  DESCRIPTION:
    !      This subroutine opens the appropriate input files for inventory import
    !      given the source category.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !      Subroutines: I/O API subroutines
    !
    !  REVISION  HISTORY:
    !       Created 4/99 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !**************************************************************************
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
    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: NMAP, MAPNAM, MAPFIL, NCOMP, VAR_FORMULA, NETCDFUNIT

    !.......  This module contains the lists of unique inventory information
    USE MODLISTS, ONLY: MEDSFLAG, NCDFLAG, APIFLAG

    !...... This module contains the cross-reference tables
    USE MODXREF, ONLY: PROC_HAPS

    !.......  This module contains data for day- and hour-specific data
    USE MODDAYHR, ONLY: DAYINVFLAG, HRLINVFLAG, FF10INVFLAG

    !.......  This module is for mobile-specific data
    USE MODMOBIL, ONLY: SCCMAPFLAG, SCCMAPLIST, EXCLSCCFLAG

    IMPLICIT NONE

    !.......   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: CATEGORY      ! source category
    INTEGER     , INTENT(OUT) :: ADEV          ! unit number for inven file
    INTEGER     , INTENT(OUT) :: DDEV          ! unit no. for day-specific file
    INTEGER     , INTENT(OUT) :: HDEV          ! unit no. for hr-specific file
    INTEGER     , INTENT(OUT) :: RDEV          ! unit no. for stack replacements
    INTEGER     , INTENT(OUT) :: SDEV          ! unit no. for optional inven in
    INTEGER     , INTENT(OUT) :: PDEV          ! unit no. for inven data table
    INTEGER     , INTENT(OUT) :: ZDEV          ! unit no. for time zones
    INTEGER     , INTENT(OUT) :: CDEV          ! unit no. for SCCs description
    INTEGER     , INTENT(OUT) :: ODEV          ! unit no. for ORIS description
    INTEGER     , INTENT(OUT) :: UDEV          ! unit no. for non-HAP inclusions/exclusions
    INTEGER     , INTENT(OUT) :: YDEV          ! unit no. for area-to-point
    CHARACTER(*), INTENT(OUT) :: ENAME         ! optional I/O API inven input
    CHARACTER(*), INTENT(OUT) :: INNAME        ! average inventory name
    CHARACTER(*), INTENT(OUT) :: IDNAME        ! day-specific inventory
    CHARACTER(*), INTENT(OUT) :: IHNAME        ! hour-specific inventory name

    !.......   EXTERNAL FUNCTIONS and their descriptionsNRAWIN

    LOGICAL, EXTERNAL :: USEEXPGEO

    !.......   Other local variables
    INTEGER       MDEV       ! tmp unit number if SCC reference map file is needed
    INTEGER       IDEV       ! tmp unit number if ENAME is map file
    INTEGER       IOS        ! i/o status
    INTEGER       I,J,L        ! counter and indices

    LOGICAL    :: CFLAG = .FALSE.      ! true: open area-to-point file
    LOGICAL    :: DFLAG = .FALSE.      ! true: open day-specific file
    LOGICAL    :: HFLAG = .FALSE.      ! true: open hour-specific file
    LOGICAL    :: IFLAG = .FALSE.      ! true: open annual/average inventory
    LOGICAL    :: NFLAG = .FALSE.      ! true: open non-HAP inclusions/exclusions
    LOGICAL    :: MFLAG = .FALSE.      ! true: treat all sources as treated

    CHARACTER(NAMLEN3) ANAME
    CHARACTER(NAMLEN3) NAMBUF          ! file name buffer
    CHARACTER(NAMLEN3) INAME           ! tmp name for inven file of unknown fmt
    CHARACTER(256)     MESG            ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'OPENINVIN'     ! program name

    !***********************************************************************
    !   begin body of subroutine OPENINVIN

    !.......  Get value of these controls from the environment

    !.......  Get setup up for later processing of a formula for adding a variable
    MESG = 'Define a formula to compute new pollutant(s)'
    CALL ENVSTR( 'SMKINVEN_FORMULA', MESG, ' ', VAR_FORMULA, IOS )
    IF ( IOS .NE. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMKINVEN_FORMULA"', 2 )
    END IF
    L = LEN_TRIM( VAR_FORMULA )
    IF( L .GT. 0 ) THEN
    !.......  Figure out how many variables there are based on the
    !               number of commas found in the string.
        NCOMP = 1
        DO I = 1, L
            IF( VAR_FORMULA( I:I ) == ',' ) NCOMP = NCOMP + 1
        ENDDO
    END IF

    MESG = 'Import average inventory data'
    IFLAG = ENVYN ( 'IMPORT_AVEINV_YN', MESG, .TRUE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "IMPORT_AVEINV_YN"', 2 )
    END IF

    MESG = 'Import day-specific data'
    DFLAG = ENVYN ( 'DAY_SPECIFIC_YN', MESG, .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "DAY_SPECIFIC_YN"', 2 )
    END IF
    DAYINVFLAG = DFLAG

    MESG = 'Import hour-specific data'
    HFLAG = ENVYN ( 'HOUR_SPECIFIC_YN', MESG, .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "HOUR_SPECIFIC_YN"', 2 )
    END IF
    HRLINVFLAG = HFLAG

    IF( CATEGORY .EQ. 'MOBILE' ) THEN
        MESG = 'Use referenced SCC activity inventory file'
        SCCMAPFLAG = ENVYN ( 'USE_REF_SCC_YN', MESG, .FALSE., IOS )

        IF( SCCMAPFLAG ) THEN
            MESG = 'Enter logical name for reference SCC input file'
            MDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'SCCXREF', PROGNAME )
            CALL RDSCCMAP( MDEV )

            MESG = 'Exclude SCCs not found in SCCXREF input file'
            EXCLSCCFLAG = ENVYN ( 'EXCLUDE_REF_SCC_YN', MESG, .FALSE., IOS )

        END IF
    END IF

    IF ( CATEGORY .EQ. 'POINT' ) THEN
        MESG = 'Import gridded MEDS-formatted inventory file'
        MEDSFLAG = ENVYN ( 'IMPORT_MEDS_YN', MESG, .FALSE., IOS )

        IF( MEDSFLAG ) THEN
            MESG = 'WARNING: MUST process daily or hourly MEDS inventory'
            IF( .NOT. DAYINVFLAG .AND. .NOT. HRLINVFLAG ) THEN
                CALL M3MESG( MESG )
            END IF
        END IF
    END IF

    IF( DAYINVFLAG .OR. HRLINVFLAG ) THEN
        MESG = 'Use a daily FF10 inventory data as an annual FF10 inventory'
        FF10INVFLAG = ENVYN ( 'FF10_AVEDAY_ANNINV_YN', MESG, .FALSE., IOS )
    END IF

    MESG = 'Define the processing method of combining haradous ' //    &
           'air pollutants with criteria VOC.'
    CALL ENVSTR( 'SMK_PROCESS_HAPS', MESG, ' ', PROC_HAPS, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK_PROCESS_HAPS"', 2 )
    END IF

    SELECT CASE( PROC_HAPS )
      CASE( 'ALL' )
        MESG = 'Treat all sources as integrate sources to compute NONHAP[VOC|TOG].'
        NFLAG = .FALSE.

      CASE( 'NONE' )
        MESG = 'Treat all sources as non-integrate sources.'
        NFLAG = .FALSE.

      CASE( 'PARTIAL' )
        MESG = 'Partially treat sources as either integrate or non-integrate sources.'
        NFLAG = .TRUE.

      CASE DEFAULT
        MESG = 'No processing of combining criteria VOC with hazardous air pollutants (HAP).'
        NFLAG = .FALSE.

    END SELECT

    CALL M3MSG2( MESG )

    IF ( CATEGORY .EQ. 'AREA' ) THEN
        MESG = 'Read and use area-to-point factors file'
        CFLAG = ENVYN ( 'SMK_ARTOPNT_YN', MESG, .FALSE., IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK_ARTOPNT_YN"', 2 )
        END IF

        MESG = 'Import gridded I/O API inventory data'
        APIFLAG = ENVYN ( 'IMPORT_GRDIOAPI_YN', MESG, .FALSE., IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "IMPORT_GRDIOAPI_YN"', 2 )
        END IF

        MESG = 'Import gridded native I/O API inventory data'
        NCDFLAG = ENVYN ( 'IMPORT_GRDNETCDF_YN', MESG, .FALSE., IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "IMPORT_GRDNETCDF_YN"', 2 )
        END IF
    END IF

    !.......  Define the unit and temporal resolution of raw NetCDF inventory files
    IF( NCDFLAG ) THEN
        NETCDFUNIT = ''
        MESG = 'Define the unit of NetCDF gridded inventory pollutant'
        CALL ENVSTR( 'NETCDF_POL_UNIT', MESG, ' ', NETCDFUNIT, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "NETCDF_POL_UNIT"', 2 )
        ELSE IF( NETCDFUNIT /= 'kg m-2 s-1' ) THEN
            MESG = 'ERROR: Unit of pollutant MUST be "kg m-2 s-1"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
    END IF

    !.......  Make sure gridded point source file is not attempted
    IF ( ( CATEGORY .EQ. 'POINT' .OR. CATEGORY .EQ. 'MOBILE'    ) .AND.    &
           ( APIFLAG    .OR.  NCDFLAG  )          ) THEN
        MESG = 'Cannot import gridded mobile or point source data.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  When gridded data are imported, override other settings
    IF( APIFLAG .OR. NCDFLAG ) THEN
        CFLAG = .FALSE.
        IFLAG = .FALSE.
        DFLAG = .FALSE.
        HFLAG = .FALSE.
    END IF

    !.......  Abort if no settings set to read data
    IF( .NOT. IFLAG .AND.    &
        .NOT. DFLAG .AND.    &
        .NOT. HFLAG .AND.    &
        .NOT. APIFLAG .AND.  &
        .NOT. NCDFLAG        ) THEN

        MESG = 'ERROR: Environment settings indicate no files are to be read'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    !.......  Find name name of raw inventory file
    J = INDEX1( CATEGORY, NCAT, CATLIST )
    IF( J .LE. 0 ) THEN
        MESG = 'INTERNAL ERROR: Do not know about category ' //    &
               TRIM( CATEGORY ) // ' in program ' // PROGNAME
        CALL M3MSG2( MESG )
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    ELSE

        INNAME = ANAMLIST( J )
        ENAME  = GNAMLIST( J )
        IDNAME = DNAMLIST( J )
        IHNAME = HNAMLIST( J )

    END IF

    !.......   Get IOAPI gridded inventory file (without SCCs)
    IF( APIFLAG ) THEN

        MESG = 'Enter logical name of the GRIDDED ' //    &
               TRIM( CATEGORY ) // ' INVENTORY ' // 'file'

        ENAME = PROMPTMFILE( MESG, FSREAD3, ENAME, PROGNAME )

    !.......   Get IOAPI gridded masking file (Country code and associated time zones)
    ELSE IF( NCDFLAG ) THEN

        MESG = 'Enter logical name of the GRID MASK input file '

        ENAME = PROMPTMFILE( MESG, FSREAD3, 'GRIDMASK', PROGNAME )

    !.......  Get ASCII file name and open average input inventory file when
    !           inventory is to be imported
    ELSE IF( IFLAG ) THEN
        MESG = 'Enter logical name of the RAW ' //    &
               TRIM( CATEGORY ) // ' AVERAGE INVENTORY ' // 'file'

        ADEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INNAME, PROGNAME )

    !.......  If SMOKE inventory already exists, open files for reading later
    ELSE

    !.......  Get input inventory file names given source category
    !.......  Use NAMBUF for HP safety
        CALL GETINAME( CATEGORY, ENAME, ANAME )

    !.......  Prompt for and open inventory file
        INAME = ENAME
        MESG = 'Enter logical name for the MAP INVENTORY file'
        IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INAME, PROGNAME )

    !.......  Open and read map file
        CALL RDINVMAP( INAME, IDEV, ENAME, ANAME, SDEV )

    END IF

    !.......  Get file name and open daily input inventory file
    IF( DFLAG ) THEN
        MESG = 'Enter logical name of the RAW ' //    &
               TRIM( CATEGORY ) // ' DAILY INVENTORY ' // 'file'

        DDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., IDNAME, PROGNAME )
    END IF

    !.......  Get file name and open daily input inventory file
    IF( HFLAG ) THEN
        MESG = 'Enter logical name of the RAW ' //    &
               TRIM( CATEGORY ) // ' HOURLY INVENTORY ' // 'file'

        HDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., IHNAME, PROGNAME )
    END IF

    IF( IFLAG ) THEN

    !.......  Open category-specific inputs
        SELECT CASE( CATEGORY )
          CASE( 'POINT' )

    !.......  Get file name for input replacement stack parameters file
            MESG = 'Enter logical name for REPLACEMENT STACK PARAMETERS file'
            RDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'PSTK', PROGNAME )

        END SELECT

    END IF

    !.......  Get file name for country, state, and county file, with time
    !           zones
    IF( .NOT. APIFLAG .AND. .NOT. USEEXPGEO() ) THEN
        ZDEV = PROMPTFFILE(    &
               'Enter logical name for COUNTRY, STATE, AND COUNTY file',    &
               .TRUE., .TRUE., 'COSTCY', PROGNAME )
    END IF

    !.......  Get SCC descriptions if needed for reporting purposes
    IF ( HFLAG .OR. CFLAG ) THEN
        CDEV = PROMPTFFILE(    &
               'Enter logical name for SCC DESCRIPTION file ',    &
               .TRUE., .TRUE., 'SCCDESC', PROGNAME )
    END IF

    !.......  Get ORIS descriptions file
    IF( HFLAG ) THEN
        ODEV = PROMPTFFILE(    &
               'Enter logical name for ORIS DESCRIPTION file ',    &
               .TRUE., .TRUE., 'ORISDESC', PROGNAME )
    END IF

    !.......  Get file name for non-HAP exclusion file
    IF( NFLAG ) THEN
        MESG = 'Enter logical name for NHAPEXCLUDE file'
        UDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'NHAPEXCLUDE', PROGNAME )
    END IF

    !.......  Get file name for area-to-point factors file
    IF( CFLAG ) THEN
        MESG = 'Enter logical name for AREA-TO-POINT FACTORS file'
        YDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'ARTOPNT', PROGNAME )
    END IF

    !.......  Get file name for inventory pollutants codes/names
    MESG = 'Enter logical name for INVENTORY DATA TABLE file'
    PDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'INVTABLE', PROGNAME )

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE OPENINVIN


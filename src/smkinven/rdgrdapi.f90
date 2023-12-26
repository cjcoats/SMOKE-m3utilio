
SUBROUTINE RDGRDAPI( FNAME, GRDNM )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine reads an I/O API NetCDF gridded inventory file.
    !
    !  PRECONDITIONS REQUIRED:
    !      Input file logical name FNAME opened
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 10/2000 by M. Houyoux
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

    !.......   MODULES for public variables
    !.......   This module is the inventory arrays
    USE MODSOURC, ONLY: POLVAL, TZONES, CIFIP, CELLID, TPFLAG, INVYR,    &
                        NPCNT, CSCC, IPOSCOD, CSOURC

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: NIPOL, NIPPA, NPPOL, CATEGORY, NEM, NDY,    &
                       EIIDX, EINAM, EANAM, EAUNIT, EADESC, NSRC

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: FNAME         ! logical name of input file
    CHARACTER(*), INTENT(OUT) :: GRDNM         ! grid name for input data

    !.......   EXTERNAL FUNCTIONS and their descriptions:

    INTEGER, EXTERNAL :: GETIFDSC

    !.......   Local allocatable arrays
    REAL, ALLOCATABLE :: INDATA( : )       ! tmp gridded input data

    !.......   Other local variables
    INTEGER         C, ES, J, K, L, S, V         !  counters and indices

    INTEGER         IDX             !  emissions array index (ann or ave day)
    INTEGER      :: INY = 0         !  tmp inventory year
    INTEGER         IOS             !  i/o status
    INTEGER      :: NCELL           !  tmp cell numbers
    INTEGER         NEDIM1          !  1st dimension for sparse emis arrays
    INTEGER         NVARS           !  no. variables in gridded input file
    INTEGER      :: TPF             !  tmp temporal adjustments setting
    INTEGER      :: TZONE = -50     !  tmp time zone
    INTEGER         WKSET           !  setting for wkly profile TPFLAG component

    LOGICAL      :: EFLAG  = .FALSE.     ! true: error occured
    LOGICAL      :: DFLAG  = .FALSE.     ! true: weekday (not full week) nrmlizr
    LOGICAL      :: TVFLAG = .FALSE.     ! true: time zone is a variable

    CHARACTER(300)  MESG            !  message buffer

    CHARACTER(CELLEN3) CCELL        ! tmp cell ID
    CHARACTER(POLLEN3) CCOD         ! character pollutant index
    CHARACTER(FIPLEN3) CFIP         ! character FIP code
    CHARACTER(NAMLEN3) VBUF         ! tmp variable name
    CHARACTER(SCCLEN3) SCCZERO      ! default source category code

    CHARACTER(16), PARAMETER :: PROGNAME =  'RDGRDAPI'     ! program name

    !***********************************************************************
    !   begin body of subroutine RDGRDAPI

    !.......  Get setting for interpreting weekly temporal profiles from the
    !           environment.
    DFLAG = .FALSE.
    MESG = 'Use weekdays only to normalize weekly profiles'
    DFLAG = ENVYN( 'WKDAY_NORMALIZE', MESG, DFLAG, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "WKDAY_NORMALIZE"', 2 )
    END IF

    !.......  Set weekly profile interpretation flag...
    !.......  Weekday normalized
    IF( DFLAG ) THEN
        WKSET = WDTPFAC
        MESG = 'NOTE: Setting inventory to use weekday '//    &
               'normalizer for weekly profiles'

    !.......  Full-week normalized
    ELSE
        WKSET = WTPRFAC
        MESG = 'NOTE: Setting inventory to use full-week '//    &
               'normalizer for weekly profiles'

    END IF

    !.......  Write message
    CALL M3MSG2( MESG )

    !.......  Set default inventory characteristics (declared in MODINFO)
    CALL INITINFO( IOGFMT )

    !.......  Read the header of the gridded I/O API file
    IF( .NOT. DESC3( FNAME ) ) THEN

        MESG = 'Could not get description of file "' // TRIM( FNAME ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    !.......  Give error if there are no variables
    IF( NVARS3D .LT. 1 ) THEN

        MESG = 'No variables found in gridded input file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    !.......  Give warning if time zone is not set in the header and it is not
    !           in the variable list
    TZONE = GETIFDSC( FDESC3D, '/TZONE/', .FALSE. )
    K = INDEX1( 'TZONES', NVARS3D, VNAME3D )
    IF( TZONE .LT. 0 .AND.    &
        K     .LE. 0       ) THEN

        WRITE( MESG,94010 ) 'WARNING: Time zone not set in ' //    &
               'header or variable list.' // CRLF() // BLANK10 //    &
               'Assuming time zone 0 (GMT) for data.'
        CALL M3MSG2( MESG )
        TZONE = 0

    ELSE IF( K .GT. 0 ) THEN
        TVFLAG = .TRUE.

    END IF

    !.......  Give warning if layers are greater than 1
    IF( NLAYS3D .GT. 1 ) THEN

        WRITE( MESG,94010 ) 'WARNING: Only the first layer out of', NLAYS3D,    &
                            'will be imported.'
        CALL M3MSG2( MESG )

    END IF

    !.......  Give warning if time steps are greater than 1
    IF( MXREC3D .GT. 1 ) THEN

        WRITE( MESG,94010 ) 'WARNING: Only the first time step out of', MXREC3D, &
                            'will be imported.'
        CALL M3MSG2( MESG )

    END IF

    !.......  Give error if year of data is not provided
    IF( SDATE3D .LT. 1900000 ) THEN

        WRITE( MESG,94010 ) 'Cannot determine year from invalid start date', SDATE3D, '.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    !.......  Otherwise set the year of the input data
    ELSE
        INY = SDATE3D / 1000

    END IF

    !.......  Set temporal flag based on time step in file
    IF( TSTEP3D .EQ. 87600000 ) THEN         ! Annual data

        TPF = MTPRFAC * WKSET
        IDX = NEM

    ELSE IF( TSTEP3D .EQ. 240000 ) THEN      ! Average day data

        TPF = WKSET
        IDX = NDY

    !.......  Give error if time step is not one of the recognized values.
    ELSE
        MESG = 'Time step needs to be 87600000 or 240000 HHMMSS ' //    &
               'to determine temporal approach'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    !.......  Set the number of cells, variables, grid name, etc.
    NCELL  = NROWS3D * NCOLS3D
    NVARS  = NVARS3D
    NIPOL  = NVARS3D
    GRDNM  = GDNAM3D

    !.......  Reset variable number if there are any special variables
    IF( TVFLAG ) THEN
        NIPOL = NIPOL - 1
    END IF

    !.......  Set dependent values (no. sources, pol/act, and dimension)
    NSRC   = NCELL
    NIPPA  = NIPOL
    NEDIM1 = NSRC * NIPPA

    !.......  Allocate memory for variables
    !.......  NOTE - Both EINAM and EANAM are created to support WRINVEMIS
    ALLOCATE( EIIDX( NIPOL ),    &
              EINAM( NIPOL ),    &
              EANAM( NIPPA ),    &
             EAUNIT( NIPPA ),    &
             EADESC( NIPPA ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EIIDX...EADESC', PROGNAME )

    J = 0
    DO V = 1, NVARS

    !.......  Skip special variable names and set output variable info
        IF( VNAME3D( V ) .NE. 'TZONES' ) THEN

            J = J + 1
            EIIDX ( J ) = V
            EINAM ( J ) = VNAME3D( V )
            EANAM ( J ) = VNAME3D( V )
            EAUNIT( J ) = UNITS3D( V )
            EADESC( J ) = VDESC3D( V )

        END IF

    END DO

    !.......  Allocate memory for (sorted) output inventory characteristics.
    !.......  The sorted arrays can be allocated right away because the only
    !           source characteristic for this type of data is the grid cell, and
    !           it is sorted already.
    CALL SRCMEM( CATEGORY, 'SORTED', .TRUE., .FALSE., NSRC, NEDIM1, NPPOL )

    CALL SRCMEM( CATEGORY, 'SORTED', .TRUE., .TRUE.,  NSRC, NEDIM1, NPPOL )

    !.......  Initialize emissions data
    POLVAL = BADVAL3        ! array

    !.......  Initialize tmp source characteristics
    CFIP    = REPEAT( '0', FIPLEN3 )
    SCCZERO = REPEAT( '0', SCCLEN3 )
    CCOD    = ' '

    !.......  Loop through cells (same as sources) and store data as sources
    ES  = 0
    DO S = 1, NSRC

        CIFIP ( S ) = CFIP
        CELLID( S ) = S
        TPFLAG( S ) = TPF
        INVYR ( S ) = INY
        NPCNT ( S ) = NIPPA
        CSCC  ( S ) = SCCZERO

        !.......  If time zone is not a variable in the input file already, then
        !         then set it based on the header
        IF( .NOT. TVFLAG ) THEN
            TZONES( S ) = TZONE
        END IF

        WRITE( CCELL, 94130 ) S

        CALL BLDCSRC( CFIP, SCCZERO, CCELL, CHRBLNK3,  &
                      CHRBLNK3, CHRBLNK3, CHRBLNK3,    &
                      CCOD, CSOURC( S ) )

    END DO

    !.......  If time zone is a variable, store it
    IF( TVFLAG ) THEN

        IF( .NOT. READ3( FNAME, 'TZONES', 1, SDATE3D,STIME3D, TZONES ) ) THEN

            !.......  Write error if data could not be read
            EFLAG = .TRUE.
            MESG = 'ERROR: Could not read "TZONES" from file ' // TRIM( FNAME ) // '".'
            CALL M3MSG2( MESG )

        END IF

    END IF

    !.......  Allocate local memory for temporary gridded data file
    ALLOCATE( INDATA( NCELL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'INDATA', PROGNAME )

    !.......  Read data from gridded file and store in appropriate data structure
    !         for use by the rest of the program
    DO V = 1, NVARS

        VBUF = EANAM( V )

        !.......  Skip special variable "TZONES"
        IF( VBUF .EQ. 'TZONES' ) CYCLE

        IF( .NOT. READ3( FNAME, VBUF, 1, SDATE3D,STIME3D, INDATA ) ) THEN

            !.......  Write error if data could not be read
            EFLAG = .TRUE.
            MESG = 'ERROR: Could not read "' // TRIM( VBUF )//    &
                   '" from file "' // TRIM( FNAME ) // '".'
            CALL M3MSG2( MESG )

        ELSE

            DO C = 1, NCELL

                ES = ( C-1 ) * NIPPA + V
                IPOSCOD( ES )     = V
                POLVAL ( ES,IDX ) = INDATA( C )

            END DO

        END IF

    END DO

    !.......  Abort if there was a reading error
    IF( EFLAG ) THEN
        MESG = 'Problem reading gridded inventory file "' //TRIM( FNAME )  // '".'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Deallocate local memory
    DEALLOCATE( INDATA )

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

94125 FORMAT( I5 )

94130 FORMAT( I8 )

END SUBROUTINE RDGRDAPI

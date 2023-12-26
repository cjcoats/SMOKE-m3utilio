
PROGRAM LAYALLOC

!***********************************************************************
!  program body starts at line
!
!  DESCRIPTION:
!
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!       I/O API
!
!  REVISION  HISTORY:
!       Version ??/???? by ???
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!***************************************************************************
!
! COPYRIGHT (C) 2008, Institute for the Environment, UNC-CH
! All Rights Reserved
!
!
!*************************************************************************
    USE M3UTILIO

    IMPLICIT NONE

!...........   INCLUDES:


    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   EXTERNAL FUNCTIONS and their descriptions:

    INTEGER, EXTERNAL :: GETFLINE
!
!...........   PARAMETERS:

    CHARACTER*16, PARAMETER :: PROGNAME = 'LAYALLOC'   !  program name
    CHARACTER*50, PARAMETER :: CVSW = '$Name SMOKEv5.0_Jun2023$' ! CVS release tag
    REAL,         PARAMETER :: CONVPA = 1.0E-2  ! conversion factor for Pa to mb

!...........   LOCAL VARIABLES and their descriptions:

    CHARACTER*80    GDESC    !  grid description

!...........   Logical file names and unit numbers

    INTEGER         LDEV    !  log-device
    INTEGER         FDEV       ! input file for user-defined vertical distribution

    CHARACTER*16    SNAME      ! grid-point src met file
    CHARACTER*16    XNAME      ! cross-point layered met file
    CHARACTER*16    ONAME      ! Output factor file
    CHARACTER*16    CNAME      ! Input emission file

!...........   Aarrays
    REAL              UBOT( MXLAYS3+1 )     ! user-defined layer bottom
    REAL              UTOP( MXLAYS3+1 )     ! user-defined layer top
    REAL             UFRAC( MXLAYS3+1 )     ! user-defined layer fraction
    REAL             LFRAC( MXLAYS3+1 )     ! model layer fractions
    REAL             VGLVLS( 0:MXLAYS3 )    ! gridded mask values to be output
    REAL     , ALLOCATABLE ::    ZZF( :,:,: )   ! layer's full height (m)
    REAL     , ALLOCATABLE :: TMPBUF( :,:,: )  ! Tmp buffer
    REAL     , ALLOCATABLE :: TMP3D ( :,:,: )

!...........   Other local variables

    INTEGER         I, J, K, M, N, C, F, S, R, NG, L, NL ! counters, subscripts
    INTEGER         IOS                     !  I/O status
    INTEGER         KREC,IREC                    !  input line (record) number
    INTEGER         JHR, JMAX
    INTEGER         IHR
    INTEGER         NVARS
    INTEGER         NFILES
    INTEGER         NSTEPS
    INTEGER         NCOLS
    INTEGER         NROWS
    INTEGER         NLAYS
    INTEGER         ULAYS           ! no of user-defined layers
    INTEGER         VGTYP
    INTEGER         LTOP
    INTEGER         LBOT, DDP

    INTEGER      :: JDATE = 0 ! Julian start date (YYYYDDD)
    INTEGER      :: MDATE = 0 ! Julian start date (YYYYDDD) for met file
    INTEGER      :: JTIME = 0 ! start time (HHMMSS)
    INTEGER      :: MTIME = 0 ! start time (HHMMSS) for met file
    INTEGER         TSTEP     ! output time step

    REAL            VGTOP
    REAL            ZBOT, ZTOP, ZFRAC, TFRAC
    REAL            PFRAC, PDIFF
    REAL            LTOT

    LOGICAL         EFLAG      ! Input error flat
    LOGICAL    ::   FIRSTIME = .TRUE.  ! first time flag

    CHARACTER*5     FFORMAT !  temporary indicator for input formats
    CHARACTER*16    CVAR, TFILE, TLAY, VNM
    CHARACTER*100   LINE
    CHARACTER*256   MESG
    CHARACTER*30    SEGMENT( 5 )

!**********************************************************************
!   begin body of program

    LDEV = INIT3()

!.........  Write out copywrite, version, web address, header info, and prompt
!           to continue running the program.

    CALL INITEM( LDEV, CVSW, PROGNAME )

!.........  Get input files

    MESG = 'Enter name for CROSS-POINT LAYERED MET file'
    XNAME = PROMPTMFILE( MESG, FSREAD3, 'MET_CRO_3D', PROGNAME )

!.........  Read description of 3d file for defining layer structure

    IF ( .NOT. DESC3( XNAME ) ) THEN
        MESG = 'Could not get description of file "' //&
        &       TRIM( XNAME ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ENDIF
    MDATE  = SDATE3D
    MTIME  = STIME3D
    NSTEPS = MXREC3D
    NLAYS  = NLAYS3D
    VGTYP  = VGTYP3D
    VGTOP  = VGTOP3D
    NCOLS  = NCOLS3D
    NROWS  = NROWS3D
    VGLVLS = VGLVS3D   ! array

!.........  Get input layer fractions file

    MESG = 'Enter logical name for gridded input file'
    FDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'LAYER_FRACTION',&
    &                    PROGNAME )
    ULAYS = GETFLINE( FDEV, 'User-defined Vertical Layer fraction file' )
    IF ( ULAYS .GT. MXLAYS3+1 ) THEN
        WRITE( MESG, '(A,I10)' )&
        &   'Max layers exceeded in LAYER_FRACTION file:', ULAYS
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Store local layer information
    ZZF    = 0.0
    LFRAC  = 0.0
    UTOP   = 0.0
    UBOT   = 0.0
    UFRAC  = 0.0

!.........  Read the layers file and store in array

    NL   = 0
    IREC = 0

    DO L = 1, ULAYS

        READ ( FDEV, 93000, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF ( IOS .GT. 0 ) THEN
            WRITE( MESG, 94020)&
            &   'I/O error', IOS, 'reading LAYERS file at line',IREC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!.............  Left adjust line
        LINE = TRIM( LINE )

!.............  Skip blank and comment lines
        IF( LINE(1:1) == '#' ) CYCLE

!.............  Get line
        CALL PARSLINE( LINE, 5, SEGMENT )
        NL = NL + 1
        UBOT ( NL ) = STR2REAL( SEGMENT( 2 ) )
        UTOP ( NL ) = STR2REAL( SEGMENT( 3 ) )
        UFRAC( NL ) = STR2REAL( SEGMENT( 4 ) )

    END DO

    ULAYS = NL

!.........  Read 2-D emissions file
    CNAME = PROMPTMFILE(&
    &        'Enter name for netCDF 2-d emissions file',&
    &        FSREAD3, 'INFILE', PROGNAME )

    IF ( .NOT. DESC3( CNAME ) ) THEN
        MESG = 'Could not get description of file "' //&
        &TRIM( CNAME ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ELSE IF ( .NOT.FILCHK3( CNAME, GRDDED3,&
    &                        NCOLS, NROWS, 1, NTHIK3D ) ) THEN
        MESG = 'Bad dimensions for file "' //&
        &TRIM( CNAME ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Allocate arrays
    ALLOCATE(   ZZF( NCOLS,NROWS,NLAYS ),&
    &          TMP3D( NCOLS,NROWS,NLAYS ),&
    &          TMPBUF( NCOLS,NROWS,1 ), STAT= IOS)
    CALL CHECKMEM( IOS, 'TMPBUF', PROGNAME )

    TMP3D  = 0.0
    TMPBUF = 0.0

!.........  Multiply layer fracitons by 2-d emissions
    JDATE = SDATE3D
    JTIME = STIME3D
    DO IHR = 1, NSTEPS

!.............  Read layer's top height (meter)
        IF( .NOT. READ3( XNAME, 'ZF', -1,&
        &                 JDATE, JTIME, ZZF ) ) THEN
            MESG = 'ERROR : Could not read ZF from file '//XNAME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

!.............  Loop through variables
        DO K = 1,NVARS3D

            TMPBUF = 0.0    ! array
            VNM = VNAME3D( K )

            IF( .NOT. READ3( CNAME, VNM, 1, JDATE, JTIME,&
            &    TMPBUF ) ) THEN
                MESG = 'ERROR : Could not read ' //VNAME3D( K ) //&
                &       ' from file ' // CNAME
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            DO R = 1, NROWS3D

                DO C = 1, NCOLS3D

!.........................  Compute layer fractions given top an bottom
!                           height of the user-defined layers
                    LFRAC = 0.0
                    TFRAC = 0.0

                    DO I = 1, ULAYS

!.............................  User-define layers
                        ZTOP  = UTOP( I )
                        ZBOT  = UBOT( I )
                        ZFRAC = UFRAC( I )

!.............................  Sum of user-defined layer fractions
                        TFRAC = TFRAC + ZFRAC

                        DO L = 1, NLAYS - 1
                            IF ( ZBOT <= ZZF( C,R,L ) ) THEN
                                LBOT = L
                                GO TO  111   ! end loop and skip reset of LBOT
                            END IF
                        END DO
                        LBOT = NLAYS           !  fallback

111                     CONTINUE                !  loop exit:  bottom found at LBOT

                        IF ( ZTOP <= ZZF( C,R,LBOT ) ) THEN  !  plume in this layer

                            PFRAC = 1.0 * ZFRAC
                            LFRAC( LBOT ) = LFRAC( LBOT ) + PFRAC
                            LTOP = LBOT

                        ELSE IF( LBOT == NLAYS ) THEN    ! plume above top layer

                            PFRAC = 1.0 * ZFRAC
                            LFRAC( LBOT ) = LFRAC( LBOT ) + PFRAC
                            LTOP = NLAYS

                        ELSE                               ! plume crosses layers

                            DO L = LBOT + 1, NLAYS
                                IF ( ZTOP <= ZZF( C,R,L ) ) THEN
                                    LTOP = L
                                    GO TO 222  ! end loop and skip reset of LTOP
                                END IF
                            END DO
                            LTOP = NLAYS

222                         CONTINUE

!.................................  Calculate between layer
                            PDIFF = ZTOP - ZBOT

!.................................  Calculate a fraction for the bottom layer
                            PFRAC = ( ( ZZF( C,R,LBOT ) - ZBOT )&
                            &        / PDIFF ) * ZFRAC
                            LFRAC( LBOT ) = LFRAC( LBOT ) + PFRAC

!.................................  Calculate a fraction for the top layer
                            PFRAC = ( ( ZTOP - ZZF( C,R,LTOP-1 ) )&
                            &        / PDIFF ) * ZFRAC
                            LFRAC( LTOP ) = LFRAC( LTOP ) + PFRAC

                            DDP = LTOP - LBOT

                            IF( DDP >= 2 ) THEN

                                DO L = LBOT+1, LTOP-1 !  layers in plume

                                    PFRAC =&
                                    &   ( ( ZZF(C,R,L) -ZZF(C,R,L-1) )&
                                    &    / PDIFF ) * ZFRAC
                                    LFRAC( L ) = LFRAC( L ) + PFRAC

                                END DO

                            ENDIF

                        END IF          !  if ztop in same layer as zbot, or not

                    ENDDO     ! loop over user-define layers

                    IF( TFRAC > 1.001 ) THEN
                        MESG = 'ERROR: The total user-defined layer fractions ' //&
                        &       ' can not be greater than 1.0'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    ENDIF

!......................... Before applying layer fractions make sure that they add to 1.0
                    LTOT = 0.0

                    DO L = 1, NLAYS
                        LTOT = LTOT + LFRAC( L )
                    ENDDO

                    IF ( LTOT < 0.999 ) THEN
                        MESG = 'ERROR: Calculated layer fractions '//&
                        &    'are less than 1.0. Emissions will be '//&
                        &    'dropped.'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    ELSE IF ( LTOT > 1.001 ) THEN
                        MESG = 'ERROR: Calculated layer fractions '//&
                        &       'are greater than 1.0. Emissions '//&
                        &       'will be increased.'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    ENDIF

!.........................  Apply layered fraction values to emissions
                    DO L = 1, NLAYS
                        TMP3D( C,R,L ) = TMPBUF( C,R,1 ) * LFRAC( L )
                    ENDDO

                ENDDO

            ENDDO

!.............  Get output file name using environment variable
            IF( FIRSTIME ) THEN

!.................  Define top layer for output file
                NLAYS3D = LTOP
                VGTYP3D = VGTYP
                VGTOP3D = VGTOP
                VGLVS3D = VGLVLS
                MESG = 'Enter logical name for output file'
                ONAME = PROMPTMFILE( MESG, FSUNKN3, 'OUTFILE',&
                &                     PROGNAME )
                FIRSTIME = .FALSE.

            ENDIF

            IF ( .NOT. WRITE3( ONAME, VNM, JDATE, JTIME, TMP3D ) )&
            &                                                    THEN
                WRITE( MESG, 93000 ) 'Could not write to "'&
                &       // TRIM( ONAME ) // '".'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ENDIF

!.................  Re-initialize tmp array
            TMP3D = 0.0

        ENDDO

        CALL NEXTIME( JDATE, JTIME, 10000 )

    ENDDO

!.........   End of program:
    CALL M3EXIT( PROGNAME, 0, 0, 'Completion of '//PROGNAME, 0 )


!******************  FORMAT  STATEMENTS   ******************************
!...........   Internal buffering formats............ 94xxx

94020 FORMAT( 10( A, I3.0 ) )
93000 FORMAT( A )

END PROGRAM LAYALLOC

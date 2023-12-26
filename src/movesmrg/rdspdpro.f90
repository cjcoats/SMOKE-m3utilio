
SUBROUTINE RDSPDPRO( SPDEV )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !       Reads hourly speed data from SPDPRO file
    !
    !  PRECONDITIONS REQUIRED:
    !       SPDEV must be opened
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:  none
    !
    !  REVISION  HISTORY:
    !       09/10: Created by C. Seppanen
    !       08/15: Modified by B.H. Baek
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
    !***********************************************************************
    USE M3UTILIO

    !.......  MODULES for public variables
    !.......  This module contains data structures and flags specific to Movesmrg
    USE MODMVSMRG, ONLY: SPDPRO

    !.......  This module contains the lists of unique source characteristics
    USE MODLISTS, ONLY: NINVIFIP, INVCFIP, NINVSCC, INVSCC

    !.......  This module is for mobile-specific data
    USE MODMOBIL, ONLY: NSCCMAP, SCCMAPFLAG, SCCMAPLIST, EXCLSCCFLAG

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    INTEGER, INTENT (IN) :: SPDEV                 ! SPDPRO file unit no.

    !.......   EXTERNAL FUNCTIONS and their descriptions:
    LOGICAL, EXTERNAL :: BLKORCMT
    LOGICAL, EXTERNAL :: CHKINT
    LOGICAL, EXTERNAL :: CHKREAL
    INTEGER, EXTERNAL :: GETFLINE

    !.......   Local arrays
    CHARACTER(20)  SEGMENT( 53 )              ! parsed input line

    !.......   Other local variables
    INTEGER         I, J, JJ, KK         ! counters and indexes
    INTEGER         IOS             ! error status
    INTEGER      :: IREC = 0        ! record counter
    INTEGER         FIPIDX          ! current FIPS index
    INTEGER         NSCC            ! no of reference SCCs
    INTEGER         SCCIDX          ! current SCC index
    INTEGER         NLINES          ! number of lines
    INTEGER         MDEV            ! open SCCXREF input file
    INTEGER     NWARN1, NWARN2, MXWARN       ! number of warning messages

    REAL            SPDVAL          ! hourly speed value

    LOGICAL         EFLAG           ! true: error found

    CHARACTER(1060)    LINE         ! line buffer
    CHARACTER(300)     MESG         ! message buffer
    CHARACTER(FIPLEN3) CFIP         ! current FIPS
    CHARACTER(SCCLEN3) SCC          ! current SCC
    CHARACTER(10)      KEYWORD      ! temperature keyword

    CHARACTER(16), PARAMETER :: PROGNAME = 'RDSPDPRO'       ! program name

    !***********************************************************************
    !   begin body of subroutine RDSPDPRO

    !.......  Get maximum number of warnings
    EFLAG  = .FALSE.
    MXWARN = ENVINT( WARNSET, ' ', 100, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble WARNSET', 2 )
    END IF

    !.......  Read cross-refrenced SCC input file
    MESG = 'Use referenced SCC activity inventory file'
    SCCMAPFLAG = ENVYN ( 'USE_REF_SCC_YN', MESG, .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "USE_REF_SCC_YN"', 2 )
    END IF

    IF( SCCMAPFLAG ) THEN
        MESG = 'Enter logical name for reference SCC input file'
        MDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'SCCXREF', PROGNAME )
        CALL RDSCCMAP( MDEV )

        MESG = 'Exclude SCCs not found in SCCXREF input file'
        EXCLSCCFLAG = ENVYN ( 'EXCLUDE_REF_SCC_YN', MESG, .FALSE., IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "EXCLUDE_REF_SCC_YN"', 2 )
        END IF

    END IF

    !.......  Allocate storage based on number of FIPs and SCCs in inventory
    ALLOCATE( SPDPRO( NINVIFIP, NINVSCC, 2, 24 ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SPDPRO', PROGNAME )
    SPDPRO = BADVAL3       ! array

    !.......  Get the number of lines in the file
    NLINES = GETFLINE( SPDEV, 'SPDPRO file' )

    !.......  Read through file and store hourly data
    NWARN1 = 0
    NWARN2 = 0
    DO I = 1, NLINES

        READ( SPDEV, 93000, END=999, IOSTAT=IOS ) LINE

        IREC = IREC + 1

        IF( IOS .NE. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010 ) 'I/O error', IOS,&
              'reading hourly speed file at line', IREC
            CALL M3MESG( MESG )
            CYCLE
        END IF

        !.......  Skip blank or comment lines
        IF( BLKORCMT( LINE ) ) CYCLE

        !.......  Parse line into fields
        CALL PARSLINE( LINE, 53, SEGMENT )

        !.......  SCC mapping loop based on SCCXREF reference input file
        SCC = ADJUSTL( SEGMENT( 2 ) )
        CALL PADZERO( SCC )

        KK = 0
        NSCC = 0
        IF( SCCMAPFLAG ) THEN
            KK   = INDEX1( SCC, NSCCMAP, SCCMAPLIST( :,1 ) )
            IF( KK > 0 ) THEN
                NSCC = STR2INT( SCCMAPLIST( KK,3 ) )
            ELSE
                IF( EXCLSCCFLAG ) CYCLE            ! drop SCCs not listed in SCCXREF file
            END IF
        END IF

        DO JJ = 0, NSCC

            IF( JJ > 0 .AND. KK > 0 ) IREC = IREC + 1
            IF( SCCMAPFLAG .AND. KK > 0 ) SCC = SCCMAPLIST( KK+JJ,2 )

            !.......  Find SCC in inventory list
            SCCIDX = FINDC( SCC, NINVSCC, INVSCC )
            IF( SCCIDX .LE. 0 ) THEN
                WRITE( MESG, 94010 ) 'NOTE: Skipping line',&
                  IREC, 'of hourly speed file because SCC ' //&
                  SCC // ' is not in the inventory.'
                NWARN1 = NWARN1 + 1
                IF( NWARN1 < MXWARN ) CALL M3MESG( MESG )
                CYCLE
            END IF

            !.......  Find county in inventory list
            CFIP = ADJUSTL( SEGMENT( 1 ) )
            CALL PADZERO( CFIP )
            FIPIDX = FINDC( CFIP, NINVIFIP, INVCFIP )

            IF( FIPIDX .LE. 0 ) THEN
                WRITE( MESG, 94010 ) 'NOTE: Skipping line',&
                  IREC, 'of hourly speed file because FIPS code '//&
                  CFIP //' is not in the inventory.'
                NWARN2 = NWARN2 + 1
                IF( NWARN2 < MXWARN ) CALL M3MESG( MESG )
                CYCLE
            END IF

            !.......  Check weekday keyword
            KEYWORD = ADJUSTL( SEGMENT( 3 ) )
            IF( KEYWORD .NE. 'WEEKDAY' ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Missing expected ' //&
                  'keyword WEEKDAY at line', IREC,&
                  'of hourly speed file.'
                CALL M3MESG( MESG )
                CYCLE
            END IF

            !.......  Check weekday speed values
            DO J = 1, 24
                IF( .NOT. CHKREAL( SEGMENT( 3 + J ) ) ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 ) 'ERROR: Bad weekday ' //&
                      'hourly speed value at line', IREC,&
                      'of hourly speed file.'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

                SPDVAL = STR2REAL( SEGMENT( 3 + J ) )
                SPDPRO( FIPIDX, SCCIDX, 2, J ) = SPDVAL
            END DO

            !.......  Check weekend keyword
            KEYWORD = ADJUSTL( SEGMENT( 28 ) )
            IF( KEYWORD .NE. 'WEEKEND' ) THEN
                EFLAG = .TRUE.
                WRITE( MESG, 94010 ) 'ERROR: Missing expected ' //&
                  'keyword WEEKEND at line', IREC,&
                  'of hourly speed file.'
                CALL M3MESG( MESG )
                CYCLE
            END IF

            !.......  Check weekend speed values
            DO J = 1, 24
                IF( .NOT. CHKREAL( SEGMENT( 28 + J ) ) ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 ) 'ERROR: Bad weekend ' //&
                      'hourly speed value at line', IREC,&
                      'of hourly speed file.'
                    CALL M3MESG( MESG )
                    CYCLE
                END IF

                SPDVAL = STR2REAL( SEGMENT( 28 + J ) )
                SPDPRO( FIPIDX, SCCIDX, 1, J ) = SPDVAL
            END DO

        END DO            ! SCCXREF loop

    END DO

    CLOSE( SPDEV )

    IF( EFLAG ) THEN
        MESG = 'Problem found in hourly speed file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    RETURN

999 MESG = 'End of file'
    MESG = 'End of file reached unexpectedly. ' //&
           'Check format of SPDPRO' // CRLF() // BLANK5 //&
           'input file.'
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Formatted file I/O formats...... 93xxx

93000 FORMAT( A )

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDSPDPRO


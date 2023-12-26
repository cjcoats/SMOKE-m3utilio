
SUBROUTINE WRSRCGRPS( VNAME, JDATE, JTIME, INPUTFLAG, INPUTEMIS )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !     This subroutine writes the emissions for each source
    !     group and grid cell into the CMAQ inline emissions file.
    !     The first time the subroutine is called, it also writes
    !     the group metadata (group number, number of sources, dummy
    !     stack parameters, etc.) into the stack groups file.
    !
    !  PRECONDITIONS REQUIRED:
    !     Stack groups and inline emissions files opened for output
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created  7/2013 by C. Seppanen
    !       Version 11/2023 by CJC:  USE M3UTILIO, ".f90" source format,  and
    !       related changes
    !***************************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !                System
    ! File: @(#)$Id$
    !
    ! COPYRIGHT (C) 2013, Environmental Modeling for Policy Development
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
    !**************************************************************************

    USE M3UTILIO

    !.........  MODULES for public variables
    !.........  This module contains the major data structure and control flags
    USE MODMERGE, ONLY: EMGGRD, NSRCGRP, NSGOUTPUT, GRPCNT,         &
                        IGRPNUM, SGINLNNAME, SRCGRPNAME,            &
                        PFLAG, PVNAME, PVSDATE, PVSTIME, ISRCGRP

    !.........  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: NGRID, NCOLS, NROWS,                         &
                       GDTYP, GRDNM, P_ALP, P_BET, P_GAM,           &
                       XCENT, YCENT, XORIG, YORIG, XCELL, YCELL

    !.........  This module contains arrays for plume-in-grid and major sources
    USE MODELEV, ONLY: NGROUP, NELEVGRPS, EMELEVGRP,                &
                       ELEVSTKGRP, ELEVSRCGRP, ELEVSTKCNT, SGFIREFLAG

    USE MODGRDLIB

    IMPLICIT NONE

    !.........  INCLUDES:
    INCLUDE 'SETDECL.h90'   !  FileSetAPI variables and functions

    !...........   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: VNAME   ! variable name to output
    INTEGER     , INTENT (IN) :: JDATE   ! Julian date to output (YYYYDDD)
    INTEGER     , INTENT (IN) :: JTIME   ! time to output (HHMMSS)
    LOGICAL     , INTENT (IN) :: INPUTFLAG      ! indicate if input emissions need to be used
    REAL        , INTENT (IN) :: INPUTEMIS( * ) ! input emissions (optional)

    !...........   Local arrays
    INTEGER     INTDATA( NGROUP )              ! generic integer data
    REAL       REALDATA( NGROUP )              ! generic real data
    INTEGER     ISTACK ( NSGOUTPUT ) ! group number
    INTEGER     STKCNT ( NSGOUTPUT ) ! num. srcs per group
    INTEGER     ROW    ( NSGOUTPUT ) ! row number
    INTEGER     COL    ( NSGOUTPUT ) ! column number
    INTEGER     LMAJOR ( NSGOUTPUT ) ! major source flag
    INTEGER     LPING  ( NSGOUTPUT ) ! PinG source flag
    REAL        XLOCA  ( NSGOUTPUT ) ! y-location at center of grid cell
    REAL        YLOCA  ( NSGOUTPUT ) ! y-location at center of grid cell
    REAL        LAT    ( NSGOUTPUT ) ! latitude of YLOCA
    REAL        LONG   ( NSGOUTPUT ) ! longitude of XLOCA
    REAL        STKDM  ( NSGOUTPUT ) ! inside stack diameter
    REAL        STKHT  ( NSGOUTPUT ) ! stack height
    REAL        STKTK  ( NSGOUTPUT ) ! stack exit temperature
    REAL        STKVE  ( NSGOUTPUT ) ! stack exit velocity
    REAL        STKFLW ( NSGOUTPUT ) ! stack exit flow rate
    REAL        ACRES  ( NSGOUTPUT ) ! acres burned for a fire
    REAL        OUTEMIS( NSGOUTPUT ) ! output emissions

    !...........   Other local variables
    INTEGER          C, G, K, IDX  ! counters and indices
    INTEGER          IOS           ! i/o status
    INTEGER          ROWNUM        ! grid cell row
    INTEGER          COLNUM        ! grid cell column
    INTEGER          ELEVIDX       ! starting index for elevated sources

    REAL             XLOCACELL     ! x-location for grid cell
    REAL             YLOCACELL     ! y-location for grid cell

    LOGICAL, SAVE :: FIRSTTIME = .TRUE. ! true: first time routine called

    CHARACTER(300)   MESG         ! message buffer

    CHARACTER(16), PARAMETER :: PNAME = 'WRSRCGRPS' ! program name

    !**********************************************************************
    !   begin body of subroutine WRSRCGRPS

    IF( FIRSTTIME ) THEN

        !.............  Output stack groups file
        ISTACK = 0      ! array
        STKCNT = 0
        ROW    = 0
        COL    = 0
        LMAJOR = 0
        LPING  = 0
        XLOCA  = BADVAL3
        YLOCA  = BADVAL3
        LAT    = BADVAL3
        LONG   = BADVAL3

        !.............  Set dummy stack parameter arrays based on environment settings
        STKDM  = ENVREAL( 'SRCGRP_STKDM',  'Stack diameter',           0.1, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PNAME,0,0, 'Bad env vble "SRCGRP_STKDM"', 2 )
        END IF
        STKHT  = ENVREAL( 'SRCGRP_STKHT',  'Stack height',             0.1, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PNAME,0,0, 'Bad env vble "SRCGRP_STKHT"', 2 )
        END IF
        STKTK  = ENVREAL( 'SRCGRP_STKTK',  'Stack exit temperature', 273.0, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PNAME,0,0, 'Bad env vble "SRCGRP_STKTK"', 2 )
        END IF
        STKVE  = ENVREAL( 'SRCGRP_STKVE',  'Stack exit velocity',      0.1, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PNAME,0,0, 'Bad env vble "SRCGRP_STKVE"', 2 )
        END IF
        STKFLW = ENVREAL( 'SRCGRP_STKFLW', 'Stack exit flow rate',     0.1, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PNAME,0,0, 'Bad env vble "SRCGRP_STKFLW"', 2 )
        END IF

        K = 0
        DO C = 1, NGRID

            !.................  Determine row and column for current grid cell
            ROWNUM = C / NCOLS           ! integer math
            IF( MOD( C, NCOLS ) .GT. 0 ) ROWNUM = ROWNUM + 1
            COLNUM = C - ( ROWNUM-1 ) * NCOLS

            !.................  Calculate x and y-position at center of grid cell
            XLOCACELL = XORIG + ( COLNUM-1 ) * XCELL + 0.5 * XCELL
            YLOCACELL = YORIG + ( ROWNUM-1 ) * YCELL + 0.5 * YCELL

            DO G = 1, NSRCGRP

                !.....................  Skip missing values
                IF( GRPCNT( C, G ) == 0 ) CYCLE

                K = K + 1
                ISTACK( K ) = IGRPNUM( G )
                STKCNT( K ) = GRPCNT( C, G )
                ROW   ( K ) = ROWNUM
                COL   ( K ) = COLNUM
                XLOCA ( K ) = XLOCACELL
                YLOCA ( K ) = YLOCACELL
            END DO
        END DO

        !.............  Convert x and y grid cell locations to lat/lon
        LAT = YLOCA
        LONG = XLOCA
        CALL CONVRTLL( K, GDTYP, GRDNM,         &
                       P_ALP, P_BET, P_GAM,     &
                       XCENT, YCENT, LONG, LAT )

        !.............  Append data for elevated source groups
        IF( PFLAG ) THEN
            ELEVIDX = K + 1

            DO G = 1, NELEVGRPS
                IDX = ELEVIDX + G - 1
                ISTACK( IDX ) = IGRPNUM( ELEVSRCGRP( G ) )
            END DO

            STKCNT( ELEVIDX:NSGOUTPUT ) = ELEVSTKCNT

            CALL INT_READ3( PVNAME, 'ROW', 1, PVSDATE, PVSTIME, INTDATA )
            DO G = 1, NELEVGRPS
                IDX = ELEVIDX + G - 1
                ROW( IDX ) = INTDATA( ELEVSTKGRP( G ) )
            END DO

            CALL INT_READ3( PVNAME, 'COL', 1, PVSDATE, PVSTIME, INTDATA )
            DO G = 1, NELEVGRPS
                IDX = ELEVIDX + G - 1
                COL( IDX ) = INTDATA( ELEVSTKGRP( G ) )
            END DO

            CALL INT_READ3( PVNAME, 'LMAJOR', 1, PVSDATE, PVSTIME, INTDATA )
            DO G = 1, NELEVGRPS
                IDX = ELEVIDX + G - 1
                LMAJOR( IDX ) = INTDATA( ELEVSTKGRP( G ) )
            END DO

            CALL INT_READ3( PVNAME, 'LPING', 1, PVSDATE, PVSTIME, INTDATA )
            DO G = 1, NELEVGRPS
                IDX = ELEVIDX + G - 1
                LPING( IDX ) = INTDATA( ELEVSTKGRP( G ) )
            END DO

            CALL REAL_READ3( PVNAME, 'XLOCA', 1, PVSDATE, PVSTIME, REALDATA )
            DO G = 1, NELEVGRPS
                IDX = ELEVIDX + G - 1
                XLOCA( IDX ) = REALDATA( ELEVSTKGRP( G ) )
            END DO

            CALL REAL_READ3( PVNAME, 'YLOCA', 1, PVSDATE, PVSTIME, REALDATA )
            DO G = 1, NELEVGRPS
                IDX = ELEVIDX + G - 1
                YLOCA( IDX ) = REALDATA( ELEVSTKGRP( G ) )
            END DO

            CALL REAL_READ3( PVNAME, 'STKDM', 1, PVSDATE, PVSTIME, REALDATA )
            DO G = 1, NELEVGRPS
                IDX = ELEVIDX + G - 1
                STKDM( IDX ) = REALDATA( ELEVSTKGRP( G ) )
            END DO

            CALL REAL_READ3( PVNAME, 'STKHT', 1, PVSDATE, PVSTIME, REALDATA )
            DO G = 1, NELEVGRPS
                IDX = ELEVIDX + G - 1
                STKHT( IDX ) = REALDATA( ELEVSTKGRP( G ) )
            END DO

            CALL REAL_READ3( PVNAME, 'STKTK', 1, PVSDATE, PVSTIME, REALDATA )
            DO G = 1, NELEVGRPS
                IDX = ELEVIDX + G - 1
                STKTK( IDX ) = REALDATA( ELEVSTKGRP( G ) )
            END DO

            CALL REAL_READ3( PVNAME, 'STKVE', 1, PVSDATE, PVSTIME, REALDATA )
            DO G = 1, NELEVGRPS
                IDX = ELEVIDX + G - 1
                STKVE( IDX ) = REALDATA( ELEVSTKGRP( G ) )
            END DO

            CALL REAL_READ3( PVNAME, 'STKFLW', 1, PVSDATE, PVSTIME, REALDATA )
            DO G = 1, NELEVGRPS
                IDX = ELEVIDX + G - 1
                STKFLW( IDX ) = REALDATA( ELEVSTKGRP( G ) )
            END DO

            !.................  If lat/lon is in existing stack groups file, append it
            IF( READ3( PVNAME, 'LATITUDE', 1, PVSDATE, PVSTIME, REALDATA ) ) THEN

                DO G = 1, NELEVGRPS
                    IDX = ELEVIDX + G - 1
                    LAT( IDX ) = REALDATA( ELEVSTKGRP( G ) )
                END DO

                !.....................  Assume longitude is available if latitude was
                CALL REAL_READ3( PVNAME, 'LONGITUDE', 1, PVSDATE, PVSTIME, REALDATA )
                DO G = 1, NELEVGRPS
                    IDX = ELEVIDX + G - 1
                    LONG( IDX ) = REALDATA( ELEVSTKGRP( G ) )
                END DO

            ELSE

                !.....................  Otherwise, convert x and y grid cell locations
                LAT( ELEVIDX:NSGOUTPUT )  = YLOCA( ELEVIDX:NSGOUTPUT )
                LONG( ELEVIDX:NSGOUTPUT ) = XLOCA( ELEVIDX:NSGOUTPUT )
                CALL CONVRTLL( NSGOUTPUT-ELEVIDX+1, GDTYP, GRDNM,       &
                               P_ALP, P_BET, P_GAM, XCENT, YCENT,       &
                               LONG( ELEVIDX:NSGOUTPUT ),               &
                               LAT(  ELEVIDX:NSGOUTPUT ) )

            END IF

        END IF

        CALL INT_WRITE3( SRCGRPNAME, 'IGROUP', JDATE, JTIME, ISTACK )
        CALL INT_WRITE3( SRCGRPNAME, 'GRPCNT', JDATE, JTIME, STKCNT )
        CALL INT_WRITE3( SRCGRPNAME, 'ROW',    JDATE, JTIME, ROW )
        CALL INT_WRITE3( SRCGRPNAME, 'COL',    JDATE, JTIME, COL )
        CALL INT_WRITE3( SRCGRPNAME, 'LMAJOR', JDATE, JTIME, LMAJOR )
        CALL INT_WRITE3( SRCGRPNAME, 'LPING',  JDATE, JTIME, LPING )

        CALL REAL_WRITE3( SRCGRPNAME, 'XLOCA', JDATE, JTIME, XLOCA )
        CALL REAL_WRITE3( SRCGRPNAME, 'YLOCA', JDATE, JTIME, YLOCA )
        CALL REAL_WRITE3( SRCGRPNAME, 'LATITUDE',  JDATE, JTIME, LAT )
        CALL REAL_WRITE3( SRCGRPNAME, 'LONGITUDE', JDATE, JTIME, LONG )
        CALL REAL_WRITE3( SRCGRPNAME, 'STKDM', JDATE, JTIME, STKDM )
        CALL REAL_WRITE3( SRCGRPNAME, 'STKHT', JDATE, JTIME, STKHT )
        CALL REAL_WRITE3( SRCGRPNAME, 'STKTK', JDATE, JTIME, STKTK )
        CALL REAL_WRITE3( SRCGRPNAME, 'STKVE', JDATE, JTIME, STKVE )
        CALL REAL_WRITE3( SRCGRPNAME, 'STKFLW', JDATE, JTIME, STKFLW )

    !.................  If acres burned is in existing stack groups file, add to output
        IF( SGFIREFLAG ) THEN

            ACRES = 0.

            CALL REAL_READ3( PVNAME, 'ACRESBURNED', 1, PVSDATE, PVSTIME, REALDATA )
            DO G = 1, NELEVGRPS
                IDX = ELEVIDX + G - 1
                ACRES( IDX ) = REALDATA( ELEVSTKGRP( G ) )
            END DO
            CALL REAL_WRITE3( SRCGRPNAME, 'ACRESBURNED', JDATE, JTIME, ACRES )

        END IF

        FIRSTTIME = .FALSE.
    END IF

    OUTEMIS = 0.  ! array

    K = 0
    DO C = 1, NGRID
        DO G = 1, NSRCGRP

            !.................  Skip missing values
            IF( GRPCNT( C, G ) == 0 ) CYCLE

            K = K + 1
            OUTEMIS( K ) = EMGGRD( C, G )

            IF( INPUTFLAG ) THEN
                OUTEMIS( K ) = OUTEMIS( K ) + INPUTEMIS( K )
            END IF
        END DO
    END DO

    !.........  Append emissions for elevated source groups
    IF( PFLAG ) THEN
        DO G = 1, NELEVGRPS

            K = K + 1
            OUTEMIS( K ) = EMELEVGRP( G )

        END DO
    END IF

    IF( .NOT. WRITESET( SGINLNNAME, VNAME, ALLFILES,    &
                        JDATE, JTIME, OUTEMIS ) ) THEN

        MESG = 'Could not write "' // VNAME // '" to file "' // SGINLNNAME // '"'
        CALL M3EXIT( PNAME, JDATE, JTIME, MESG, 2 )

    END IF

    RETURN

    !*****************  INTERNAL SUBPROGRAMS  ******************************

CONTAINS

    !.............  This internal subprogram reads real data from an
    !               I/O API file, and aborts if not successful.
    SUBROUTINE REAL_READ3( FILNAM, VARNAM, LAYER, RDATE, RTIME, REALBUF )

    !.............  Subprogram arguments
        CHARACTER(*) FILNAM       ! logical file name
        CHARACTER(*) VARNAM       ! variable name
        INTEGER      LAYER        ! layer number
        INTEGER      RDATE        ! read Julian date
        INTEGER      RTIME        ! read time
        REAL         REALBUF(*)   ! real data buffer

        !----------------------------------------------------------------------

        IF ( GET_VTYPE( FILNAM, VARNAM, RDATE, RTIME ) .NE. M3REAL ) THEN
            MESG = 'Type for variable "' // TRIM( VARNAM ) //       &
                   '" in "'// TRIM( FILNAM ) // '" not M3REAL'
            CALL M3EXIT( PNAME, RDATE, RTIME, MESG, 2 )
        ELSE IF ( .NOT. READ3( FILNAM, VARNAM, LAYER,               &
                               RDATE, RTIME, REALBUF ) ) THEN

            MESG = 'Could not read "' // TRIM( VARNAM ) //          &
                   '" from file "' // TRIM( FILNAM ) // '"'
            CALL M3EXIT( PNAME, RDATE, RTIME, MESG, 2 )

        END IF

        RETURN

    END SUBROUTINE REAL_READ3

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------

    !.............  This internal subprogram reads integer data from an
    !               I/O API file, and aborts if not successful.
    SUBROUTINE INT_READ3( FILNAM, VARNAM, LAYER, RDATE, RTIME, INTBUF )

    !.............  Subprogram arguments
        CHARACTER(*) FILNAM       ! logical file name
        CHARACTER(*) VARNAM       ! variable name
        INTEGER      LAYER        ! layer number
        INTEGER      RDATE        ! read Julian date
        INTEGER      RTIME        ! read time
        INTEGER      INTBUF(*)    ! integer data buffer

    !----------------------------------------------------------------------

        IF ( GET_VTYPE( FILNAM, VARNAM, RDATE, RTIME ) .NE. M3INT ) THEN

            MESG = 'Type for variable "' // TRIM( VARNAM ) //   &
                   '" in "'              // TRIM( FILNAM ) //   &
                   '" not M3INT'
            CALL M3EXIT( PNAME, RDATE, RTIME, MESG, 2 )

        ELSE IF ( .NOT. READ3( FILNAM, VARNAM, LAYER,           &
                               RDATE, RTIME, INTBUF ) ) THEN

            MESG = 'Could not read "' // TRIM( VARNAM ) //      &
                   '" from file "'    // TRIM( FILNAM ) // '"'
            CALL M3EXIT( PNAME, RDATE, RTIME, MESG, 2 )

        END IF

        RETURN

    END SUBROUTINE INT_READ3

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------

    !.............  This internal subprogram writes real data to an
    !               I/O API file, and aborts if not successful.
    SUBROUTINE REAL_WRITE3( FILNAM, VARNAM,             &
                            WDATE, WTIME, REALBUF )

    !.............  Subprogram arguments
        CHARACTER(*) FILNAM       ! logical file name
        CHARACTER(*) VARNAM       ! variable name
        INTEGER      WDATE        ! write Julian date
        INTEGER      WTIME        ! write time
        REAL         REALBUF(*)   ! real data buffer

    !----------------------------------------------------------------------

        IF ( GET_VTYPE( FILNAM, VARNAM, WDATE, WTIME ) .NE. M3REAL ) THEN

            MESG = 'Type for variable "' // TRIM( VARNAM ) //   &
                   '" in "'              // TRIM( FILNAM ) //   &
                   '" not M3REAL'
            CALL M3EXIT( PNAME, WDATE, WTIME, MESG, 2 )

        ELSE IF ( .NOT. WRITE3( FILNAM, VARNAM, WDATE, WTIME, REALBUF ) ) THEN

            MESG = 'Could not write "' // TRIM( VARNAM ) //         &
                   '" to file "'       // TRIM( FILNAM ) // '"'
            CALL M3EXIT( PNAME, WDATE, WTIME, MESG, 2 )

        END IF

        RETURN

    END SUBROUTINE REAL_WRITE3

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------

    !.............  This internal subprogram writes integer data to an
    !               I/O API file, and aborts if not successful.
    SUBROUTINE INT_WRITE3( FILNAM, VARNAM, WDATE, WTIME, INTBUF )

    !.............  Subprogram arguments
        CHARACTER(*) FILNAM       ! logical file name
        CHARACTER(*) VARNAM       ! variable name
        INTEGER      WDATE        ! write Julian date
        INTEGER      WTIME        ! write time
        INTEGER      INTBUF(*)    ! integer data buffer

    !----------------------------------------------------------------------

        IF ( GET_VTYPE( FILNAM, VARNAM, WDATE, WTIME ) .NE. M3INT ) THEN

            MESG = 'Type for variable "' // TRIM( VARNAM ) //   &
                   '" in "'              // TRIM( FILNAM ) //   &
                   '" not M3INT'
            CALL M3EXIT( PNAME, WDATE, WTIME, MESG, 2 )

        ELSE IF ( .NOT. WRITE3( FILNAM, VARNAM, WDATE, WTIME, INTBUF ) ) THEN

            MESG = 'Could not write "' // TRIM( VARNAM ) //     &
                   '" to file "'       // TRIM( FILNAM ) // '"'
            CALL M3EXIT( PNAME, WDATE, WTIME, MESG, 2 )

        END IF

        RETURN

    END SUBROUTINE INT_WRITE3

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------

    !.............  This internal subprogram writes integer data to an
    !               I/O API file, and aborts if not successful.
    INTEGER FUNCTION GET_VTYPE( FILNAM, VARNAM, JDATE, JTIME )

        CHARACTER(LEN=*), INTENT(IN   ) :: FILNAM, VARNAM
        INTEGER         , INTENT(IN   ) :: JDATE, JTIME

        INTEGER     K

        IF ( .NOT.DESC3( FILNAM ) ) THEN
            MESG = 'File "' // TRIM( FILNAM ) // '" not available'
            CALL M3EXIT( PNAME, JDATE, JTIME, MESG, 2 )
        END IF

        K = INDEX1( VARNAM, NVARS3D, VNAME3D )
        IF ( K .LE. 0 ) THEN
            MESG = 'Variable "'          // TRIM( VARNAM ) //&
                   '" not available in "'// TRIM( FILNAM ) // '"'
            CALL M3EXIT( PNAME, JDATE, JTIME, MESG, 2 )
        END IF

        GET_VTYPE = VTYPE3D( K )
        RETURN

    END FUNCTION GET_VTYPE

END SUBROUTINE WRSRCGRPS

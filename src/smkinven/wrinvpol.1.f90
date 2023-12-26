
SUBROUTINE WRINVPOL( CATEGORY, INPATH, INNAM, NREC,&
&                     NPVAR, VBUF, SRCID, POLBUF, EFLAG )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine writes the inventory pollutant-based variables
!      and/or VMT to the I/O API files in sparse storage format
!
!  PRECONDITIONS REQUIRED:
!      Number of sources NSRC defined correctly
!      Pollutant count IPCNT and output names POLNAM defined correctly
!      Output array POLVAL populated
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!      Subroutines: I/O API subroutine, BLDENAMS
!      Functions: I/O API subroutine
!
!  REVISION  HISTORY:
!       Created 11/98 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!***************************************************************************
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

!...........   MODULES for public variables
    USE MODFILESET

    IMPLICIT NONE

!...........   INCLUDES

    INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
    INCLUDE 'SETDECL.EXT'   !  FileSetAPI function declarations

!.............  Subroutine arguments
    CHARACTER(*), INTENT (IN) :: CATEGORY        ! source category
    CHARACTER(*), INTENT (IN) :: INPATH          ! path for output file
    CHARACTER(*), INTENT (IN) :: INNAM           ! name for output file
    INTEGER     , INTENT (IN) :: NREC            ! size of output arrays
    INTEGER     , INTENT (IN) :: NPVAR           ! number variables per pol/act
    CHARACTER(*), INTENT (IN) :: VBUF            ! names of pols/act
    INTEGER     , INTENT (IN) :: SRCID( NREC )   ! source IDs
    REAL        , INTENT (IN) :: POLBUF( NREC,NPVAR )  ! pol or act data
    LOGICAL     , INTENT(OUT) :: EFLAG           ! true: error found

!.............  Arrays for sparse output
    INTEGER     INTBUF ( NREC )

!.............  Local variables
    INTEGER      J, L, N
    INTEGER      IOS      ! i/o status

    LOGICAL       :: MSGFLAG  = .FALSE.  ! true: need to output error message

    CHARACTER(16)      :: FNAME = 'IOAPI_DAT' ! tmp logical name for outputs
    CHARACTER(256)     :: MESG                ! message buffer
    CHARACTER(PHYLEN3) :: FPHYS = ' '         ! full physical file name

    CHARACTER(16) :: PROGNAME = 'WRINVPOL' !  program name

!***********************************************************************
!   begin body of program WRINVPOL

    FPHYS = TRIM( INPATH ) // TRIM( INNAM )

!.........  Set output logical file name
    IF( .NOT. SETENVVAR( FNAME, FPHYS ) ) THEN
        EFLAG = .TRUE.
        MESG = 'ERROR: Could not set logical file name ' //&
        &       'for "' // TRIM( VBUF ) // '" file:'//&
        &       CRLF()// BLANK10// TRIM( FPHYS )
        CALL M3MSG2( MESG )

!.........  Open I/O API file
    ELSE IF( .NOT. OPENSET( FNAME, FSNEW3, PROGNAME )) THEN
        EFLAG = .TRUE.
        MESG = 'ERROR: Could not open I/O API "' // TRIM( VBUF )//&
        &       '" file for file name:' // CRLF() // BLANK10 //&
        &       TRIM( FPHYS )
        CALL M3MSG2( MESG )

    END IF

!.........  Write source ID index
    IF( .NOT. WRITESET( FNAME, 'SRCID', ALLFILES,&
    &                    0, 0, SRCID               )) THEN
        EFLAG = .TRUE.
        MESG = 'Could not write "SRCID" to file:'//&
        &        CRLF()//BLANK10// TRIM( FPHYS )
        CALL M3MSG2( MESG )
    END IF

!.........  Write all variables for each pollutant or activity
    DO J = 1, NPVAR

        MSGFLAG = .FALSE.
        IF( VTYPESET( J+1 ) .EQ. M3REAL ) THEN
            IF( .NOT. WRITESET( FNAME, VNAMESET( J+1 ),&
            &    ALLFILES, 0, 0, POLBUF( 1,J )   ) ) MSGFLAG = .TRUE.

        ELSE
            INTBUF( 1:NREC ) = INT( POLBUF( 1:NREC, J ) )  ! array
            IF( .NOT. WRITESET( FNAME, VNAMESET( J+1 ),&
            &    ALLFILES, 0, 0, INTBUF          ) ) MSGFLAG = .TRUE.
        END IF

        IF( MSGFLAG ) THEN
            EFLAG = .TRUE.
            MESG = 'Could not write "'// TRIM( VNAMESET(J+1) )&
            &       // '" to file:'// CRLF()// BLANK10//&
            &       TRIM( FPHYS )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

    END DO

!........  Close output file for this variable
    IF( .NOT. CLOSESET( FNAME ) ) THEN
        MESG = 'Could not close file:'//CRLF()//BLANK10//&
        &               TRIM( FPHYS )
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    RETURN

END SUBROUTINE WRINVPOL

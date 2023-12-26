
SUBROUTINE WMRGEMIS( VNAME, JDATE, JTIME )

!***********************************************************************
!  subroutine WMRGEMIS body starts at line
!
!  DESCRIPTION:
!      The purpose of this subroutine is to write out all NetCDF files
!      coming from the merge program.  It is expected that this routine
!      will be called for each time step.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 3/99 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
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

!.........  MODULES for public variables
!.........  This module contains the major data structure and control flags
    USE MODMERGE, ONLY: AFLAG, BFLAG, MFLAG, PFLAG, XFLAG, SFLAG,&
    &                    ANMSPC, BNMSPC, MNMSPC, PNMSPC,&
    &                    AEMNAM, BEMNAM, MEMNAM, PEMNAM,&
    &                    ANIPOL,         MNIPPA, PNIPOL,&
    &                    AEINAM,         MEANAM, PEINAM,&
    &                    AONAME, BONAME, MONAME, PONAME, TONAME,&
    &                    AEMGRD, BEMGRD, MEMGRD, PEMGRD, TEMGRD,&
    &                    LGRDOUT, PINGFLAG, PINGNAME,&
    &                    INLINEFLAG, INLINENAME

!.........  This module contains arrays for plume-in-grid and major sources
    USE MODELEV, ONLY: PGRPEMIS,ELEVEMIS

    IMPLICIT NONE

!.........  INCLUDES:

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
    INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

!.........  SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: VNAME   ! variable name to output
    INTEGER     , INTENT (IN) :: JDATE   ! Julian date to output (YYYYDDD)
    INTEGER     , INTENT (IN) :: JTIME   ! time to output (HHMMSS)

!.........  Other local variables
    INTEGER         J

    LOGICAL      :: AOUTFLAG = .FALSE.  ! true: output area sources
    LOGICAL      :: BOUTFLAG = .FALSE.  ! true: output biogenic sources
    LOGICAL      :: MOUTFLAG = .FALSE.  ! true: output mobile sources
    LOGICAL      :: POUTFLAG = .FALSE.  ! true: output point sources
    LOGICAL      :: IOUTFLAG = .FALSE.  ! true: output point Ping or inline

    CHARACTER(NAMLEN3) FILNAM       ! tmp logical file name
    CHARACTER(300)     MESG         ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'WMRGEMIS' ! program name

!***********************************************************************
!   begin body of subroutine WMRGEMIS

!.........  Initialize output flags
    AOUTFLAG = .FALSE.
    BOUTFLAG = .FALSE.
    MOUTFLAG = .FALSE.
    POUTFLAG = .FALSE.
    IOUTFLAG = .FALSE.

!.........  Determine the source categories that are valid for output
!.........  If the subroutine call is for speciated output, use different
!           indicator arrays for determining output or not.
    IF( SFLAG ) THEN

        IF( LGRDOUT .AND. AFLAG ) THEN
            J = INDEX1( VNAME, ANMSPC, AEMNAM )
            AOUTFLAG = ( J .GT. 0 )
        END IF

        IF( LGRDOUT .AND. BFLAG ) THEN
            J = INDEX1( VNAME, BNMSPC, BEMNAM )
            BOUTFLAG = ( J .GT. 0 )
        END IF

        IF( LGRDOUT .AND. MFLAG ) THEN
            J = INDEX1( VNAME, MNMSPC, MEMNAM )
            MOUTFLAG = ( J .GT. 0 )
        END IF

        IF( LGRDOUT .AND. PFLAG ) THEN
            J = INDEX1( VNAME, PNMSPC, PEMNAM )
            POUTFLAG = ( J .GT. 0 )
        END IF

        IF( ( PINGFLAG .OR. INLINEFLAG ) .AND. PFLAG ) THEN
            J = INDEX1( VNAME, PNMSPC, PEMNAM )
            IOUTFLAG = ( J .GT. 0 )
        END IF

!.........  Non-speciated (it's not possible to have biogenics w/o speciation)
    ELSE

        IF( LGRDOUT .AND. AFLAG ) THEN
            J = INDEX1( VNAME, ANIPOL, AEINAM )
            AOUTFLAG = ( J .GT. 0 )
        END IF

        IF( LGRDOUT .AND. MFLAG ) THEN
            J = INDEX1( VNAME, MNIPPA, MEANAM )
            MOUTFLAG = ( J .GT. 0 )
        END IF

        IF( LGRDOUT .AND. PFLAG ) THEN
            J = INDEX1( VNAME, PNIPOL, PEINAM )
            POUTFLAG = ( J .GT. 0 )
        END IF

        IF( ( PINGFLAG .OR. INLINEFLAG) .AND. PFLAG ) THEN
            J = INDEX1( VNAME, PNIPOL, PEINAM )
            IOUTFLAG = ( J .GT. 0 )
        END IF

    END IF

!.........  For area sources, output file...
    IF( AOUTFLAG ) CALL SAFE_WRITE3( AONAME, AEMGRD )

!.........  For biogenic, output file...
    IF( BOUTFLAG )&
    &    CALL SAFE_WRITE3( BONAME, BEMGRD )

!.........  For mobile sources, output file...
    IF( MOUTFLAG ) CALL SAFE_WRITE3( MONAME, MEMGRD )

!.........  For point sources, output file...
    IF( POUTFLAG ) CALL SAFE_WRITE3( PONAME, PEMGRD )

!.........  For plume-in-grid, output file...
    IF( IOUTFLAG .AND. PINGFLAG )&
    &    CALL SAFE_WRITE3( PINGNAME, PGRPEMIS )

!.........  For plume-in-grid, output file...
    IF( IOUTFLAG .AND. INLINEFLAG )&
    &    CALL SAFE_WRITE3( INLINENAME, ELEVEMIS )

!.........  For multiple source categories, output totals file...
    IF( LGRDOUT .AND. XFLAG ) CALL SAFE_WRITE3( TONAME, TEMGRD )

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats.............94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

!*****************  INTERNAL SUBPROGRAMS  ******************************

CONTAINS

!.............  This internal subprogram uses WRITE3 and exits gracefully
!               if a write error occurred
    SUBROUTINE SAFE_WRITE3( FILNAM, EMDATA )

        INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

        CHARACTER(*), INTENT (IN) :: FILNAM
        REAL        , INTENT (IN) :: EMDATA( * )

!.......................................................................

        IF( .NOT. WRITESET( FILNAM, VNAME, ALLFILES,&
        &                    JDATE, JTIME, EMDATA ) ) THEN

            MESG = 'Could not write "' // VNAME //&
            &       '" to file "'// FILNAM // '"'
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

        END IF

        RETURN

    END SUBROUTINE SAFE_WRITE3

END SUBROUTINE WMRGEMIS

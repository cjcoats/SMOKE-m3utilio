
SUBROUTINE WRRMAT( NSRC, NMSPC, FDEV, FILNAM, INDX, REPEM,&
&                   PRJFAC, MKTPEN, RMTX, CSCC,&
&                   SPROF, VNAMES )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine writes the reactivity matrices
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 3/99 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!****************************************************************************/
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

    IMPLICIT NONE

!...........   INCLUDES

    INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters
    INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

!.........  Subroutine arguments and their descriptions:
    INTEGER     , INTENT (IN) :: NSRC           ! no. of source
    INTEGER     , INTENT (IN) :: NMSPC          ! no. model species
    INTEGER     , INTENT (IN) :: FDEV           ! supplement file unit no.
    CHARACTER(*), INTENT (IN) :: FILNAM         ! I/O API file name
    INTEGER     , INTENT (IN) :: INDX  ( NSRC ) ! pollutant count
    REAL        , INTENT (IN) :: REPEM ( NSRC ) ! replacement emissions
    REAL        , INTENT (IN) :: PRJFAC( NSRC ) ! projection factors
    REAL        , INTENT (IN) :: MKTPEN( NSRC ) ! market penetration
    REAL        , INTENT (IN) :: RMTX  ( NSRC,NMSPC ) ! speciation data
    CHARACTER(*), INTENT (IN) :: CSCC  ( NSRC ) ! reactivity SCCs
    CHARACTER(*), INTENT (IN) :: SPROF ( NSRC ) ! speciation profile IDs
    CHARACTER(*), INTENT (IN) :: VNAMES( NMSPC ) ! variable names for spcs

!...........   Other local variables
    INTEGER                 I, L, L2, S

    LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first time subroutine called

    CHARACTER(NAMLEN3) VAR     !  tmp variable names

    CHARACTER(100) OUTFMT           !  format buffer
    CHARACTER(300) MESG             !  message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'WRRMAT' !  program name

!***********************************************************************
!   begin body of program WRRMAT

!.........  The first time this routine is called, write out the supplement file
    IF( FIRSTIME ) THEN

!.............  Generate format
        WRITE( OUTFMT, 94020 ) SCCLEN3, SPNLEN3

!.............  Write out header
        WRITE( FDEV, 93010 ) 3, OUTFMT( 1:LEN_TRIM( OUTFMT ) )
        WRITE( FDEV, 93000 ) 'SRCID   Source index'
        WRITE( FDEV, 93000 ) 'CSCC    Source category code '
        WRITE( FDEV, 93000 ) 'SPROF   Speciation profile number'

!.............  Write out records for the ASCII supplement file
        DO S = 1, NSRC
            WRITE( FDEV, OUTFMT ) INDX( S ), CSCC( S ), SPROF( S )
        END DO

        FIRSTIME = .FALSE.

    END IF

!.........  Initialize message to use in case there is an error

    MESG = 'Problem writing to output file "' // TRIM( FILNAM )//'"'

    L = LEN_TRIM( MESG )

!.........  Write the I/O API variables for the non-speciation data

    IF( .NOT. WRITESET( FILNAM,'SRCID',ALLFILES, 0,0,INDX ) ) THEN
        MESG = MESG( 1:L ) // ' for variable "SRCID"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF( .NOT. WRITESET( FILNAM,'REPEMIS',ALLFILES,0,0,REPEM ) ) THEN
        MESG = MESG( 1:L ) // ' for variable "REPEMIS"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF( .NOT. WRITESET( FILNAM,'PRJFAC',ALLFILES,0,0,PRJFAC ) ) THEN
        MESG = MESG( 1:L ) // ' for variable "PRJFAC"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF( .NOT. WRITESET( FILNAM,'MKTPEN',ALLFILES,0,0,MKTPEN ) ) THEN
        MESG = MESG( 1:L ) // ' for variable "MKTPEN"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Write the I/O API variables for the speciation factors

    DO I = 1, NMSPC

        IF( .NOT. WRITESET( FILNAM, VNAMES( I ), ALLFILES,&
        &                    0, 0, RMTX( 1,I )     ) ) THEN
            MESG = MESG( 1:L  ) // ' for variable "' //&
            &       TRIM( VNAMES( I ) ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

    END DO

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

93010 FORMAT( I2, A )

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

94020 FORMAT( '(I6,1X,A', I2.2, ',1X,A', I2.2, ')' )

END SUBROUTINE WRRMAT

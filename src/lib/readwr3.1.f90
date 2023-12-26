
SUBROUTINE READWR3( INFILE, OUTFILE, VNAME, LAYSVAL,&
&                    JDATE, JTIME, VTYPE, NDIM, STATUS )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      Reads variable from INFILE and write variable to OUTFILE
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Version ??/???? by ???
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
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

    IMPLICIT NONE

!...........   SUBROUTINE ARGUMENTS
    CHARACTER(*)   INFILE           ! Name of file being read
    CHARACTER(*)   OUTFILE          ! Name of file being written
    CHARACTER(*)   VNAME            ! Variable name being read/written
    INTEGER        LAYSVAL          ! layer number or value of ALLAYS3
    INTEGER        JDATE            ! Julian date being read/written
    INTEGER        JTIME            ! Julian time being read/written
    INTEGER        VTYPE            ! Integer code for variable type
    INTEGER        NDIM             ! Dimension of var being read/written
    INTEGER        STATUS           ! Exit status

!...........   Allocatable arrays
    INTEGER         INTVAL ( NDIM )  !  Integer value
    REAL            REALVAL( NDIM )  !  Real value

!...........   Other local variables
    CHARACTER(300)  MESG

    CHARACTER(16) :: PROGNAME = 'READWR3' ! program name

!***********************************************************************
!   begin body of subroutine READWR3

    STATUS = 0

    IF( VTYPE .EQ. M3INT ) THEN

        IF( .NOT. READ3( INFILE, VNAME, LAYSVAL,&
        &    JDATE, JTIME, INTVAL ) ) THEN
            STATUS = 1
            MESG = 'ERROR: Could not read "' //&
            &       TRIM( VNAME ) // '" from file.'
            CALL M3MESG( MESG )
        END IF

        IF( .NOT. WRITE3( OUTFILE, VNAME, LAYSVAL,&
        &    JDATE, JTIME, INTVAL ) ) THEN
            STATUS = 1
            MESG = 'ERROR: Could not write "' //&
            &       TRIM( VNAME ) // '" from file.'
            CALL M3MESG( MESG )
        END IF

    ELSEIF( VTYPE .EQ. M3REAL ) THEN

        IF( .NOT. READ3( INFILE, VNAME, LAYSVAL,&
        &    JDATE, JTIME, REALVAL ) ) THEN
            STATUS = 1
            MESG = 'ERROR: Could not read "' //&
            &       TRIM( VNAME ) // '" from file.'
            CALL M3MESG( MESG )
        END IF

        IF( .NOT. WRITE3( OUTFILE, VNAME, LAYSVAL,&
        &    JDATE, JTIME, REALVAL ) ) THEN
            STATUS = 1
            MESG = 'ERROR: Could not write "' //&
            &       TRIM( VNAME ) // '" from file.'
            CALL M3MESG( MESG )
        END IF

    ENDIF

    RETURN

END SUBROUTINE READWR3


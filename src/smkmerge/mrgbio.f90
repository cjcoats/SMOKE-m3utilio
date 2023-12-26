
SUBROUTINE MRGBIO( VNAME, FNAME, JDATE, JTIME, NGRID,    &
                   UNITFAC, BIOARR, ALLARR )

    !***********************************************************************
    !  subroutine MRGBIO body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine reads the
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 11/99 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***********************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !       System
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

    IMPLICIT NONE

    !.....   INCLUDES:

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.....  SUBROUTINE ARGUMENTS

    CHARACTER(*), INTENT  (IN) :: VNAME           ! variable name to read
    CHARACTER(*), INTENT  (IN) :: FNAME           ! file name to read
    INTEGER     , INTENT  (IN) :: JDATE           ! Julian date
    INTEGER     , INTENT  (IN) :: JTIME           ! time
    INTEGER     , INTENT  (IN) :: NGRID           ! no. grid cells
    REAL        , INTENT  (IN) :: UNITFAC         ! units conv factor
    REAL        , INTENT (OUT) :: BIOARR( NGRID )     ! biogenic emissions
    REAL        , INTENT (OUT) :: ALLARR( NGRID )     ! merged emissions

    !.....   Other local variables
    INTEGER         C    !  counters and indices

    CHARACTER(300)  MESG     ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'MRGBIO'     ! program name

    !***********************************************************************
    !   begin body of subroutine MRGBIO

    !.....  Read biogenic emissions
    IF( .NOT. READ3( FNAME, VNAME, ALLAYS3,JDATE,JTIME, BIOARR      ) ) THEN
        MESG = 'Could not read "'// TRIM( VNAME ) // '" from ' // FNAME
        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

    END IF       ! if read3() failed

    !.....  Update all emissions data array with biogenics
    DO C = 1, NGRID
        BIOARR( C ) = BIOARR( C ) * UNITFAC
        ALLARR( C ) = ALLARR( C ) + BIOARR( C )
    END DO

    RETURN

END SUBROUTINE MRGBIO

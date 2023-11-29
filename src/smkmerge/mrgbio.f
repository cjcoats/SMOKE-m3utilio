
        SUBROUTINE MRGBIO( VNAME, FNAME, JDATE, JTIME, NGRID,
     &                     UNITFAC, BIOARR, ALLARR )

C***********************************************************************
C  subroutine MRGBIO body starts at line
C
C  DESCRIPTION:
C      This subroutine reads the
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 11/99 by M. Houyoux
C       Version 11/2023 by CJC:  USE M3UTILIO and related changes
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$
C
C****************************************************************************
        USE M3UTILIO

        IMPLICIT NONE

C...........   INCLUDES:

        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........  SUBROUTINE ARGUMENTS

        CHARACTER(*), INTENT  (IN) :: VNAME           ! variable name to read
        CHARACTER(*), INTENT  (IN) :: FNAME           ! file name to read
        INTEGER     , INTENT  (IN) :: JDATE           ! Julian date
        INTEGER     , INTENT  (IN) :: JTIME           ! time
        INTEGER     , INTENT  (IN) :: NGRID           ! no. grid cells
        REAL        , INTENT  (IN) :: UNITFAC         ! units conv factor
        REAL        , INTENT (OUT) :: BIOARR( NGRID ) ! biogenic emissions
        REAL        , INTENT (OUT) :: ALLARR( NGRID ) ! merged emissions

C...........   Other local variables
        INTEGER         C    !  counters and indices

        CHARACTER(300)         :: MESG ! message buffer

        CHARACTER(16) :: PROGNAME = 'MRGBIO' ! program name

C***********************************************************************
C   begin body of subroutine MRGBIO

C.........  Read biogenic emissions
        IF( .NOT. READ3( FNAME, VNAME, ALLAYS3,
     &                   JDATE, JTIME, BIOARR      ) ) THEN
            MESG = 'Could not read "'// TRIM( VNAME ) //
     &             '" from ' // FNAME
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

        END IF   ! if read3() failed

C.........  Update all emissions data array with biogenics
        DO C = 1, NGRID
            BIOARR( C ) = BIOARR( C ) * UNITFAC
            ALLARR( C ) = ALLARR( C ) + BIOARR( C )
        END DO

        RETURN

        END SUBROUTINE MRGBIO

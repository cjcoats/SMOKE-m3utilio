
        SUBROUTINE PROCINVSRCS( NRAWSRCS )

C**************************************************************************
C  subroutine body starts at line 114
C
C  DESCRIPTION:
C      This subroutine sorts and stores the unique source information.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C      Created 2/03 by C. Seppanen based on procinven.f
C
C       Version 11/2023 by CJC:  USE M3UTILIO and related changes
C**************************************************************************
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
C***************************************************************************
        USE M3UTILIO

C...........   MODULES for public variables
C...........   This module is the inventory arrays
        USE MODSOURC, ONLY: SRCIDA, CSOURCA, CSOURC, CIFIP,
     &                      CSCC, XLOCA, YLOCA, CELLID, IRCLAS,
     &                      IVTYPE, CLINK, CVTYPE

C.........  This module contains the information about the source category
        USE MODINFO, ONLY: CATEGORY, NSRC

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters

C...........   SUBROUTINE ARGUMENTS
        INTEGER , INTENT (IN) :: NRAWSRCS ! no. raw srcs

C...........   Other local variables
        INTEGER         I, J, S            ! counter and indices
        INTEGER         IOS                ! I/O status

        CHARACTER(ALLLEN3) TSRC        !  tmp source information

        CHARACTER(16) :: PROGNAME = 'PROCINVSRCS' ! program name

C***********************************************************************
C   begin body of subroutine PROCINVSRCS

C.........  Allocate memory for sorted inventory arrays
        ALLOCATE( CIFIP( NSRC ),
     &             CSCC( NSRC ),
     &           CSOURC( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CSOURC', PROGNAME )

        SELECT CASE( CATEGORY )
        CASE( 'AREA' )
            ALLOCATE( XLOCA( NSRC ),
     &                YLOCA( NSRC ),
     &               CELLID( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CELLID', PROGNAME )
        CASE( 'MOBILE' )
            ALLOCATE( IRCLAS( NSRC ),
     &                IVTYPE( NSRC ),
     &                 CLINK( NSRC ),
     &                CVTYPE( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CVTYPE', PROGNAME )
        CASE( 'POINT' )
        END SELECT

C.........  Loop through sources to store sorted arrays
C           for output to I/O API file.
C.........  Keep case statement outside the loops to speed processing
        SELECT CASE ( CATEGORY )
        CASE( 'AREA' )

             DO I = 1, NRAWSRCS

                 S = SRCIDA( I )
                 TSRC = CSOURCA( I )

                 CIFIP( S ) = TSRC( 1:FIPLEN3 )
                 CSCC( S )  = TSRC( SCCPOS3:SCCPOS3+SCCLEN3-1 )
                 CSOURC( S ) = TSRC

            END DO

            XLOCA = BADVAL3   ! array
            YLOCA = BADVAL3   ! array
            CELLID = 0        ! array

        CASE( 'MOBILE' )

            DO I = 1, NRAWSRCS

                S = SRCIDA( I )
                TSRC = CSOURCA( I )

                CIFIP( S )  = TSRC( 1:FIPLEN3 )
                IRCLAS( S ) =
     &              STR2INT( TSRC( RWTPOS3:RWTPOS3+RWTLEN3-1 ) )
                IVTYPE( S ) =
     &              STR2INT( TSRC( VIDPOS3:VIDPOS3+VIDLEN3-1 ) )
                CSCC( S )   = TSRC( MSCPOS3:MSCPOS3+SCCLEN3-1 )
                CLINK( S )  = TSRC( LNKPOS3:LNKPOS3+LNKLEN3-1 )
                CSOURC( S ) = TSRC
                CVTYPE( S ) = 'MOVES'

            END DO

        CASE( 'POINT' )

            DO I = 1, NRAWSRCS

                S = SRCIDA( I )
                TSRC = CSOURCA( I )

                CIFIP( S ) = TSRC( 1:FIPLEN3 )
                CSCC( S ) = TSRC( CH4POS3:CH4POS3+SCCLEN3-1 )

                CSOURC( S ) = TSRC
            END DO

        END SELECT

C.........  Deallocate per-source unsorted arrays
        DEALLOCATE( CSOURCA )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE PROCINVSRCS

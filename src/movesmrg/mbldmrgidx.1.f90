
SUBROUTINE BLDMRGIDX

!***********************************************************************
!  subroutine BLDMRGIDX body starts at line
!
!  DESCRIPTION:
!      The purpose of this subroutine is to allocate and populate indicator
!      arrays that say which pollutants and species are present for each
!      type of input file (inventory, speciation matrix, multiplicative
!      control matrix, reactivity matrix, etc.) for each source category
!      (area, biogenic, mobile, point).  The indicator arrays store the
!      position in the list of i/o api file variables that the given pollutant
!      or species matches.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 2/99 by M. Houyoux
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
    USE MODMERGE, ONLY: NIPPA,&
    &                    NSMATV, NMSPC,&
    &                    TSVDESC,&
    &                    SIINDEX, SPINDEX,&
    &                    EMNAM, EANAM

!.........  This module contains data structures and flags specific to Movesmrg
    USE MODMVSMRG, ONLY: EANAMREP

    IMPLICIT NONE

!...........   INCLUDES:

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   Call allocated arrays
!...........   Group index counter for each source-category-specific list of
!              pollutants and activities.
    INTEGER     KM( NIPPA )    !  mobile

!...........   Allocatable arrays
    INTEGER     EANAMIDX( NIPPA )  ! index from EANAM to TSVDESC

!...........   Other local variables
    INTEGER         J, K, L1, L2, V    !  counters and indices

    INTEGER         IOS      ! i/o error status
    INTEGER         MCNT     ! mobile src var counter
    INTEGER         TGRP     ! tmp group number

    LOGICAL      :: EFLAG = .FALSE.  ! true: error found
    LOGICAL         NEXTGRP  ! true: time to increment group number cntr

    CHARACTER(300)     :: MESG ! message buffer
    CHARACTER(IODLEN3) :: CBUF ! tmp pol-to-species buffer
    CHARACTER(NAMLEN3) :: CPOL ! tmp pollutant buffer
    CHARACTER(NAMLEN3) :: CSPC ! tmp species buffer
    CHARACTER(NAMLEN3) :: PSPC ! tmp previous species
    CHARACTER(NAMLEN3) :: VBUF ! tmp variable name buffer
    CHARACTER(PLSLEN3) :: SVBUF ! tmp speciation name buffer

    CHARACTER(16) :: PROGNAME = 'BLDMRGIDX' ! program name

!***********************************************************************
!   begin body of subroutine BLDMRGIDX

!.........  Allocate memory for building reporting flag array
    ALLOCATE( EANAMREP( NSMATV ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EANAMREP', PROGNAME )
    EANAMIDX = 0
    EANAMREP = .FALSE.

!.........  Allocate memory for variable to pollutant and variable to species
!           indexes.
    ALLOCATE( SIINDEX( NSMATV, 1 ),&
    &          SPINDEX( NSMATV, 1 ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SIINDEX,SPINDEX', PROGNAME )
    SIINDEX = 0  ! array
    SPINDEX = 0  ! array

!.........  Loop through pol-to-species names
    PSPC = ' '
    DO V = 1, NSMATV

!.............  Extract pollutant name and species name
        CBUF = TSVDESC( V )
        L1 = INDEX   ( CBUF, SPJOIN )
        L2 = LEN_TRIM( CBUF )
        CPOL = CBUF(    1:L1-1 )
        CSPC = CBUF( L1+1:L2   )

!.............  Determine location in the list of pollutants/activities
        K = INDEX1( CPOL, NIPPA, EANAM )
        J = INDEX1( CSPC, NMSPC, EMNAM )

!.............  Store pol/act and species indices
        SIINDEX( V, 1 ) = K             ! store pol/act index
        SPINDEX( V, 1 ) = J             ! store species index

!.............  Save index into TSVDESC array from EANAM
        IF( EANAMIDX( K ) .EQ. 0 ) THEN
            EANAMIDX( K ) = V
        END IF
        PSPC  = CSPC

    END DO

!.........  Build array to indicate when report values should be saved
    DO V = 1, NIPPA
        IF ( EANAMIDX( V ) > 0 ) THEN
            EANAMREP( EANAMIDX( V ) ) = .TRUE.
        END IF
    END DO

    RETURN

END SUBROUTINE BLDMRGIDX

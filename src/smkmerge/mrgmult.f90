
SUBROUTINE MRGMULT( NSRC, NG, NL, NMAT1, NMAT2,         &
                    KEY1, KEY2, KEY4, ISPC, FG, FT,     &
                    EMSRC, RINFO, CUMATX, SMATX,        &
                    NX, IX, GMATX, ICNY, GOUT1, GOUT2,  &
                    COUT1, COUT2, COUT4, COUT5 )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine multiplies a source-emissions vector with a gridding
    !      matrix and optionally a speciation array and multiplicative control
    !      array.  Which matrices are applied depend on the setting of the
    !      keys in the subroutine call.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version ??/???? by ???
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related cha
    !************************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !        System
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

    !.....  MODULES for public variables
    !.....  This module contains the major data structure and control flags
    USE MODMERGE, ONLY: NPSRC, ELEVFLAG, ELEVADJ, LFRAC, PINGFLAG,    &
         INLINEFLAG, SRCGRPFLAG, ISRCGRP, EMGGRD, SUBSECFLAG, IGRPNUM, IUGRPIDX

    !.....  This module contains arrays for plume-in-grid and major sources
    USE MODELEV, ONLY: ELEVFLTR

    !.....  This module contains the arrays for state and county summaries
    USE MODSTCY, ONLY: NCOUNTY

    IMPLICIT NONE

    !.....  SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN) :: NSRC        ! number of source
    INTEGER     , INTENT (IN) :: NG          ! (local) number of grid cells
    INTEGER     , INTENT (IN) :: NL          ! (local) number of layers
    INTEGER     , INTENT (IN) :: NMAT1       ! dim 1 for gridding matrix
    INTEGER     , INTENT (IN) :: NMAT2       ! dim 2 for gridding matrix
    INTEGER     , INTENT (IN) :: KEY1        ! inven emissions index
    INTEGER     , INTENT (IN) :: KEY2        ! mult controls index
    INTEGER     , INTENT (IN) :: KEY4        ! speciation matrix index
    INTEGER     , INTENT (IN) :: ISPC        ! output species index
    REAL        , INTENT (IN) :: FG          ! gridded output units conv
    REAL        , INTENT (IN) :: FT          ! st/co tot output units conv
    REAL        , INTENT (IN) :: EMSRC ( NSRC,* )     ! source-based emissions
    REAL        , INTENT (IN) :: RINFO ( NSRC,2 )     ! reactivity information
    REAL        , INTENT (IN) :: CUMATX( NSRC,* )     ! mult control factors
    REAL        , INTENT (IN) :: SMATX ( NSRC,* )     ! speciation factors
    INTEGER     , INTENT (IN) :: NX    ( NG )     ! no. of sources per cell
    INTEGER     , INTENT (IN) :: IX    ( NMAT1 )      ! list of sources per cell
    REAL        , INTENT (IN) :: GMATX ( NMAT2 )      ! gridding coefficients
    INTEGER     , INTENT (IN) :: ICNY  ( NSRC )       ! county index by source
    REAL        , INTENT(OUT) :: GOUT1 ( NG, NL )     ! one-time gridded emis
    REAL     , INTENT(IN OUT) :: GOUT2 ( NG, NL )     ! cumulative gridded emis
    REAL     , INTENT(IN OUT) :: COUT1 ( NCOUNTY, * )    ! no control county
    REAL     , INTENT(IN OUT) :: COUT2 ( NCOUNTY, * )    ! multiplv cntl county
    REAL     , INTENT(IN OUT) :: COUT4 ( NCOUNTY, * )    ! reactivity cntl cnty
    REAL     , INTENT(IN OUT) :: COUT5 ( NCOUNTY, * )    ! all control cntl cnty

    !.....  Other local variables
    INTEGER         C, J, K, L, S       ! counters and indicies
    INTEGER         IDX             ! index to list of counties in grid
    INTEGER         GIDX            ! index to source group
    INTEGER         GRPNUM          ! source group number
    REAL            GFAC            ! tmp gridding factor
    REAL            FG0             ! gridding conv fac div. totals conv fac
    REAL(8)         SUM1            ! sum for GOUT1
    REAL(8)         SUM2            ! sum for GOUT2
    REAL(8)         MULT            ! tmp value with multiplictv controls
    REAL(8)         REAC            ! tmp value with reactivity controls
    REAL(8)         VAL             ! tmp value
    REAL(8)         VMP             ! tmp market penetration value

    CHARACTER(16), PARAMETER :: PROGNAME = 'MRGMULT'     ! program name

    !***********************************************************************
    !   begin body of subroutine MRGMULT

    FG0 = FG / FT

    !.....  If the sources are point sources and there are elevated sources,
    !       transfer ELEVFLTR to ELEVADJ
    IF( NSRC .EQ. NPSRC .AND. (ELEVFLAG .OR. PINGFLAG .OR. INLINEFLAG) ) THEN

        ELEVADJ( 1:NPSRC ) = ELEVFLTR( 1:NPSRC )

    ELSE

        ELEVADJ = 0.

    END IF

    !.....  Check if this is a valid inventory pollutant for this call, and
    !       if the number of layers is one.
    IF( NL .LE. 1 .AND. KEY1 .GT. 0 ) THEN

        !..... If multiplicative controls  speciation
        IF( KEY2 .GT. 0 .AND. KEY4 .GT. 0 ) THEN

            K = 0
            DO C = 1, NG

                SUM1 = GOUT1( C,1 )
                SUM2 = GOUT2( C,1 )

                DO J = 1, NX( C )
                    K = K + 1
                    S = IX( K )
                    IDX  = ICNY( S )
                    GFAC = GMATX( K ) * FT

                    VAL = EMSRC( S,KEY1 )* SMATX( S,KEY4 )* GFAC
                    COUT1( IDX,ISPC ) = COUT1( IDX,ISPC ) + VAL

                    MULT = VAL * CUMATX( S,KEY2 )
                    COUT2( IDX,ISPC ) = COUT2( IDX,ISPC ) + MULT

                    VMP  = RINFO( S,2 )
                    REAC = ( VAL * (1.-VMP) + RINFO( S,1 ) * VMP ) * GFAC
                    COUT4( IDX,ISPC ) = COUT4( IDX,ISPC ) + REAC

                    VAL = MULT * (1.-VMP) + RINFO( S,1 )* VMP* GFAC

                    COUT5( IDX,ISPC ) = COUT5( IDX,ISPC ) + VAL

                    VAL  = VAL * ( 1.-ELEVADJ( S ) )
                    SUM1 = SUM1 + VAL * FG0
                    SUM2 = SUM2 + VAL * FG0

                    IF( SRCGRPFLAG .OR. SUBSECFLAG ) THEN
                        GIDX = 0
                        IF( SRCGRPFLAG ) THEN
                            GIDX = ISRCGRP( S )
                        END IF
                        IF( SUBSECFLAG ) THEN
                            GRPNUM = IGRPNUM( ISRCGRP( S ) )
                            IF( GRPNUM .NE. 0 ) THEN
                                GIDX = IUGRPIDX( GRPNUM )
                            END IF
                        END IF
                        IF( GIDX .NE. 0 ) THEN
                            EMGGRD( C,GIDX ) = EMGGRD( C,GIDX ) + VAL * FG0
                        END IF
                    END IF
                END DO

                GOUT1( C,1 ) = SUM1
                GOUT2( C,1 ) = SUM2

            END DO

        !..... If multiplicative controls only
        ELSE IF( KEY2 .GT. 0 ) THEN

            K = 0
            DO C = 1, NG

                SUM1 = GOUT1( C,1 )
                SUM2 = GOUT2( C,1 )

                DO J = 1, NX( C )
                    K = K + 1
                    S = IX( K )
                    IDX  = ICNY( S )
                    GFAC = GMATX( K ) * FT

                    VAL = EMSRC( S,KEY1 )*GFAC
                    COUT1( IDX,KEY1 ) = COUT1( IDX,KEY1 ) + VAL

                    MULT = VAL * CUMATX( S,KEY2 )
                    COUT2( IDX,KEY1 ) = COUT2( IDX,KEY1 ) + MULT
                    COUT4( IDX,KEY1 ) = COUT4( IDX,KEY1 ) + VAL
                    COUT5( IDX,KEY1 ) = COUT5( IDX,KEY1 ) + MULT

                    VAL  = MULT
                    VAL  = VAL * ( 1.-ELEVADJ( S ) )
                    SUM1 = SUM1 + VAL * FG0
                    SUM2 = SUM2 + VAL * FG0

                    IF( SRCGRPFLAG .OR. SUBSECFLAG ) THEN
                        GIDX = 0
                        IF( SRCGRPFLAG ) THEN
                            GIDX = ISRCGRP( S )
                        END IF
                        IF( SUBSECFLAG ) THEN
                            GRPNUM = IGRPNUM( ISRCGRP( S ) )
                            IF( GRPNUM .NE. 0 ) THEN
                                GIDX = IUGRPIDX( GRPNUM )
                            END IF
                        END IF
                        IF( GIDX .NE. 0 ) THEN
                            EMGGRD( C,GIDX ) = EMGGRD( C,GIDX ) + VAL * FG0
                        END IF
                    END IF
                END DO

                GOUT1( C,1 ) = SUM1
                GOUT2( C,1 ) = SUM2

            END DO

        !.....  If speciation only
        ELSE IF( KEY4 .GT. 0 ) THEN

            K = 0
            DO C = 1, NG

                SUM1 = GOUT1( C,1 )
                SUM2 = GOUT2( C,1 )

                DO J = 1, NX( C )
                    K = K + 1
                    S = IX( K )
                    IDX  = ICNY( S )
                    GFAC = GMATX( K ) * FT

                    VAL = EMSRC( S,KEY1 ) * SMATX( S,KEY4 ) * GFAC
                    COUT1( IDX,ISPC )= COUT1( IDX,ISPC ) + VAL
                    COUT2( IDX,ISPC )= COUT2( IDX,ISPC ) + VAL

                    VMP = RINFO( S,2 )
                    VAL = VAL * (1.-VMP) + RINFO( S,1 ) * VMP * GFAC

                    COUT4( IDX,ISPC )= COUT4( IDX,ISPC ) + VAL
                    COUT5( IDX,ISPC )= COUT5( IDX,ISPC ) + VAL

                    VAL  = VAL * ( 1.-ELEVADJ( S ) )
                    SUM1 = SUM1 + VAL * FG0
                    SUM2 = SUM2 + VAL * FG0

                    IF( SRCGRPFLAG .OR. SUBSECFLAG ) THEN
                        GIDX = 0
                        IF( SRCGRPFLAG ) THEN
                            GIDX = ISRCGRP( S )
                        END IF
                        IF( SUBSECFLAG ) THEN
                            GRPNUM = IGRPNUM( ISRCGRP( S ) )
                            IF( GRPNUM .NE. 0 ) THEN
                                GIDX = IUGRPIDX( GRPNUM )
                            END IF
                        END IF
                        IF( GIDX .NE. 0 ) THEN
                            EMGGRD( C,GIDX ) = EMGGRD( C,GIDX ) + VAL * FG0
                        END IF
                    END IF
                END DO

                GOUT1( C,1 ) = SUM1
                GOUT2( C,1 ) = SUM2

            END DO

        !.....  If inventory pollutant only
        ELSE
            K = 0
            DO C = 1, NG

                SUM1 = GOUT1( C,1 )
                SUM2 = GOUT2( C,1 )

                DO J = 1, NX( C )
                    K = K + 1
                    S = IX( K )
                    IDX = ICNY( S )

                    VAL = EMSRC( S,KEY1 ) * GMATX( K ) * FT
                    COUT1( IDX,KEY1 )= COUT1( IDX,KEY1 ) + VAL
                    COUT2( IDX,KEY1 )= COUT2( IDX,KEY1 ) + VAL
                    COUT4( IDX,KEY1 )= COUT4( IDX,KEY1 ) + VAL
                    COUT5( IDX,KEY1 )= COUT5( IDX,KEY1 ) + VAL

                    VAL  = VAL * ( 1.-ELEVADJ( S ) )
                    SUM1 = SUM1 + VAL * FG0
                    SUM2 = SUM2 + VAL * FG0

                    IF( SRCGRPFLAG .OR. SUBSECFLAG ) THEN
                        GIDX = 0
                        IF( SRCGRPFLAG ) THEN
                            GIDX = ISRCGRP( S )
                        END IF
                        IF( SUBSECFLAG ) THEN
                            GRPNUM = IGRPNUM( ISRCGRP( S ) )
                            IF( GRPNUM .NE. 0 ) THEN
                                GIDX = IUGRPIDX( GRPNUM )
                            END IF
                        END IF
                        IF( GIDX .NE. 0 ) THEN
                            EMGGRD( C,GIDX ) = EMGGRD( C,GIDX ) + VAL * FG0
                        END IF
                    END IF
                END DO

                GOUT1( C,1 ) = SUM1
                GOUT2( C,1 ) = SUM2

            END DO

        END IF      ! End which of controls and speciation

    !.....  If we need to use layer fractions...

    ELSE IF( KEY1 .GT. 0 ) THEN

        !..... If multiplicative controls  speciation  layer fractions
        IF( KEY2 .GT. 0 .AND. KEY4 .GT. 0 ) THEN

            DO L = 1, NL

                K = 0
                DO C = 1, NG

                    SUM1 = GOUT1( C,L )
                    SUM2 = GOUT2( C,L )

                    DO J = 1, NX( C )
                        K = K + 1
                        S = IX( K )
                        IDX  = ICNY( S )
                        GFAC = GMATX( K ) * LFRAC( S,L ) * FT

                        VAL = EMSRC( S,KEY1 )* SMATX( S,KEY4 )* GFAC
                        COUT1( IDX,ISPC ) = COUT1( IDX,ISPC ) + VAL

                        MULT = VAL * CUMATX( S,KEY2 )
                        COUT2( IDX,ISPC ) = COUT2( IDX,ISPC ) + MULT

                        VMP  = RINFO( S,2 )
                        REAC = ( VAL * (1.-VMP) + RINFO(S,1)* VMP )* GFAC
                        COUT4( IDX,ISPC ) = COUT4( IDX,ISPC ) + REAC

                        VAL = MULT* (1.-VMP) + RINFO(S,1)* VMP* GFAC

                        COUT5( IDX,ISPC ) = COUT5( IDX,ISPC ) + VAL

                        VAL  = VAL * ( 1.-ELEVADJ( S ) )
                        SUM1 = SUM1 + VAL * FG0
                        SUM2 = SUM2 + VAL * FG0

                    END DO

                    GOUT1( C,L ) = SUM1
                    GOUT2( C,L ) = SUM2

                END DO
            END DO

        !..... If multiplicative controls and layer fractions
        ELSE IF( KEY2 .GT. 0 ) THEN

            DO L = 1, NL

                K = 0
                DO C = 1, NG

                    SUM1 = GOUT1( C,L )
                    SUM2 = GOUT2( C,L )

                    DO J = 1, NX( C )
                        K = K + 1
                        S = IX( K )
                        IDX  = ICNY( S )
                        GFAC = GMATX( K ) * LFRAC( S,L ) * FT

                        VAL = EMSRC( S,KEY1 )*GFAC
                        COUT1( IDX,KEY1 ) = COUT1( IDX,KEY1 ) + VAL

                        MULT = VAL * CUMATX( S,KEY2 )
                        COUT2( IDX,KEY1 ) = COUT2( IDX,KEY1 ) + MULT
                        COUT4( IDX,KEY1 ) = COUT4( IDX,KEY1 ) + VAL
                        COUT5( IDX,KEY1 ) = COUT5( IDX,KEY1 ) + MULT

                        VAL  = MULT
                        VAL  = VAL * ( 1.-ELEVADJ( S ) )
                        SUM1 = SUM1 + VAL * FG0
                        SUM2 = SUM2 + VAL * FG0

                    END DO

                    GOUT1( C,L ) = SUM1
                    GOUT2( C,L ) = SUM2

                END DO
            END DO

        !..... If speciation and layer fraction
        ELSE IF( KEY4 .GT. 0 ) THEN

            DO L = 1, NL

                K = 0
                DO C = 1, NG

                    SUM1 = GOUT1( C,L )
                    SUM2 = GOUT2( C,L )

                    DO J = 1, NX( C )
                        K = K + 1
                        S = IX( K )
                        IDX  = ICNY( S )
                        GFAC = GMATX( K ) * LFRAC( S,L ) * FT

                        VAL = EMSRC( S,KEY1 )* SMATX( S,KEY4 )* GFAC
                        COUT1( IDX,ISPC ) = COUT1( IDX,ISPC ) + VAL
                        COUT2( IDX,ISPC ) = COUT2( IDX,ISPC ) + VAL

                        VMP  = RINFO( S,2 )
                        VAL  = ( VAL*(1.-VMP) + RINFO( S,1 ) * VMP * GFAC )

                        COUT4( IDX,ISPC ) = COUT4( IDX,ISPC ) + VAL
                        COUT5( IDX,ISPC ) = COUT5( IDX,ISPC ) + VAL

                        VAL  = VAL * ( 1.-ELEVADJ( S ) )
                        SUM1 = SUM1 + VAL * FG0
                        SUM2 = SUM2 + VAL * FG0
                    END DO

                    GOUT1( C,L ) = SUM1
                    GOUT2( C,L ) = SUM2

                END DO
            END DO

        !.....  If inventory pollutant and layer fractions
        ELSE

            DO L = 1, NL

                K = 0
                DO C = 1, NG

                    SUM1 = GOUT1( C,L )
                    SUM2 = GOUT2( C,L )

                    DO J = 1, NX( C )
                        K = K + 1
                        S = IX( K )
                        IDX = ICNY( S )

                        VAL = LFRAC( S,L ) * EMSRC( S,KEY1 ) * GMATX( K ) * FT
                        COUT1( IDX,KEY1 ) = COUT1( IDX,KEY1 ) + VAL
                        COUT2( IDX,KEY1 ) = COUT2( IDX,KEY1 ) + VAL
                        COUT4( IDX,KEY1 ) = COUT4( IDX,KEY1 ) + VAL
                        COUT5( IDX,KEY1 ) = COUT5( IDX,KEY1 ) + VAL

                        VAL  = VAL * ( 1.-ELEVADJ( S ) )
                        SUM1 = SUM1 + VAL * FG0
                        SUM2 = SUM2 + VAL * FG0
                    END DO

                    GOUT1( C,L ) = SUM1
                    GOUT2( C,L ) = SUM2

                END DO
            END DO

        END IF      ! End which of controls and speciation

    END IF      ! End if no inventory emissions or L > 1

    RETURN

END SUBROUTINE MRGMULT

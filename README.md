# M3UTILIO/F90-ized SMOKE

## Notice

**This is a development/technology-transfer repository, not a production
repository.** It contains a GitHub Oct. 28, 2023 version of
SMOKE (`.0.f`), a `M3UTILIO`ized, I/O API 3.2-ized  (`.f`) of
SMOKE, an *findent* (`.1.f90`) free-source-format reference version of
SMOKE constructed from the `.f`, and a cleaned-up `.f90` version of
SMOKE.

**This code version has not been tested** (I don't have the facilities
to do that), but it does compile successfully with Intel *ifort* version
14, GNU *gfortran* version 4.9.2, and Sun *sunf95* version 8.6.  (Due to
the idiosyncratic interpretation of the Fortran Standard in GNU
*gfortran* versions 10 and later, those compiler-versions need testing,
and I don't have a system that has them.)

In several places, copyright, attribution, and license notices that had
been stripped from the code have been replaced.

SMOKE should not have its own private, out-of-date copy of I/O API
routine `NEXTIME`.  This has been removed.

Where appropriate (and easy), codes have been brought up to I/O API 3.2.
For example, *make*-symbol `E132` now gives the compile-flags for the
(non-Standard) "extended-132-column fixed" Fortran source format used by
SMOKE (eliminating the previously-necessary compiler-dependency in
*SMOKE/src/Makeinclude*).

Argument-list bugs that were found, as well as the `FORMAT` bugs in
*src/movesmrg/rdmrclist* and *lib/efsetup* (`MCREFIDX` and `VOLNAM` are
`CHARACTER`, not `INTEGER` as `FORMAT 94010` demands; see 
https://forum.cmascenter.org/t/error-running-movesmerg-for-rpd/4606/4) 
have been fixed.  There is a greatly-simplified  *lib/efsetup.2.f90*.


## Introduction

Fortran [`MODULE
M3UTILIO`](https://cjcoats.github.io/ioapi/M3UTILIO.html) does the
"proper F90 way" of encapsulating the Models-3 I/O API parameter, data
structure, and function-declarations, and has done so since I/O API 3.0
(2002).  In addition to the I/O routines like `WRITE3`, it also provides
`INTERFACE` blocks for almost all the routines in the I/O API, which
provides for compile-time checking of argument lists, etc.:  so far, the
conversion from I/O API-2 style code to `USE M3UTILIO` has found
argument-list errors in every code to which it has been applied.  This
now includes eliminating argument-list bugs in the current (GitHub Oct.
28, 2023) version of SMOKE, upon which this code distribution is based.

Note that `MODULE M3UTILIO` also provides generic/polymorphic interfaces
for a number of routines:  for example,
[`ENVGET`](https://cjcoats.github.io/ioapi/ENVGETS.html) provides a
type-safe name  for calling any of `ENVINT(), ENVINT8(), ENVDBLE(),
ENVREAL(), ENVSTR(), ENVYN(), BENVINT(), BENVINT8(), BENVDBLE(),` and
`BENVREAL()` (These have not been introduced into this version of SMOKE.)

The resulting code has then been converted to Fortran-90 Standard (*.f90*)
free source format, which should have been done easily and quickly as an
automated process using Willem Vermin's [*findent*
program](https://github.com/wvermin/findent), with command-line options:
<pre>
        findent -i4 -k6 -ofree  -RR  !*
</pre>
The resulting reference code was named `.1.f90`.  The final `.f90` codes
needed fix-up on loop-nest and comment indentation, etc., as well as
fixing up numerous botches in the original code, such as 
<pre>
    "( /", "/ )", ". AND.", ". LT."
</pre>
that should *never* have had embedded blanks. The corresponding free
source-format `INCLUDE`-files are named `.h90`.

`Makefile`s were updated to include missing dependencies where the
changes to them had not kept up with the changes to the source-code
base.


## Structure

The directory structure matches that of the default *SMOKE/src* subtree
from GitHub.

There is a reference read-only `.0.f` copy of each of the original GitHub
source files, e.g., *src/biog/czangle.0.f*.

New "fixed-132 codes go by the "standard" SMOKE naming, e.g., 
*src/biog/czangle.f* Intermediate "scratch" or "improved" versions of
the codes have other in-fixes, e.g., *src/biog/tmpbeis4.1.f* and 
*src/biog/tmpbeis4.2.f*. Read-only reference-copy *findent* outputs use
a `.1` in-fix, e.g., *src/biog/czangle.1.f90*, and final "free
.f90-style" codes go by `.f90` (without an in-fix).

For comparing these code-versions, note that programs
[*xxdiff*](https://github.com/blais/xxdiff),
[*tkdiff*](https://sourceforge.net/projects/tkdiff/), and
[*kompare*](https://apps.kde.org/kompare/) are graphical "difference"
programs that display the codes side-by-side with highlighted
differences, making it easy to see what the changes are between all
these file versions, e.g.,
<pre>
xxdiff temporal/wrtsup.0.f temporal/wrtsup.f
</pre>

Note that **the `.f` and `.f90` codes require different *Makefile* and
*Makeinclude* files**.  Files `.f.Makefile` and `.f.Makeinclude` are for
the former, and `.f90.Makefile` and `.f90.Makeinclude` are for the
latter.  To build SMOKE, in the *src* directory copy the relevant pair
of these to the standard *Makefile* and *Makeinclude* names, and then
*make*.  For reference, the original `Makefile` can be found in `Makefile.0`


## Changes

### M3UTILIO

All relevant codes have been `M3UTILIO`-ized.  This entailed removing
`INCLUDE` statements for the I/O API include-files `PARMS3.EXT,
IODECL3.EXT, FDESC3.EXT` as well as external-function declarations for
I/O API fuctions that were called.  Other external function declarations
were put into Fortran-90 style, e.g.
<pre>
        INTEGER, EXTERNAL :: GETFLINE
</pre>
Argument-list bugs were found and fixed for certain calls to `SORTI2`,
`SORTR2`, and `INDEX1`.


### ENVINT, M3ERR, TRIMLEN, etc.

Error-checking was missing for almost all of the `ENVINT, ENVREAL,
ENVYN` calls.  This error-checking has now been added, so if the
user does something inappropriate in a script like assigning a
non-integer where an environment variable should be an integer:

        setenv IFOO dingbats

Then the new code now error-checks to detect this invalid value, report
the problem, and exit (whereas the old code would  have allowed the code
to continue inappropriately with a potentially-bad value).

I/O API routine `M3ERR` was deprecated (replaced by `M3EXIT` and `M3WARN`)
for I/O API Prototype 0.5 of July 1992, and `TRIMLEN` was deprecated
(replaced by Fortran-90 intrinsic `LEN_TRIM`) for I/O API version 3.0
(2002).  Occurrences of both of these have been replaced.


### ALLOCATEs

`ALLOCATE`s were grouped together for improved readability (and improved
performance:  each `ALLOCATE` makes a separate system-call, and these
are expensive).  For example, the more-readable
<pre>
        ALLOCATE(    WSAT( NCOLS, NROWS ),
     &                 RN( NCOLS, NROWS ),
     &                 RC( NCOLS, NROWS ),
     &              TEMPG( NCOLS, NROWS ),
     &             RADYNI( NCOLS, NROWS ),
     &               RATM( NCOLS, NROWS ),
     &                LAI( NCOLS, NROWS ),
     &                 Q2( NCOLS, NROWS ),
     &             RSTOMI( NCOLS, NROWS ),
     &              RSTOM( NCOLS, NROWS ),
     &              USTAR( NCOLS, NROWS ),
     &              RGRND( NCOLS, NROWS ),
     &               PRES( NCOLS, NROWS ),
     &           RAINFALL( NCOLS, NROWS, RHOURS ),
     &              PTYPE( NCOLS, NROWS ),
     &          PULSEDATE( NCOLS, NROWS ),
     &          PULSETIME( NCOLS, NROWS ),
     &              EMPOL( NCOLS, NROWS, NSEF ),
     &              EMISL( NCOLS, NROWS, MSPCS ),
     &              EMISS( NCOLS, NROWS, MSPCS ),
     &              SEMIS( NCOLS, NROWS, NSEF ),
     &            NONAGNO( NCOLS, NROWS ),
     &          NGROWAGNO( NCOLS, NROWS ),
     &           GROWAGNO( NCOLS, NROWS ),
     &               SLAI( NCOLS, NROWS, NLAI ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WSAT...SLAI', PROGNAME )
</pre>
replaces the far more cluttered
<pre>
        ALLOCATE( WSAT( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WSAT', PROGNAME )
        ALLOCATE( RN( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RN', PROGNAME )
        ALLOCATE( RC( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RC', PROGNAME )
        ALLOCATE( TEMPG(NCOLS, NROWS), STAT=IOS )
        CALL CHECKMEM( IOS, 'TEMPG', PROGNAME )
        ALLOCATE( RADYNI(NCOLS, NROWS), STAT=IOS )
        CALL CHECKMEM( IOS, 'RADYNI', PROGNAME )
        ALLOCATE( RATM(NCOLS, NROWS), STAT=IOS )
        CALL CHECKMEM( IOS, 'RATM', PROGNAME )
        ALLOCATE( LAI(NCOLS, NROWS), STAT=IOS )
        CALL CHECKMEM( IOS, 'LAI', PROGNAME )	
        ALLOCATE( Q2(NCOLS, NROWS), STAT=IOS )
        CALL CHECKMEM( IOS, 'Q2', PROGNAME )
        ALLOCATE( RSTOMI(NCOLS, NROWS), STAT=IOS )
        CALL CHECKMEM( IOS, 'RSTOMI', PROGNAME )
        ALLOCATE( RSTOM(NCOLS, NROWS), STAT=IOS )
        CALL CHECKMEM( IOS, 'RSTOM', PROGNAME )
        ALLOCATE( USTAR(NCOLS, NROWS), STAT=IOS )
        CALL CHECKMEM( IOS, 'USTAR', PROGNAME )
        ALLOCATE( RGRND(NCOLS, NROWS), STAT=IOS )
        CALL CHECKMEM( IOS, 'RGRND', PROGNAME )
        ALLOCATE( RAINFALL( NCOLS, NROWS, RHOURS ),STAT=IOS )
        CALL CHECKMEM( IOS, 'RAINFALL', PROGNAME )
        ALLOCATE( PTYPE( NCOLS, NROWS), STAT=IOS )
        CALL CHECKMEM( IOS, 'PTYPE', PROGNAME )
        ALLOCATE( PULSEDATE( NCOLS, NROWS), STAT=IOS )
        CALL CHECKMEM( IOS, 'PULSEDATE', PROGNAME )
        ALLOCATE( PULSETIME( NCOLS, NROWS), STAT=IOS )
        CALL CHECKMEM( IOS, 'PULSETIME', PROGNAME )
        ALLOCATE( EMPOL( NCOLS, NROWS, NSEF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMPOL', PROGNAME )
        ALLOCATE( PRES( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRES', PROGNAME )
        ALLOCATE( EMISL( NCOLS, NROWS, MSPCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMISL', PROGNAME )
        ALLOCATE( EMISS( NCOLS, NROWS, MSPCS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EMISS', PROGNAME )
        ALLOCATE( SEMIS( NCOLS, NROWS, NSEF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SEMIS', PROGNAME )
        ALLOCATE( NONAGNO( NCOLS, NROWS), STAT=IOS )
        CALL CHECKMEM( IOS, 'NONAGNO', PROGNAME )
        ALLOCATE( NGROWAGNO( NCOLS, NROWS), STAT=IOS )
        CALL CHECKMEM( IOS, 'NGROWAGNO', PROGNAME )
        ALLOCATE( GROWAGNO( NCOLS, NROWS), STAT=IOS )
        CALL CHECKMEM( IOS, 'GROWAGNO', PROGNAME )
        ALLOCATE( SLAI( NCOLS, NROWS, NLAI ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SLAI', PROGNAME )
</pre>


### Numerics Etc.

It was recognized early-on in Models-3 development that grid-coordinate
related calculations needed to be done in `DOUBLE PRECISION`, i.e.,
`REAL*8` &mdash; the critical situation being exactly the situation
dealt with in *lib/ingrid*.  This routine has been corrected to do so.

A number of loop nests (especially in biogenics) were found to be in
nest-order that is as cache-hostile as possible. These have been
replaced by the cache-friendly versions.  This should result in improved
performance (where these nests occur), especially for large-grid scenarios. 
See [Optimizing Environmental Models for Microprocessor Based Systems
&mdash; *The Easy
Stuff*](https://cjcoats.github.io/optimization/efficient_models.html).

It was recognized from the very beginning of the Models-3 project in
1991 (and documented as such) that because `REAL` arithmetic is always
subject to machine dependent round-off problems [^1] and therefore
`REAL` values should **never** be tested for exact equality, a robust
scheme was needed for the detection of "missing" for `REAL` values. 
Therefore, (with over-kill) the I/O API provided parameters
<pre>
        REAL, PARAMETER :: BADVAL3 =  -9.999E36  !  for "missing"
        REAL, PARAMETER :: AMISS3  =  -9.000E36  !  for testing
</pre>
with the intention that the code for testing `REAL` "missing" should
look like
<pre>
        REAL    X
        ...
        X = BADVAL3                 !  set X  to "missing"
        ...
        IF ( X .LT. AMISS3 ) THEN   ! test X for "missing"
        ...
</pre>
Even on systems with the most extreme round-off difficulty properties
(e.g., IBM 360 or POWER, Intel x86/x87, early versions of SunOS or
IRIX), this is still a robust test.  In many places in SMOKE (especially
program `UAM2NCF`), the tests were for equality to `AMISS3`, for which
the programmer should consider himself or herself lucky if the answer
comes out "right".  In this new edition of SMOKE, these have all been
chased down and fixed.

The I/O API version of the `FLTERR` and `DBLERR` "these `REAL`s are
certainly unequal, up to reasonable round-off" functions were carefully
tuned with tolerances that handle  the numerical inadequacies of all of
the reasonable compute platforms (including WMO GRIB, which is *quite*
bad).  In several places, the tolerances had been tightened unreasonably
by those un-knowing of the actual problems.  These have been fixed.

In the biogenics, the photolysis process needs the tangent of the zenith
angle.  The code from the original biogenics authors computed the cosine
`CZEN` of the zenith angle, then took the inverse cosine
`ZEN=ARCCOS(CZEN)`, and finally used the tangent `TAN(ZEN)` of the
zenith angle, with the effect of both added computational costs due to
the use of extremely-expensive trigonometric and inverse-trigonometric
functions, and additional round-off error.  The original SMOKE biogenics
took advantage of the high-school trigonometry Pythagorean identities to
compute this as
<pre>
        TAN(ZEN) = SQRT( (1/CZEN)^2 - 1 )
</pre>
with both lower round-off error and much-reduced computational cost.
Someone in the mean-time replaced the improved version with the (lower
quality) original.  This has now been fixed.


### Scratch arrays, `PARAMETER`s, `TRIM`, etc.

Fortran-90 provides very simple and flexible array structure for "auto"
local-variable arrays (which, BTW, provide "leak-proof re-use of
memory).  SMOKE systematically avoids this simplicity, using `ALLOCATE`
and `DEALLOCATE` for what should be "auto" local-variable arrays. Note
that  `DEALLOCATE` does not necessarily reverse the effects of the
matching `ALLOCATE`:  it generally causes "holes" in the program's
memory-map, and if performed repeated constitutes a memory leak. In a
number of places, `ALLOCATE / DEALLOCATE` arrays were turned into auto
arrays.  (Doing this to the maximum extent possible would have been an
extremely tedious task ;-( ) [See `PERMUTI`, below.]

Many constants were moved from `DATA` statements to `PARAMETER` statements.

In many places, constructs involving `CHARACTER`-string lengths of the
form `FOO( 1 : LEN_TRIM( FOO ) )` are replaced by the cheaper,
much-simpler, and  more readable but equivalent `TRIM( FOO )`.  Likeise,
some "stupid" double-assignments like:
<pre>
        ALLOCATE( ENDLEN( MXPTCHR3 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ENDLEN', PROGNAME )
        ENDLEN = 1      ! array
        ENDLEN( 1:MXPTCHR3 ) = PTENDL3( 1:MXPTCHR3 )
</pre>
have been simplified to:
<pre>
        ALLOCATE( ENDLEN( MXPTCHR3 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ENDLEN', PROGNAME )
        ENDLEN( 1:MXPTCHR3 ) = PTENDL3( 1:MXPTCHR3 )
</pre>

Numerous embedded tab-characters were found (these make code-indentation
problematical, depending upon one's editor-settings); these were
replaced by the appropriate numbers of blanks.

In many places, the original code-authors introduced potential bugs due
to a failure to understand how declaration-time initialization for
variables works.  The following sort of declaration
<pre>
        INTEGER :: NFOUND = 0
</pre>
should  almost always instead have a declaration, followed by a separate
initialization-statement at the beginning of the body of the routine:
</pre>
        INTEGER :: NFOUND
        ...
        NFOUND = 0
</pre>
A number of these were fixed.  (It is conceivable that the potential
bugs were never triggered, because the relevant routines were called
only a single time.  However, the routines do not error-check to ensure
that they are not called repeatedly.)


### IOVLEN3 etc.

SMOKE parameters `IOVLEN3`, `IOULEN3`, and `IODLEN3`, which must of necessity
match I/O API parameters `NAMLEN3` and `MXDLEN3` (but were introduced by
a failure to understand how to do this sort of thing properly in Fortran-90),
have been replaced by the matching I/O API parameters, now obtained from
`M3UTILIO`.  This results in the effective elimination of
*src/inc/IOSTRG3.EXT* and *src/inc/IOPRVT3.EXT*.

`PARAMETER ALLFILES` belongs in `MODULE MODFILESETAPI`, not scattered
through the *filesetapi* source code and additionally in `INCLUDE`-file
`SETDECL`.  This has been implemented for the `.f90` version of the
code.


### `'f90` Free source format.

SMOKE was converted to Fortran-90 Standard (*.f90*) free source format,
(from the `M3UTILIO`ized code, using Willem Vermin's [*findent*
program](https://github.com/wvermin/findent), with command-line options:
<pre>
        findent -i4 -k6 -ofree  -RR  !*
</pre>
creating the `.1.f90` version of the code. From this, the final `.f90`
code was generated, fixing up continuation and comment-indentation
styles.

Lines were re-folded to allow for the wider effective source-width
created by the new free-format source indentation, especially when doing
so eliminated "breaks" within phrases and enhanced readability.


### Next Steps

The set of supported map projections should be expanded to include at
least the Mercator and polar map projection types, and ideally to
include the full set of I/O API map projections.  This can easily be
done using routines [`GRID2XY`](https://cjcoats.github.io/ioapi/GRID2XY.html)
and [`XY2XY`](https://cjcoats.github.io/ioapi/XY2XY.html) from I/O API
[`MODULE MODGCTP`](https://cjcoats.github.io/ioapi/MODGCTP.html).

File timestep-consistency checks in the `TMPBEIS*` programs are bogus:
they should be replaced by checks that ask "does this file contain the
data" (formulated in terms, perhaps, of I/O API routine `JSTEP3()`)
instead of "does this file have *exactly* the time steps I want?".

Likewise, hard-coded one-hour-timestep assumptions are inconsistent with
many kinds of modeling and should be removed.  They are so pervasive
that this effort will be quite tedious.

The whole allocation-and-sorting system can and should be simplified
enormously, using a variant of *lib/getfline.f* that also returns the
number of non-comment lines in ASCII files, together with the new I/O
API generic routine [`PERMUTI`](https://cjcoats.github.io/ioapi/PERMUTI.html)
<pre>
        PERMUTI(N, INDX, ARR1 [, ARR2 [, ARR3]] )
</pre>
that sorts its array-arguments `ARR1` etc. in-place, on the basis of
the `INDX` array returned by the `SORTI`, without the need for extra
temporary scratch-arrays.  Note that retrofitting `PERMUTI` would be a
major task.

The "generate new files in a file-set" constructs ensuring that file
sizes do not exceed the (netCDF-2, 1990's) 2GB file size limit have not
been needed since the introduction of netCDF-3 in the late 1990's. 
These constructs are no longer needed and should go away.  Possibly
(since I/O API version 3.2 supports up to 2048 variables and I/O API
3.2-large supports up to 16384 variables), **all** of *src/filesetapi*
should be eliminated, using standard I/O API routines instead. (Note
that netCDF-3 supports file sizes larger than SMOKE's 2 GB assumption,
with a 2GB-per-timestep limit prior to netCDF-3.6 (2002) which does away
with that limit...)

In this matter, *lib/ioapi_grd_size* needs to be changed from type
`INTEGER` to `INTEGER*8`, since with large grids and variable-counts,
timestep size in files may well exceed 2 GB.

Since SMOKE scripts always set environment variable `PROMPTFLAG` to
`NO`, the interactive-style prompt-for-file routines `PROMPTMFILE` and 
`PROMPTFFILE` should not have been used; these should be replaced by the
script-style [*OPEN3*](https://cjcoats.github.io/ioapi/OPEN3.html) and
[*GETEFILE*](https://cjcoats.github.io/ioapi/GETEFILE.html)

I/O API functions
[`STR2INT,STR2REAL,STR2DBLE`](https://cjcoats.github.io/ioapi/STR2S.html)
were originally developed for use in SMOKE (starting with the SMOKE
prototpe version 0.2 (1993), where they are used to read
`INTEGER,REAL,REAL*8` numbers from character strings with error
checking, returning `IMISS3` or `BADVAL3` in case of errors.  As such,
the  new-fangled SMOKE `CHKINT` and `CHKREAL` are superfluous; such
sequences
as 
<pre>
    IF( .NOT. CHKINT( SEGMENT( 2 ) ) ... ) THEN
        ...
    ELSE
        IMON = STR2INT( SEGMENT( 2 ) )
        ...
    END IF
</pre>
are more complicated **and much more expensive** than
<pre>
    IMON = STR2INT( SEGMENT( 2 ) )
    IF( IMON .EQ. IMISS3 ) THEN
        ...
    END IF
</pre>
The original coding should be preferred.



### NOTES  (*btw, GitHub has a Markdown-formatting bug that shows up here*)

[^1]: The Fortran Standard explicitly **refuses** to dictate the quality
of how round-off behaves.  As an extreme example, *no two of the
following* are guaranteed to be equal (a problem found especially on
Cray vector, IBM mainfreame and POWER, Intel x86/x87, and SGI platforms),
although the naive impression is that they should all be equal:
<pre>
        INTEGER, PARAMETER :: A = 1.0 / 3.0
        INTEGER, PARAMETER :: B = 0.333333333333333333  ! with extra-digit overkill
        INTEGER, PARAMETER :: C = FLOAT( 1 ) / FLOAT( 3 )
        INTEGER     D, E, F, G
        ...
        D = 1.0 / 3.0
        E = FLOAT( 1 ) / FLOAT( 3 )
        F = 0.333333333333333333
        READ( '0.333333333333333333', * )  G
        ...
</pre>

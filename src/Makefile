#.........................................................................
# Version "$Id$"
# Copyright (C) 2009 Baron Advanced Meteorological Systems, LLC., and
# (C) 2016-2023 UNC Institute for the Environment
# Distributed under the GNU GENERAL PUBLIC LICENSE version 2
# See file "GPL.txt" for conditions of use.
#.........................................................................
#  Environment Variables:
#       BIN     machine/OS/compiler/mode type. Shows up as suffix
#               for "Makeinclude.${BIN}" to determine compilation
#               flags, and in ${OBJDIR} and $(INSTALL) to determine
#               binary directories
#       INSTALL installation-directory root, used for "make install":
#               SMOKE executables will be installed in $(INSTALL)/${BIN}
#.........................................................................
#  Directories:
#       ${BASEDIR}  is the root directory for the SMOKE source and
#                   the (machine/compiler/flag-specific) binary
#                   object/library/executable directories.
#.........................................................................
#
#       ---------------     Definitions:   -------------------------

include Makeinclude

# find source, object, and module files here:

VPATH = ${OBJDIR}:\
${BASEDIR}/emmod:${BASEDIR}/lib:${BASEDIR}/filesetapi:${BASEDIR}/biog:\
${BASEDIR}/cntlmat:${BASEDIR}/emqa:${BASEDIR}/emutil:${BASEDIR}/grdmat:\
${BASEDIR}/movesmrg:${BASEDIR}/point:${BASEDIR}/smkinven:${BASEDIR}/smkmerge:\
${BASEDIR}/spcmat:${BASEDIR}/temporal

#  emmod:

MODSRC = \
 modar2pt.f90   modbeis3.f90    modbiog.f90         \
 modcntrl.f90   moddayhr.f90    modelev.f90         \
 modemfac.f90   modgrid.f90     modinfo.f90         \
 modlists.f90   modmbset.f90    modmerge.f90        \
 modmet.f90     modmobil.f90    modmvsmrg.f90       \
 modrepbn.f90   modreprt.f90    modsourc.f90        \
 modspro.f90    modstcy.f90     modsurg.f90         \
 modtag.f90     modtmprl.f90    modxref.f90         \
 modgrdlib.f90 

MOD  = $(MODSRC:.f90=.mod)

MODOBJ = $(MODSRC:.f90=.o)

#  lib:

LIBSRC = \
 alocatbl.f90     alocchrt.f90      alocctbl.f90          \
 alocetbl.f90     alocgtbl.f90      alocmtbl.f90          \
 alocptbl.f90     alocstbl.f90      alocttbl.f90          \
 applreac.f90     applumat.f90      bldcsrc.f90           \
 bldenams.f90     blkorcmt.f90      calcrelhum.f90        \
 checkmem.f90     chkcpvar.f90      chkemepi.f90          \
 chkexpscc.f90    chkexpsic.f90     chkgrid.f90           \
 chkint.f90       chkisiz.f90       chkmetem.f90          \
 chkptdef.f90     chkreal.f90       chksrcno.f90          \
 cvtrdtype.f90                                            \
 cvtvehtype.f90   dscm3grd.f90      dscm3lay.f90          \
 dscsprof.f90     efsetup.f90       evalcrit.f90          \
 fillatbl.f90     fillchrt.f90      fillctbl.f90          \
 filletbl.f90     fillgtbl.f90      fillmtbl.f90          \
 fillptbl.f90     fillstbl.f90      fillttbl.f90          \
 find1first.f90   flterr.f90        fltrneg.f90           \
 fltrxref.f90     fmtcsrc.f90       genptcel.f90          \
 genptvcel.f90    genuslst.f90      getbasyr.f90          \
 getcfdsc.f90     getctgry.f90      getdysav.f90          \
 getfline.f90     getformt.f90      getidasz.f90          \
 getifdsc.f90     getiname.f90      getinvyr.f90          \
 getm3epi.f90     getnlist.f90      getntisz.f90          \
 getrfdsc.f90     getsinfo.f90      gettzone.f90          \
 grd2cnty.f90     grdfips.f90       hdrmiss3.f90          \
 ingrid.f90       initem.f90        ioapi_grd_size.f90    \
 rdsprofdsc.f90   multunit.f90                            \
 normtpro.f90     openphys.f90      padnzero.f90          \
 padzero.f90      parscsrc.f90      parsline.f90          \
 pdsetup.f90      polmesg.f90       prclinrc.f90          \
 procspro.f90     progdesc.f90      rd3mask.f90           \
 rdar2pt.f90      rdascc.f90        rdchrscc.f90          \
 rdcodnam.f90     rddates.f90       rdeproc.f90           \
 rdgeocodes.f90   rdgmat.f90        rdgref.f90            \
 rdhdays.f90      rdinvchr.f90      rdinvmap.f90          \
 rdinvpol.f90     rdlines.f90       rdmactdsc.f90         \
 rdmapmask.f90    rdmappol.f90      rdnaicsdsc.f90        \
 rdorsdsc.f90     rdpelv.f90        rdrmat.f90            \
 rdsccdsc.f90     rdsccmap.f90      rdsconv.f90           \
 rdsetmask.f90    rdsicdsc.f90      rdsmat.f90            \
 rdspdprof.f90    rdspdref.f90      rdsprof.f90           \
 rdsrcgrps.f90    rdsref.f90        rdsrg.f90             \
 rdsrgdesc.f90    rdsrghdr.f90      rdstcy.f90            \
 rdtzone.f90      rdumat.f90        rdmxref.f90           \
 rdxclude.f90     readwr3.f90       rmcommnt.f90          \
 setscctype.f90   setsrcdy.f90      srcgrpcnt.f90         \
 tagtable.f90                       unitfac.f90           \
 unitmatch.f90    useexpgeo.f90     verchar.f90           \
 wrchrscc.f90     wrdaymsg.f90      wrorlout.f90          \
 wrsrcgrps.f90    xfractbl.f90      xreftbl.f90

LIBOBJ = $(LIBSRC:.f90=.o)

#  filesetapi:  both ".f90" and ".f90"
#  special case because of both .f90 and .f90 in ${FILSRC}

FILSRC = \
 appendname.f90        chkfileset.f90         \
 chksetdesc.f90        cleanup.f90            \
 closeset.f90          createset.f90          \
 descset.f90           modfileset.f90         \
 openset.f90           promptset.f90          \
 readset.F90           writeset.F90

FILOBJ = \
 appendname.o        chkfileset.o         \
 chksetdesc.o        cleanup.o            \
 closeset.o          createset.o          \
 descset.o           modfileset.o         \
 openset.o           promptset.o          \
 readset.o           writeset.o

FIL  = $(FILSRC:.f90=.mod)

#  biog:

BIOSRC = \
 czangle.f90     getpar.f90      getparb.f90                                \
 hrbeis3.f90     hrbeis_361.f90  hrbeis4.f90      hrno.f90  hrno_beis4.f90  \
 normbeis3.f90   normbeis314.f90 normbeis361.f90  normbeis4.f90             \
 prebmet.f90     rdb3fac.f90     rdb4fac.f90      rdnormbeis4_efs.f90       \
 rdb4fac_csv.f90 rdbpro.f90      tmpbeis3.f90                               \
 tmpbeis314.f90  tmpbeis361.f90  tmpbeis4.f90     normbeis370.f90
 
BIOOBJ =  $(BIOSRC:.f90=.o)

#  cntlmat:

CTLSRC = \
 alocpkts.f90     asgncntl.f90      cntlmat.f90       \
 errpkts.f90      fillcdat.f90      fillcntl.f90      \
 genmultc.f90     genproj.f90       genreact.f90      \
 opencmat.f90     openctmp.f90      openpmat.f90      \
 openrmat.f90     pktloop.f90       procpkts.f90      \
 rdpacket.f90     wcntlrep.f90      wrctmp.f90        \
 wrrmat.f90
 
CTLOBJ =  $(CTLSRC:.f90=.o)

#  emqa (smkreport):

RPTSRC = \
 asgnbins.f90    bldrepidx.f90    genrprt.f90  \
 openrepin.f90   openrepout.f90   qarepin.f90  \
 rdgrps.f90      rdrepin.f90      rdrprts.f90  \
 repmrggrd.f90   repunits.f90     scanrepc.f90 \
 selectsrc.f90   smkreport.f90    wrrephdr.f90 \
 wrrepout.f90    rdssup.f90
 
RPTOBJ =  $(RPTSRC:.f90=.o)

#  emutil:

UTISRC = \
 aggwndw.f90   beld3to2.f90       saregroup.f90          \
 cemscan.f90   extractida.f90     gcntl4carb.f90         \
 gentpro.f90   geofac.f90         invsplit.f90           \
 layalloc.f90  metcombine.f90     metscan.f90            \
 pktreduc.f90  smk2emis.f90       surgtool.f90           \
 uam2ncf.f90   inlineto2d.f90     
 
UTIOBJ =  $(UTISRC:.f90=.o)

#  grdmat:

GRDSRC = \
 asgnsurg.f90  genagmat.f90       genggmat.f90    \
 genlgmat.f90  genmgmat.f90       genpgmat.f90    \
 grdmat.f90    opengmat.f90       rdsrg4grd.f90   \
 setfrac.f90   sizgmat.f90
 
GRDOBJ =  $(GRDSRC:.f90=.o)

#  movesmrg:

MOVSRC = \
 bldprocidx.f90     bldsrccell.f90      hourmet.f90     \
 mbldmrgidx.f90     met4moves.f90       mgetmrgev.f90   \
 minitstcy.f90      mmrgonams.f90       mmrgunits.f90   \
 mmrgvnams.f90      mopenmrgin.f90      mopenmrgout.f90 \
 movesmrg.f90       mwrmrgrep.f90       rdcfpro.f90     \
 rdfmref.f90        rdmrclist.f90       rdspdist.f90    \
 rdrpdemfacs.f90    rdrphemfacs.f90     rdrppemfacs.f90 \
 rdrpvemfacs.f90    rdspdpro.f90        setrefcnty.f90  \
 wrtemprof.f90      
 
MOVOBJ =  $(MOVSRC:.f90=.o)

#  point:

PNTSRC = \
 asgngrps.f90       elevpoint.f90       fire_plmris.f90     \
 fire_postplm.f90   fire_preplm.f90     laypoint.f90        \
 mxgrpemis.f90      openeout.f90        openlayout.f90      \
 plmris.f90         plsprd.f90          plumesti.f90        \
 postplm.f90        preplm.f90          rpelvcfg.f90        \
 wpingstk.f90
 
PNTOBJ =  $(PNTSRC:.f90=.o)

#  smkinven:

INVSRC = \
 adjustinv.f90    asgnar2pt.f90    asgnnhapx.f90    chklstfl.f90     \
 fixstk.f90       formlist.f90     genmedsout.f90   \
 genpdout.f90     gethdr.f90       getpdinfo.f90    \
 grwinven.f90     initinfo.f90     opengrwout.f90   \
 openinvin.f90    openinvout.f90   openpdout.f90    \
 procar2pt.f90    procinven.f90    procinvsrcs.f90  \
 rdcempd.f90      rdcemsum.f90     rddataff10ar.f90 \
 rddataff10mb.f90 rddataff10pt.f90 rddatamedspt.f90 \
 rddatantiar.f90  rddatantifr.f90  rddatantimb.f90  \
 rddatantinp.f90  rddatantipt.f90  rdemspd.f90      \
 rdff10pd.f90     rdgrdapi.f90     rdgrdncf.f90     \
 rdinvdata.f90    rdinvsrcs.f90    rdlooppd.f90     \
 rdmedsinfo.f90   rdmedspd.f90     rdorlfr.f90      \
 rdsrcff10ar.f90  rdsrcff10mb.f90  rdsrcff10pt.f90  \
 rdsrcmedspt.f90  rdsrcntiar.f90   rdsrcntifr.f90   \
 rdsrcntimb.f90   rdsrcntinp.f90   rdsrcntipt.f90   \
 setnonhap.f90    smkinven.f90     srcmem.f90       \
 wrepinven.f90    wrinvchr.f90     wrinvemis.f90    \
 wrinvpol.f90     wrpdemis.f90     wrptref.f90

INVOBJ =  $(INVSRC:.f90=.o)

#  smkmerge:

MRGSRC = \
 allocmrg.f90    bldmrgidx.f90    getmrgev.f90     \
 initstcy.f90    mrgasciielev.f90 mrgbio.f90       \
 mrgelev.f90     mrggrid.f90      mrgmult.f90      \
 mrgonams.f90    mrgunits.f90     mrgvnams.f90     \
 openmrgin.f90   openmrgout.f90   rdmrginv.f90     \
 setoutdate.f90  smkmerge.f90     wmrgelev.f90     \
 wmrgemis.f90    wrmrgrep.f90     mrgpt.f90

MRGOBJ =  $(MRGSRC:.f90=.o)

#  spcmat:

SPCSRC = \
 asgnspro.f90       asgntag.f90        chknonhap.f90      \
 opensmat.f90       rdcombo.f90        rdtag.f90          \
 spcmat.f90

SPCOBJ =  $(SPCSRC:.f90=.o)

#  temporal:

TMPSRC = \
 asgntpro.f90    genhemis.f90     mktmat.f90       \
 opentmp.f90     opentmpin.f90    proctpro.f90     \
 setdaylt.f90    temporal.f90     tmnamunt.f90     \
 updtmat.f90     wrtsup.f90

TMPOBJ =  $(TMPSRC:.f90=.o)

OBJ  =  $(MODOBJ) $(FILOBJ) $(LIBOBJ) $(BIOOBJ) $(CTLOBJ) $(RPTOBJ) $(UTIOBJ) \
        $(GRDOBJ) $(MOVOBJ) $(PNTOBJ) $(INVOBJ) $(MRGOBJ) $(SPCOBJ) $(TMPOBJ)

MLIB  = libemmod.a

SLIB  = libsmoke.a

FLIB  = libfileset.a

UTIL  = aggwndw     beld3to2    cemscan     extractida  \
        geofac      invsplit    metcombine  metscan     pktreduc    \
        smk2emis    surgtool    uam2ncf     gentpro     gcntl4carb  \
        layalloc    inlineto2d  saregroup

EXE   = normbeis3   normbeis4   tmpbeis3    tmpbeis4    cntlmat     \
        smkreport   grdmat      \
        movesmrg    met4moves   laypoint    elevpoint   smkinven    \
        grwinven    smkmerge    mrggrid     mrgelev     spcmat      \
        temporal    mrgpt                                           \
        ${UTIL}


######################################################################
#                           High-level targets:

all:    lib exe

lib:    ${MLIB} ${SLIB} ${FLIB}

exe:    ${EXE}

dir:    ${OBJDIR}

emmod:  ${MOD} ${MODOBJ} ${MLIB}

smklib: ${SLIB}

biog:   normbeis3 normbeis4 tmpbeis3 tmpbeis4

cntl:   cntlmat

report: smkreport

grd  :  grdmat

moves:  movesmrg    met4moves

point:  laypoint    elevpoint

inven:  smkinven    grwinven

merge:  smkmerge    mrggrid     mrgelev      mrgpt

spc  :  spcmat

tmprl:  temporal

util :  ${UTIL}

clean:
	cd ${OBJDIR}; rm -f ${OBJ} ${EXE} ${MOD} ${FIL} ${MLIB} ${SLIB} ${FLIB}

rmexe:
	cd ${OBJDIR}; rm ${EXE}

relink:
	make BIN=${BIN} -i rmexe ; make BIN=${BIN} all

install:  ${INSTDIR} ${EXE}
	echo "Installing SMOKE in ${INSTDIR}"
	cd ${OBJDIR};  cp ${EXE} ${INSTDIR}


#      -----------------------   RULES:   -------------------------

%.o : %.mod        #  Disable "gmake"s obnoxious implicit Modula-2 rule !!
%.o : %.f
%.o : %.F

.f90.o .f90.mod:
	cd ${OBJDIR};  ${FC} ${FFLAGS} -c $<

.F90.o .F90.mod:
	cd ${OBJDIR};  ${FC} ${FFLAGS} -c $<


#  ---------------------------  Dependencies:  --------------------------
# libraries must be up-to-date to build executables

${EXE} :  ${MLIB} ${SLIB} ${FLIB}       # libraries must be up-to-date to build executables

# modules must be up-to-date to build ".o"-files

alocatbl.o    : $(MODXREF)
alocchrt.o    : $(MODXREF)
alocctbl.o    : $(MODXREF)
alocetbl.o    : $(MODXREF)
alocgtbl.o    : $(MODXREF)
alocmtbl.o    : $(MODXREF)
alocptbl.o    : $(MODXREF)
alocstbl.o    : $(MODXREF)
alocttbl.o    : $(MODXREF)
applumat.o    : $(MODGRID)
bldenams.o    : $(MODLISTS)
chkgrid.o     : $(MODGRID)
chkmetem.o    : $(MODGRDLIB)
chkptdef.o    : $(MODINFO)
dscsprof.o    : $(MODINFO)
dscsprof.o    : $(MODSPRO)
efsetup.o     : $(MODEMFAC)
fillatbl.o    : $(MODXREF)
fillchrt.o    : $(MODINFO)
fillctbl.o    : $(MODXREF)
filletbl.o    : $(MODXREF)
fillgtbl.o    : $(MODXREF)
fillmtbl.o    : $(MODXREF)
fillptbl.o    : $(MODXREF)
fillstbl.o    : $(MODXREF)
fillttbl.o    : $(MODXREF)
fltrxref.o    : $(MODSOURC) $(MODLISTS) $(MODINFO)
fmtcsrc.o     : $(MODINFO)
genptvcel.o   : $(MODGRID)  $(MODGRDLIB)
genuslst.o    : $(MODSOURC) $(MODLISTS) $(MODINFO)
getbasyr.o    : $(MODSOURC)
getctgry.o    : $(MODINFO)
getntisz.o    : $(MODLISTS)
getsinfo.o    : $(MODFILESET) $(MODINFO)
gettzone.o    : $(MODSTCY)
grd2cnty.o    : $(MODSTCY) $(MODSURG) $(MODGRID)
grdfips.o     : $(MODSURG) $(MODGRID) $(MODMET)
normtpro.o    : $(MODTMPRL)
openphys.o    : $(MODFILESET)
procar2pt.o   : $(MODSOURC) $(MODXREF) $(MODINFO) $(MODAR2PT)
prclinrc.o    : $(MODREPRT) $(MODINFO)
procspro.o    : $(MODSPRO)
rdar2pt.o     : $(MODXREF) $(MODLISTS) $(MODAR2PT)
rdcodnam.o    : $(MODLISTS)
rddates.o     : $(MODTMPRL)
rdeproc.o     : $(MODEMFAC) $(MODINFO)
rdgeocodes.o  : $(MODSTCY)
rdgref.o      : $(MODXREF) $(MODINFO)
rdhdays.o     : $(MODTMPRL)
rdinvchr.o    : $(MODSOURC)
rdinvmap.o    : $(MODINFO) $(MODFILESET)
rdinvpol.o    : $(MODFILESET)
rdmactdsc.o   : $(MODLISTS)
rdmappol.o    : $(MODINFO) $(MODFILESET)
rdnaicsdsc.o  : $(MODLISTS)
rdorsdsc.o    : $(MODSTCY)
rdpelv.o      : $(MODELEV)
rdsccdsc.o    : $(MODLISTS)
rdsccmap.o    : $(MODMOBIL)
rdsconv.o     : $(MODLISTS) $(MODSPRO) $(MODINFO) $(MODAR2PT)
rdsetmask.o   : $(MODFILESET)
rdsicdsc.o    : $(MODLISTS)
rdsprofdsc.o  : $(MODSOURC) $(MODSOURC)
rdsmat.o      : $(MODFILESET)
rdspdprof.o   : $(MODMOBIL)
rdspdref.o    : $(MODXREF)
rdsprof.o     : $(MODSPRO)  $(MODLISTS)
rdsrcgrps.o   : $(MODMERGE) $(MODLISTS) $(MODGRID) $(MODELEV) $(MODSURG) $(MODSOURC)
rdsref.o      : $(MODXREF)  $(MODINFO)
rdsrgdesc.o   : $(MODSURG)
rdsrg.o       : $(MODSURG) $(MODGRID)
rdstcy.o      : $(MODSTCY)
rdtref.o      : $(MODXREF) $(MODINFO) $(MODTMPRL)
rdxclude.o    : $(MODXREF) $(MODINFO)
setscctype.o  : $(MODINFO)
setsrcdy.o    : $(MODINFO)
srcgrpcnt.o   : $(MODMERGE) $(MODGRID) $(MODELEV)
tagtable.o    : $(MODXREF)  $(MODSPRO) $(MODTAG) $(MODINFO)
wrchrscc.o    : $(MODLISTS) $(MODINFO)
wrorlout.o    : $(MODSOURC) $(MODSTCY) $(MODLISTS) $(MODINFO)
wrsrcgrps.o   : $(MODMERGE) $(MODGRID) $(MODELEV) $(MODGRDLIB) $(MODFILESET)
xfractbl.o    : $(MODSPRO)
xreftbl.o     : $(MODXREF) $(MODINFO)
wrrmat.o      : $(MODFILESET)

modgrdlib.o   : $(MODGRID)

chkfileset.o  : $(MODFILESET)
chksetdesc.o  : $(MODFILESET)
cleanup.o     : $(MODFILESET)
closeset.o    : $(MODFILESET)
createset.o   : $(MODFILESET)
descset.o     : $(MODFILESET)
openset.o     : $(MODFILESET)
promptset.o   : $(MODFILESET)
readset.o     : $(MODFILESET)
writeset.o    : $(MODFILESET)

normbeis314.o : $(MODBEIS3) $(MODGRDLIB)
normbeis361.o : $(MODBEIS3) $(MODGRDLIB)
normbeis370.o : $(MODBEIS3) $(MODGRDLIB)
normbeis4.o   : $(MODBEIS3) $(MODGRDLIB)

tmpbeis314.o  : $(MODSPRO) $(MODBEIS3)
tmpbeis361.o  : $(MODSPRO) $(MODBEIS3)
tmpbeis4.o    : $(MODSPRO) $(MODBEIS3)

alocpkts.o    : $(MODCNTRL)
asgncntl.o    : $(MODSOURC) $(MODXREF)
asgncntl.o    : $(MODINFO)  $(MODCNTRL)
cntlmat.o     : $(MODINFO)  $(MODFILESET)
fillcdat.o    : $(MODCNTRL) $(MODXREF)  $(MODINFO)
genmultc.o    : $(MODSOURC) $(MODCNTRL) $(MODINFO) $(MODFILESET)
genproj.o     : $(MODSOURC) $(MODCNTRL) $(MODINFO) $(MODFILESET)
genreact.o    : $(MODSOURC) $(MODCNTRL) $(MODINFO) $(MODSPRO)
opencmat.o    : $(MODCNTRL) $(MODINFO)  $(MODFILESET)
opencmat.o    : $(MODFILESET)
openpmat.o    : $(MODINFO)  $(MODCNTRL) $(MODFILESET)
openrmat.o    : $(MODINFO)  $(MODFILESET)
pktloop.o     : $(MODXREF)  $(MODINFO)
procpkts.o    : $(MODSOURC) $(MODXREF) $(MODCNTRL) $(MODINFO)
rdpacket.o    : $(MODINFO)
wcntlrep.o    : $(MODSOURC) $(MODCNTRL) $(MODINFO)
wrctmp.o      : $(MODINFO)

asgnbins.o    : $(MODSOURC) $(MODLISTS) $(MODREPRT) $(MODREPBN) $(MODGRID)
asgnbins.o    : $(MODELEV)  $(MODSTCY)  $(MODINFO)  $(MODREPRT) $(MODREPBN) $(MODTMPRL)
bldrepidx.o   : $(MODCNTRL) $(MODINFO)  $(MODFILESET)
genrprt.o     : $(MODSOURC) $(MODREPRT) $(MODREPBN) $(MODTMPRL) $(MODCNTRL) $(MODINFO) $(MODFILESET)
openrepin.o   : $(MODREPRT) $(MODTMPRL) $(MODREPBN) $(MODCNTRL) $(MODINFO)  $(MODFILESET)
openrepout.o  : $(MODREPRT)
qarepin.o     : $(MODREPRT) $(MODSTCY)  $(MODINFO)
rdgrps.o      : $(MODREPRT) $(MODGRID)  $(MODSTCY)
rdrepin.o     : $(MODSOURC) $(MODREPRT) $(MODREPBN) $(MODCNTRL) $(MODELEV) $(MODFILESET)
rdrepin.o     : $(MODLISTS) $(MODGRID)  $(MODINFO)  $(MODFILESET)
rdrprts.o     : $(MODREPRT) $(MODINFO)
rdssup.o      : $(MODREPRT) $(MODINFO)  $(MODSOURC)
repmrggrd.o   : $(MODSOURC) $(MODREPRT) $(MODREPBN) $(MODGRID)
repunits.o    : $(MODREPRT) $(MODREPBN) $(MODTMPRL) $(MODINFO)
scanrepc.o    : $(MODREPRT) $(MODINFO)
selectsrc.o   : $(MODSOURC) $(MODREPRT) $(MODREPBN) $(MODELEV)  $(MODINFO) $(MODSPRO)
smkreport.o   : $(MODREPRT) $(MODREPBN) $(MODGRID)  $(MODINFO)
wrrephdr.o    : $(MODSOURC) $(MODLISTS) $(MODREPRT) $(MODREPBN) $(MODSTCY) $(MODGRID) $(MODINFO)
wrrepout.o    : $(MODSOURC) $(MODLISTS) $(MODREPRT) $(MODREPBN) $(MODSTCY) $(MODINFO)

aggwndw.o     : $(MODGRID)
gentpro.o     : $(MODINFO) $(MODSOURC) $(MODGRID) $(MODSURG) $(MODSTCY)
geofac.o      : $(MODFILESET)
inlinto2d.o   : $(MODGRID) $(MODGRDLIB)
invsplit.o    : $(MODINFO)
smk2emis.o    : $(MODINFO)
surgtool.o    : $(MODSURG) $(MODGRID)

asgnsurg.o    : $(MODSOURC) $(MODXREF) $(MODINFO) $(MODSURG)
genagmat.o    : $(MODSOURC) $(MODGRID) $(MODINFO) $(MODSURG) $(MODLISTS) $(MODXREF) $(MODGRDLIB)
genggmat.o    : $(MODSOURC) $(MODGRID) $(MODINFO) $(MODGRDLIB)
genlgmat.o    : $(MODSOURC) $(MODGRID) $(MODINFO) $(MODGCTP) $(M3UTILIO)
genmgmat.o    : $(MODSOURC) $(MODGRID) $(MODINFO) $(MODSURG) $(MODLISTS) $(MODXREF)
grdmat.o      : $(MODSOURC) $(MODGRID) $(MODINFO) $(MODXREF) $(MODSURG)  $(MODFILESET) $(MODGRDLIB)
opengmat.o    : $(MODINFO)
rdsrg4grd.o   : $(MODSURG)  $(MODGRID)
setfrac.o     : $(MODSURG)  $(MODXREF)
sizgmat.o     : $(MODSURG)  $(MODGRID) $(MODGRDLIB)

mbldmrgidx.o  : $(MODMERGE) $(MODMVSMRG)
bldprocidx.o  : $(MODMERGE) $(MODMVSMRG)
bldsrccell.o  : $(MODMVSMRG)
mgetmrgev.o   : $(MODMERGE) $(MODINFO) $(MODMVSMRG)
hourmet.o     : $(MODMBSET) $(MODMET)  $(MODSURG)
minitstcy.o   : $(MODMERGE) $(MODSTCY) $(MODSOURC) $(MODMVSMRG) $(MODLISTS)
met4moves.o   : $(MODMBSET) $(MODGRID) $(MODSURG)  $(MODSTCY)   $(MODMET)
movesmrg.o    : $(MODMERGE) $(MODGRID) $(MODSTCY)  $(MODLISTS)  $(MODMBSET) $(MODSOURC) $(MODMVSMRG) $(MODFILESET)
mmrgonams.o   : $(MODMERGE)
mmrgunits.o   : $(MODMERGE) $(MODMVSMRG)
mmrgvnams.o   : $(MODMERGE) $(MODLISTS) $(MODMVSMRG)
mopenmrgin.o  : $(MODMERGE) $(MODINFO)  $(MODSOURC)  $(MODEMFAC) $(MODGRID) $(MODMVSMRG) $(MODFILESET)
mopenmrgout.o : $(MODMERGE) $(MODGRID)  $(MODMVSMRG) $(MODFILESET)
rdcfpro.o     : $(MODMERGE) $(MODLISTS) $(MODMBSET)  $(MODMVSMRG)
rdfmref.o     : $(MODMBSET)
rdmrclist.o   : $(MODMBSET) $(MODMVSMRG)
rdmxref.o     : $(MODMBSET)
rdrpdemfacs.o : $(MODMBSET)  $(MODMVSMRG) $(MODMERGE) $(MODLISTS)
rdrphemfacs.o : $(MODMBSET)  $(MODMVSMRG) $(MODMERGE) $(MODLISTS)
rdrppemfacs.o : $(MODMBSET)  $(MODMVSMRG) $(MODMERGE) $(MODLISTS)
rdrpvemfacs.o : $(MODMBSET)  $(MODMVSMRG) $(MODMERGE) $(MODLISTS)
rdspdpro.o    : $(MODMVSMRG) $(MODLISTS)  $(MODMOBIL)
rdspdist.o    : $(MODMVSMRG) $(MODLISTS)  $(MODMOBIL)
setrefcnty.o  : $(MODMVSMRG) $(MODINFO)   $(MODSOURC) $(MODLISTS) $(MODMBSET)
mwrmrgrep.o   : $(MODMERGE)  $(MODSTCY)   $(MODLISTS) $(MODMVSMRG)
wrmrgrep.o    : $(MODGRID)
wrtemprof.o   : $(MODMET)  $(MODMBSET)

asgngrps.o    : $(MODSOURC) $(MODELEV) $(MODINFO) $(MODGRDLIB)
elevpoint.o   : $(MODSOURC) $(MODELEV) $(MODINFO) $(MODGRID) $(MODGRDLIB)
laypoint.o    : $(MODSOURC) $(MODELEV) $(MODINFO) $(MODGRID) $(MODGRDLIB) $(MODFILESET)
mxgrpemis.o   : $(MODSOURC) $(MODELEV) $(MODINFO) $(MODGRID) $(MODGRDLIB) $(MODFILESET)
openeout.o    : $(MODELEV)  $(MODINFO)
openlayout.o  : $(MODELEV)  $(MODINFO) $(MODGRID)
rpelvcfg.o    : $(MODREPRT) $(MODELEV) $(MODINFO)
wpingstk.o    : $(MODELEV)

asgnar2pt.o   : $(MODSOURC) $(MODELEV) $(MODINFO)
asgnnhapx.o   : $(MODSOURC) $(MODELEV) $(MODINFO)
fixstk.o      : $(MODSOURC)
formlist.o    : $(MODLISTS) $(MODINFO)
genmedsout.o  : $(MODINFO)  $(MODDAYHR)
genpdout.o    : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR) $(MODSTCY)
gethdr.o      : $(MODLISTS) $(MODSTCY)  $(MODINFO)
getpdinfo.o   : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR)
grwinven.o    : $(MODSOURC) $(MODINFO)  $(MODFILESET)
initinfo.o    : $(MODINFO)
opengrwout.o  : $(MODINFO)  $(MODFILESET)
openinvin.o   : $(MODINFO)  $(MODLISTS) $(MODXREF) $(MODDAYHR) $(MODMOBIL)
openinvout.o  : $(MODINFO)  $(MODFILESET)
openpdout.o   : $(MODINFO)
procar2pt.o   : $(MODSOURC) $(MODXREF)  $(MODINFO) $(MODAR2PT)
procinven.o   : $(MODSOURC) $(MODLISTS) $(MODINFO)
procinvsrcs.o : $(MODSOURC) $(MODINFO)
rdcempd.o     : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR)
rdcemsum.o    : $(MODDAYHR)
rddataff10ar.o: $(MODINFO) $(MODDAYHR)
rddataff10mb.o: $(MODINFO) $(MODLISTS)
rddataff10pt.o: $(MODINFO) $(MODDAYHR)
rddatamedspt.o: $(MODINFO) $(MODSOURC)
rddatantiar.o : $(MODINFO)
rddatantifr.o : $(MODINFO)
rddatantimb.o : $(MODINFO)
rddatantinp.o : $(MODINFO)
rddatantipt.o : $(MODINFO)
rdemspd.o     : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR)
rdff10pd.o    : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR) $(MODSTCY)
rdgrdapi.o    : $(MODSOURC) $(MODINFO)
rdgrdncf.o    : $(MODSOURC) $(MODINFO)
rdinvdata.o   : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODMOBIL)
rdinvsrcs.o   : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODMOBIL)
rdlooppd.o    : $(MODINFO)  $(MODLISTS)
rdmedsinfo.o  : $(MODSOURC)
rdmedspd.o    : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR) $(MODSTCY)
rdorlfr.o     : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR)
rdsrcff10ar.o : $(MODLISTS) $(MODDAYHR)
rdsrcff10mb.o : $(MODLISTS)
rdsrcff10pt.o : $(MODLISTS) $(MODDAYHR)
rdsrcmedspt.o : $(MODSOURC)
rdsrcntiar.o  : $(MODLISTS)
rdsrcntifr.o  : $(MODINFO) $(MODLISTS)
rdsrcntimb.o  : $(MODLISTS)
rdsrcntinp.o  : $(MODLISTS)
rdsrcntipt.o  : $(MODLISTS)
setnonhap.o   : $(MODSOURC) $(MODLISTS) $(MODXREF) $(MODINFO)
smkinven.o    : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR)
srcmem.o      : $(MODSOURC)
wrepinven.o   : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODAR2PT)
wrinvchr.o    : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODFILESET)
wrinvemis.o   : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODFILESET)
wrinvpol.o    : $(MODFILESET)
wrpdemis.o    : $(MODSOURC) $(MODLISTS) $(MODINFO) $(MODDAYHR)
wrptref.o     : $(MODINFO)

allocmrg.o    : $(MODMERGE) $(MODELEV) $(MODCNTRL) $(MODSTCY) $(MODGRID)
bldmrgidx.o   : $(MODMERGE)
getmrgev.o    : $(MODMERGE) $(MODINFO)
initstcy.o    : $(MODMERGE) $(MODSTCY)
mrgelev.o     : $(MODMERGE) $(MODCNTRL) $(MODELEV)
mrggrid.o     : $(MODGRID)  $(MODFILESET)
mrgmult.o     : $(MODMERGE) $(MODELEV)  $(MODSTCY)
mrgonams.o    : $(MODMERGE)
mrgunits.o    : $(MODMERGE)
mrgvnams.o    : $(MODMERGE) $(MODLISTS) 
openmrgin.o   : $(MODMERGE) $(MODINFO)  $(MODELEV) $(MODGRID) $(MODFILESET)
openmrgout.o  : $(MODMERGE) $(MODELEV)  $(MODGRID) $(MODFILESET)
rdmrginv.o    : $(MODMERGE) $(MODLISTS) $(MODSOURC)
smkmerge.o    : $(MODMERGE) $(MODCNTRL) $(MODELEV) $(MODLISTS) $(MODSTCY) $(MODSURG) $(MODGRID) $(MODFILESET)
wmrgelev.o    : $(MODMERGE) $(MODSOURC) $(MODELEV) $(MODGRID)
wmrgemis.o    : $(MODMERGE) $(MODELEV)  $(MODFILESET)
wrmrgrep.o    : $(MODMERGE) $(MODSTCY)  $(MODGRID)

asgnspro.o    : $(MODLISTS) $(MODXREF)  $(MODSPRO) $(MODINFO)
asgntag.o     : $(MODSOURC) $(MODTAG)   $(MODINFO)
chknonhap.o   : $(MODLISTS) $(MODSPRO)
opensmat.o    : $(MODINFO)  $(MODTAG)   $(MODFILESET)
rdcombo.o     : $(MODSPRO)  $(MODLISTS)
rdtag.o       : $(MODXREF)  $(MODSPRO)  $(MODINFO) $(MODLISTS) $(MODTAG)
spcmat.o      : $(MODSPRO)  $(MODEMFAC) $(MODINFO) $(MODLISTS) $(MODTAG) $(MODFILESET)

asgntpro.o    : $(MODSOURC) $(MODTMPRL) $(MODINFO)
genhemis.o    : $(MODSOURC) $(MODXREF)  $(MODINFO) $(MODTMPRL) $(MODDAYHR)
mktmat.o      : $(MODSOURC) $(MODTMPRL) $(MODDAYHR)
opentmp.o     : $(MODINFO)  $(MODTMPRL) $(MODFILESET)
opentmpin.o   : $(MODINFO)
proctpro.o    : $(MODSOURC) $(MODTMPRL) $(MODINFO)
setdaylt.o    : $(MODSOURC) $(MODSTCY)  $(MODINFO)
temporal.o    : $(MODSOURC) $(MODTMPRL) $(MODINFO) $(MODLISTS) $(MODEMFAC) $(MODDAYHR) $(MODMBSET) $(MODFILESET)
tmnamunt.o    : $(MODINFO)
updtmat.o     : $(MODSOURC) $(MODTMPRL) $(MODDAYHR)
wrtsup.o      : $(MODSOURC) $(MODTMPRL) $(MODINFO)


######################################################################
#                            builds:

${OBJDIR}:
	echo "Making " ${OBJDIR} ; mkdir -p ${OBJDIR}

##############################  Libraries

${MLIB}:  ${MODOBJ}
	cd ${OBJDIR}; $(AR) ${ARFLAGS} ${MLIB} ${MODOBJ}

${SLIB}: ${LIBOBJ}
	cd ${OBJDIR}; $(AR) ${ARFLAGS} ${SLIB} ${LIBOBJ}

${FLIB}: ${FILOBJ}
	cd ${OBJDIR}; echo ${FILOBJ}; $(AR) ${ARFLAGS} ${FLIB} ${FILOBJ}

##############################  SMOKE programs

tmpbeis3:  czangle.o getpar.o getparb.o hrbeis_361.o hrbeis3.o prebmet.o rdbpro.o \
 tmpbeis3.o tmpbeis314.o tmpbeis361.o hrno.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

tmpbeis4: czangle.o getparb.o hrbeis4.o  prebmet.o rdbpro.o tmpbeis4.o hrno_beis4.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

normbeis3: normbeis3.o rdb3fac.o rdb4fac.o rdb4fac_csv.o normbeis314.o normbeis361.o normbeis370.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

normbeis4: normbeis4.o rdnormbeis4_efs.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

cntlmat: $(CTLOBJ)
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

smkreport: $(RPTOBJ)
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

grdmat: $(GRDOBJ) 
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

movesmrg: movesmrg.o mbldmrgidx.o bldprocidx.o bldsrccell.o \
 mgetmrgev.o minitstcy.o mmrgonams.o mmrgunits.o  mmrgvnams.o   \
 mopenmrgin.o mopenmrgout.o rdfmref.o  rdmrclist.o rdmxref.o    \
 rdrpdemfacs.o rdrppemfacs.o rdrpvemfacs.o rdrphemfacs.o rdspdpro.o     \
 rdcfpro.o setrefcnty.o mwrmrgrep.o rdspdist.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) $^ ${LIBS} -o $@

met4moves: met4moves.o rdmxref.o rdfmref.o hourmet.o wrtemprof.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) $^ ${LIBS} -o $@

laypoint: laypoint.o openlayout.o preplm.o plmris.o postplm.o plsprd.o \
          fire_preplm.o fire_plmris.o fire_postplm.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

elevpoint: elevpoint.o asgngrps.o openeout.o plumesti.o wpingstk.o  \
           mxgrpemis.o rpelvcfg.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

smkinven: adjustinv.o asgnar2pt.o asgnnhapx.o chklstfl.o fixstk.o formlist.o \
 genmedsout.o genpdout.o gethdr.o getpdinfo.o initinfo.o      \
 openinvin.o openinvout.o openpdout.o procar2pt.o procinven.o \
 procinvsrcs.o rdcempd.o rdcemsum.o rddataff10ar.o rddataff10mb.o rddataff10pt.o \
 rddatamedspt.o rddatantiar.o rddatantifr.o rddatantimb.o rddatantinp.o rddatantipt.o \
 rdemspd.o rdff10pd.o rdgrdapi.o rdgrdncf.o rdinvdata.o rdinvsrcs.o \
 rdlooppd.o rdmedsinfo.o rdmedspd.o rdorlfr.o rdsrcff10ar.o rdsrcff10mb.o \
 rdsrcff10pt.o rdsrcmedspt.o rdsrcntiar.o rdsrcntifr.o rdsrcntimb.o rdsrcntinp.o \
 rdsrcntipt.o setnonhap.o smkinven.o srcmem.o wrepinven.o wrinvchr.o \
 wrinvemis.o wrinvpol.o wrpdemis.o wrptref.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

grwinven: opengrwout.o srcmem.o wrinvchr.o wrinvpol.o grwinven.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

smkmerge: smkmerge.o allocmrg.o bldmrgidx.o getmrgev.o mrgmult.o mrgonams.o \
 mrgvnams.o openmrgin.o openmrgout.o wmrgemis.o wrmrgrep.o initstcy.o \
 rdmrginv.o mrgelev.o mrgunits.o mrgbio.o wmrgelev.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

mrggrid: mrggrid.o setoutdate.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

mrgelev: mrgasciielev.o setoutdate.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

mrgpt: mrgpt.o setoutdate.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)
 
spcmat: $(SPCOBJ)
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

temporal: $(TMPOBJ)
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

##############################  emutil programs:

aggwndw: aggwndw.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

beld3to2: beld3to2.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

cemscan: cemscan.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

extractida: extractida.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

gcntl4carb: gcntl4carb.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

gentpro: gentpro.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

geofac: geofac.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

inlineto2d: inlineto2d.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

invsplit: invsplit.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

layalloc: layalloc.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

metcombine: metcombine.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

metscan: metscan.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

pktreduc: pktreduc.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

saregroup: saregroup.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

smk2emis: smk2emis.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

surgtool: surgtool.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)

uam2ncf: uam2ncf.o
	cd ${OBJDIR}; $(FC) $(FFLAGS) -o $(@) $^ $(LIBS)



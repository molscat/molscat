# Copyright (C) 2018 J. M. Hutson & C. R. Le Sueur
# Distributed under the GNU General Public License, version 3
#
# This makefile written by CR Le Sueur 2015-2018
#
# set up directories =========================================================
OBJDIR  = .
EXECDIR = .
##############################################################################
#
# set flags for linking ======================================================
#
# (we are using gfortran for these tests)
#
# possible alternatives:
# 1) gfortran
Compiler = gfortran
LIBS     = -llapack -lblas
#
# 2) ifort
#Compiler = ifort
#LIBS     =
#
# 3) pgf90
#Compiler = pgf90
#LIBS     =

##############################################################################
#
# set up recipe rules ========================================================
#
# set flags for compilation
COMPILE.f = $(Compiler) -c  -O0 -g
LINK.f    = $(Compiler)
#
# pattern matching rules for generating object files from source code ========
#
$(OBJDIR)/%.o : %.f
	$(COMPILE.f) $< -o $@
$(OBJDIR)/%.o: %.f90
	$(COMPILE.f) $< -o $@
#
# specific recipe for object files that have other dependencies ==============
#
$(OBJDIR)/pot_data_Tiemann-%.o: pot_data_Tiemann-%.f pot_data_tiemann.mod
	$(COMPILE.f) $< -o $@
##############################################################################
#
# recipe rules for module files (to ensure they are remade when necessary) ===
#
physical_constants.mod: physical_constants.f $(OBJDIR)/physical_constants.o
	@true
sizes.mod:              sizes_module.f       $(OBJDIR)/sizes_module.o
	@true
angles.mod:             angles_module.f      $(OBJDIR)/angles_module.o
	@true
efvs.mod:               efvs_module.f        $(OBJDIR)/efvs_module.o
	@true
potential.mod:          potential_module.f   $(OBJDIR)/potential_module.o
	@true
basis_data.mod:         basis_module.f       $(OBJDIR)/basis_module.o
	@true
base9_suite.mod:        base9-skeleton.f     $(OBJDIR)/base9-skeleton.o
	@true
##############################################################################
#
# lists of object codes that have a dependency on a particular module =======
#
# physical constants module =================================================
#
PHYS_CONST_DEPS = base9-1S_3Sigma_cpld.o base9-alk_alk.o \
                  bd.driver.o chckru.o eavg.o ecnv.o extpot-Rg_HX.o \
                  fld.driver.o mol.driver.o pot_data_Tiemann-Rb2-2010.o \
                  surbas.o

(addprefix $(OBJDIR)/,$(PHYS_CONST_DEPS)): physical_constants.mod
#
# sizes module ==============================================================
#
SIZES_DEPS =      angles_module.o bd.driver.o basis_module.o eavg.o \
                  fld.driver.o iosbin.o iosout.o j6to4.o locpol.o \
                  mol.driver.o prbr.o restrt.o setefv.o potential_module.o

(addprefix $(OBJDIR)/,$(SIZES_DEPS)): sizes.mod
#
# angles module =============================================================
#
ANGLES_DEPS =     iosbin.o iosout.o potenl.o vrtp-He_H2O.o vrtp-Rg_HX-eta.o \
                  vrtp-H2_H2.o

(addprefix $(OBJDIR)/,$(ANGLES_DEPS)): angles.mod
#
# efvs module ===============================================================
#
EFVS_DEPS =       base.o base9-1S_3Sigma_cpld.o base9-alk_alk.o \
                  bd.driver.o efv9-dummy.o fld.driver.o header.o ldrdwt.o \
                  locpol.o mol.driver.o output.o prnsum.o rdsigu.o \
                  scalam-default.o setefv.o sread.o wvhead.o cheint.o

(addprefix $(OBJDIR)/,$(EFVS_DEPS)): efvs.mod
#
# potential module ==========================================================
#
POTL_DEPS =       base.o base9-1S_3Sigma_cpld.o base9-alk_alk.o bd.driver.o \
                  calca.o chckru.o cheint.o cntrct.o dermat.o ev_ord.o \
                  fld.driver.o getlev.o mol.driver.o odprop.o output.o \
                  potenl.o restrt.o rmprop.o scctrl.o scwave.o setefv.o \
                  thresh.o vlchk.o vstar-Tiemann.o wavmat.o wavvec.o ytrans.o

(addprefix $(OBJDIR)/,$(POTL_DEPS)): potential.mod
#
# basis_data module =========================================================
#
BASIS_DEPS =      asrot.o base.o base9-1S_3Sigma_cpld.o base9-alk_alk.o \
                  getlev.o iosbin.o iosout.o iospb.o \
                  j6to4.o mol.driver.o output.o prbr.o restrt.o set4.o set6.o \
                  set6c.o setbas.o surbas.o thresh.o \

(addprefix $(OBJDIR)/,$(BASIS_DEPS)): basis_data.mod
##############################################################################
#
# Now set aliases to shorten lists of dependencies ===========================
#
MODULES =   physical_constants.o efvs_module.o sizes_module.o potential_module.o \
            basis_module.o angles_module.o
CLOCK =     gclock.o gdate.o
DIAG =      diagvc.o diagvl.o
LA_UTILS =  dgemul.o syminl-alloc.o dsyneg.o syminv.o
#
# set up dependencies ======================================================
#
TRANS =     ytrans.o multop.o ev_ord.o chnsrt.o wvcalc.o
LDPROPS =   ldprop.o mdprop.o mgprop.o maprop.o trnsfm.o trnsp.o \
            aiprop.o potent.o corr.o maxmgv.o outmat.o spropn.o waveig.o \
            airymp.o scairy.o yinit.o dermat.o
BASE =      mcgcpl.o cpl21.o cpl22.o cpl23.o cpl4.o qsymtp.o rsymtp.o cpl25.o \
            gsymtp.o set6.o ipasym.o esymtp.o fsymtp.o couple.o cpl3.o \
            surbas.o find.o order.o cplout.o idpart.o setbas.o set4.o check6.o \
            set6c.o asrot.o j6to4.o
UTILITIES = dmsym.o dsyfil.o matprn.o chkstr.o
LIBUTILS  =
ANG_MOM =   j3j000.o j6j.o j9j.o threej.o sixj.o xninej.o parsgn.o thrj.o
WMAT =      wavmat.o wavvec.o pertrb.o scapot-default.o
ALLDRV =    base.o brent.o calck.o chckru.o cheint.o drcalc.o drset.o ecnv.o \
            ecnvx.o efprop.o findrm.o idpchk.o progms.o propst.o setefv.o \
            stsrch.o thresh.o thrlst.o vlchk.o wvhead.o \
            $(TRANS) $(LDPROPS) $(BASE) $(UTILITIES) $(LIBUTILS) $(ANG_MOM) $(WMAT)
IOSBIN =    chck6i.o plm.o set6i.o get102.o gaussp.o gasleg.o zbes.o \
            yrr.o sig6.o iosclc.o iosout.o iospb.o isutp.o
PRBR =      prbr.o dasize.o prbr3.o
RESTRT =    sread.o epsum.o restrt.o
PROPS =     dvprop.o stabil.o rmprop.o sgnchk.o vvprop.o \
            pert1.o pert2.o delrd.o odprop.o wkb.o
YTOS =      ytok.o rbesjy.o rbessk.o ytokg.o ksym.o ktos.o
#
# set up dependencies needed for each program ===============================
#
FLD = bdctrl.o bdpset.o bd.iosbin.o bd.prpini.o evmtch.o prnsum.o $(ALLDRV)
BD =  cntrct.o $(FLD)
MOL = calca.o convrg.o eavg.o findrx.o getlev.o header.o iosbin.o iosdrv.o \
      ldrdwt.o locpol.o mol.prpini.o nexte.o output.o rdsigu.o psiasy.o \
      scctrl.o scpset.o scwave.o \
      $(ALLDRV) $(IOSBIN) $(PRBR) $(RESTRT) $(PROPS) $(YTOS)
#
# set up sets of dependencies needed for POTENL ============================
#
POTUTILS    = gaussp.o gasleg.o zbes.o gaushp.o herm.o hrecur.o yrr.o plm.o
GP-POTENL   = potenl.o       $(POTUTILS)
RGH2-POTENL = potenl-Rg_H2.o $(POTUTILS)
#
# dependencies needed for system-specific potentials ========================
#
BASIC-POTENL   = $(GP-POTENL) vstar-dummy.o       vrtp-dummy.o
H2H2-POTENL    = $(GP-POTENL) vstar-H2_H2.o       vrtp-H2_H2.o
ARHX-POTENL    = $(GP-POTENL) vstar-dummy.o       vrtp-Rg_HX-eta.o   extpot-Rg_HX.o \
                              sumleg.o pleg.o pasleg.o assleg.o
ARCH4-POTENL   = $(GP-POTENL) vstar-Rg_CH4_Buck.o vrtp-extpot_2ang.o extpot-Ar_CH4.o
ARCO2-POTENL   = $(GP-POTENL) vstar-dummy.o       vrtp-extpot_1ang.o extpot-Ar_CO2.o \
                              sumleg.o
MGNH-POTENL    = $(GP-POTENL) vstar-Mg_NH.o       vrtp-dummy.o
#
# note that modules must come before object codes that use them
RB2-POTENL     = pot_data_Tiemann-Rb2-2010.o \
                 $(GP-POTENL) vstar-Tiemann.o     vrtp-dummy.o
#
# dependencies needed for plug-in basis-set suites ==========================
#
UNUSED-BASE9 = base9-skeleton.o \
               efv9-dummy.o degen9-nondegenerate.o thrsh9-dummy.o
ALK2-BASE9   = base9-alk_alk.o \
               efv9-dummy.o degen9-nondegenerate.o
MGNH-BASE9   = base9-1S_3Sigma_cpld.o \
               efv9-dummy.o degen9-nondegenerate.o potin9-example.o
# definitions for MOLSCAT post-processors
#
# set up 'normal' set of object files dependencies ==========================
#
BASIC_MOL = $(MODULES) main1M.o mol.driver.o $(CLOCK) $(DIAG) $(LA_UTILS) $(MOL)
BASIC_BND = $(MODULES) main1M.o  bd.driver.o $(CLOCK) $(DIAG) $(LA_UTILS) $(BD)
BASIC_FLD = $(MODULES) main1M.o fld.driver.o $(CLOCK) $(DIAG) $(LA_UTILS) $(FLD)
##############################################################################
#
# list of targets (which are all executables) ================================
#
PROGS =     molscat-basic \
            molscat-H2_H2 molscat-Rg_HX  molscat-Rg_CO2 \
            molscat-Rg_H2 molscat-Ar_CH4 molscat-Mg_NH \
            molscat-Rb2 \
            bound-basic \
            bound-Rg_HX   bound-Rb2      bound-Ar_CH4 \
            bound-H2_H2   bound-Rg_H2    bound-Mg_NH \
            bound-Rg_CO2  \
            field-basic \
            field-Mg_NH   field-Rb2
##############################################################################
#
# rules for making executables ===============================================
#
.PHONY: $(PROGS) # $(TESTPROGS)
#$(PROGS) : %: $(EXECDIR)/%

$(EXECDIR)/molscat-basic     : $(addprefix $(OBJDIR)/,$(BASIC_MOL) \
                                                      $(UNUSED-BASE9) \
                                                      $(BASIC-POTENL))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/molscat-H2_H2     : $(addprefix $(OBJDIR)/,$(BASIC_MOL) \
                                                      $(UNUSED-BASE9) \
                                                      $(H2H2-POTENL))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/molscat-Rg_HX     : $(addprefix $(OBJDIR)/,$(BASIC_MOL) \
                                                      $(UNUSED-BASE9) \
                                                      $(ARHX-POTENL))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/molscat-Rg_CO2    : $(addprefix $(OBJDIR)/,$(BASIC_MOL) \
                                                      $(UNUSED-BASE9) \
                                                      $(ARCO2-POTENL))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/molscat-Rg_H2     : $(addprefix $(OBJDIR)/,$(BASIC_MOL) \
                                                      $(UNUSED-BASE9) \
                                                      $(RGH2-POTENL))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/molscat-Ar_CH4    : $(addprefix $(OBJDIR)/,$(BASIC_MOL) \
                                                      $(UNUSED-BASE9) \
                                                      $(ARCH4-POTENL))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/molscat-Mg_NH     : $(addprefix $(OBJDIR)/,$(BASIC_MOL) \
                                                      $(MGNH-BASE9) \
                                                      $(MGNH-POTENL))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/molscat-Rb2       : $(addprefix $(OBJDIR)/,$(BASIC_MOL) \
                                                      $(ALK2-BASE9) \
                                                      $(RB2-POTENL))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/bound-basic       : $(addprefix $(OBJDIR)/,$(BASIC_BND) \
                                                      $(UNUSED-BASE9) \
                                                      $(BASIC-POTENL))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/bound-Rg_HX       : $(addprefix $(OBJDIR)/,$(BASIC_BND) \
                                                      $(UNUSED-BASE9) \
                                                      $(ARHX-POTENL))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/bound-Rb2         : $(addprefix $(OBJDIR)/,$(BASIC_BND)  \
                                                      $(ALK2-BASE9) \
                                                      $(RB2-POTENL))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/bound-H2_H2       : $(addprefix $(OBJDIR)/,$(BASIC_BND) \
                                                      $(UNUSED-BASE9) \
                                                      $(H2H2-POTENL))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/bound-Rg_H2       : $(addprefix $(OBJDIR)/,$(BASIC_BND) \
                                                      $(UNUSED-BASE9) \
                                                      $(RGH2-POTENL))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/bound-Ar_CH4      : $(addprefix $(OBJDIR)/,$(BASIC_BND) \
                                                      $(UNUSED-BASE9) \
                                                      $(ARCH4-POTENL))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/bound-Rg_CO2      : $(addprefix $(OBJDIR)/,$(BASIC_BND) \
                                                      $(UNUSED-BASE9) \
                                                      $(ARCO2-POTENL))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/bound-Mg_NH       : $(addprefix $(OBJDIR)/,$(BASIC_BND) \
                                                      $(MGNH-BASE9) \
                                                      $(MGNH-POTENL))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/field-basic       : $(addprefix $(OBJDIR)/,$(BASIC_FLD) \
                                                      $(UNUSED-BASE9) \
                                                      $(BASIC-POTENL))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/field-Rb2         : $(addprefix $(OBJDIR)/,$(BASIC_FLD) \
                                                      $(ALK2-BASE9) \
                                                      $(RB2-POTENL))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/field-Mg_NH       : $(addprefix $(OBJDIR)/,$(BASIC_FLD) \
                                                      $(MGNH-BASE9) \
                                                      $(MGNH-POTENL))
	$(LINK.f) $^ $(LIBS) -o $@

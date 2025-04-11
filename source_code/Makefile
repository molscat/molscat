# Copyright (C) 2020 J. M. Hutson & C. R. Le Sueur
# Distributed under the GNU General Public License, version 3
#
# This makefile written by CR Le Sueur 2015-2020
#
# set up directories ==========================================================
OBJDIR  = .
EXECDIR = .
###############################################################################
#
# set flags for linking =======================================================
#
# (we are using gfortran for these tests)
#
# possible alternatives:
# 1) gfortran
Compiler = gfortran
LIBS     = -llapack -lblas
# 2) ifort
#Compiler = ifort
#LIBS     =
#
# 3) pgf90
#Compiler = pgf90
#LIBS     =

###############################################################################
#
# set up recipe rules =========================================================
#
# set flags for compilation
COMPILE.f = $(Compiler) -std=legacy -c  -O0 -g
LINK.f    = $(Compiler)
#
# pattern matching rules for generating object files from source code =========
#
$(OBJDIR)/%.o : %.f
	$(COMPILE.f) $< -o $@
$(OBJDIR)/%.o: %.f90
	$(COMPILE.f) $< -o $@
#
# specific recipe for object files that have other dependencies ===============
#
$(OBJDIR)/pot_data_Tiemann-%.o: pot_data_Tiemann-%.f pot_data_tiemann.mod
	rm $(OBJDIR)/vstar-Tiemann*.o ; $(COMPILE.f) $< -o $@
###############################################################################
#
# recipe rules for module files (to ensure they are remade when necessary) ====
#
###############################################################################
#
# lists of the object files that depend on a particular module ================
#
# angles module ===============================================================
#
ANGLES_DEPS = iosbin.o iosout.o potenl.o vrtp-extpot_1ang.o \
              vrtp-extpot_2ang.o vrtp-H2_H2.o vrtp-Rg_HX-eta.o  
(addprefix \$(OBJDIR)/,$(ANGLES_DEPS)): angles.mod
#
angles.mod:         angles_module.f
	@true
#
# basis_data module ===========================================================
#
BASIS_DATA_DEPS = asrot.o base9-1S_3Sigma_cpld.o base9-alk_alk_ucpld.o base.o \
                  bd.driver.o fld.driver.o getlev.o iosbin.o iosout.o \
                  iospb.o j6to4.o mol.driver.o prbr.o restrt.o \
                  set4.o set6c.o set6.o setbas.o sproc.o \
                  surbas.o thresh.o   
(addprefix \$(OBJDIR)/,$(BASIS_DATA_DEPS)): basis_data.mod
#
basis_data.mod:         basis_data_module.f
	@true
#
# efvs module =================================================================
#
EFVS_DEPS = base9-1S_3Sigma_cpld.o base9-alk_alk_ucpld.o base.o \
            bckgrd-dummy.o bd.driver.o efv9-dummy.o fld.driver.o header.o \
            ldrdwt.o locpol.o mol.driver.o prnsum.o rdsigu.o restrt.o \
            scapot-default.o setefv.o sproc.o wvhead.o 
(addprefix \$(OBJDIR)/,$(EFVS_DEPS)): efvs.mod
#
efvs.mod:         efvs_module.f
	@true
#
# physical_constants module ===================================================
#
PHYSICAL_CONSTANTS_DEPS = base9-1S_3Sigma_cpld.o base9-alk_alk_ucpld.o \
                          bd.driver.o chckru.o eavg.o ecnv.o extpot-Ar_CH4.o \
                          extpot-Ar_CO2.o extpot-Rg_HX.o fld.driver.o \
                          mol.driver.o pot_data_Tiemann-Rb2-2010.o surbas.o \
                          vstar-Tiemann.o
(addprefix \$(OBJDIR)/,$(PHYSICAL_CONSTANTS_DEPS)): physical_constants.mod
#
physical_constants.mod:         physical_constants_module.f
	@true
#
# potential module ============================================================
#
POTENTIAL_DEPS = base9-1S_3Sigma_cpld.o base9-alk_alk_ucpld.o \
                 base9-skeleton.o base.o bd.driver.o chckru.o cheint.o \
                 cntrct.o dermat.o fld.driver.o getlev.o gnham.o hammat.o \
                 mol.driver.o multop.o odprop.o potenl.o potenl-Rg_H2.o \
                 restrt.o scwave.o setefv.o sproc.o sumlam.o thresh.o \
                 vstar-Tiemann.o ytrans.o   
(addprefix \$(OBJDIR)/,$(POTENTIAL_DEPS)): potential.mod
#
potential.mod:         potential_module.f
	@true
#
# sizes module ================================================================
#
SIZES_DEPS = angles_module.o basis_data_module.o bd.driver.o eavg.o \
             fld.driver.o iosbin.o iosout.o j6to4.o locpol.o \
             mol.driver.o potential_module.o prbr.o restrt.o setefv.o
(addprefix \$(OBJDIR)/,$(SIZES_DEPS)): sizes.mod
#
sizes.mod:         sizes_module.f
	@true
###############################################################################
#
# Set aliases to shorten dependency lists =====================================
###############################################################################
#
# Set aliases to shorten dependency lists =====================================
#
MODULES =   physical_constants_module.o sizes_module.o efvs_module.o \
            potential_module.o basis_data_module.o angles_module.o
#MODULES =   physical_constants_2014.o sizes_module.o efvs_module.o \
#            potential_module.o basis_data_module.o angles_module.o
CLOCK =     gclock.o gdate.o
DIAG =      diagvc.o diagvl.o
LA_UTILS =  dgemul.o syminl-alloc.o dsyneg.o syminv.o
#
# set up dependencies =========================================================
#
YTRANS      = chnsrt.o ev_ord.o gtwvec.o multop.o ytrans.o
PROPS_LD    = aiprop.o airymp.o corr.o ldprop.o maprop.o maxmgv.o mdprop.o \
              mgprop.o outmat.o potent.o scairy.o spropn.o trnsfm.o trnsp.o \
              waveig.o yinit.o
PROPS_OTHER = delrd.o dvprop.o pert1.o pert2.o odprop.o rmprop.o sgnchk.o \
              stabil.o vvprop.o wkb.o
BASE        = asrot.o base.o check6.o couple.o cpl21.o cpl22.o cpl23.o \
              cpl25.o cpl3.o cpl4.o cplout.o dmsym.o esymtp.o find.o fsymtp.o \
              gsymtp.o idpart.o ipasym.o j6to4.o mcgcpl.o order.o qsymtp.o \
              rsymtp.o set4.o set6.o set6c.o setbas.o surbas.o
UTILITIES   = chkstr.o dsyfil.o matprn.o
LIBUTILS    =
ANG_MOM     = j3j000.o j6j.o j9j.o parsgn.o sixj.o threej.o thrj.o xninej.o
W_MATRIX    = dermat.o hammat.o sumlam.o
CORE_ALL    = brent.o calck.o chckru.o cheint.o drcalc.o drset.o ecnv.o \
              ecnvx.o findrm.o gnham.o idpchk.o prprop.o progms.o propst.o \
              setefv.o stsrch.o thresh.o thrlst.o wvbksb.o wvhead.o \
              $(YTRANS) $(PROPS_LD) $(BASE) $(UTILITIES) $(LIBUTILS) \
              $(ANG_MOM) $(W_MATRIX)
IOS         = chck6i.o dasize.o get102.o gaussp.o gasleg.o iosbin.o iosclc.o \
              iosdrv.o iosout.o iospb.o isutp.o plm.o set6i.o sig6.o yrr.o \
              zbes.o
PRBR        = prbr.o prbr3.o
RESTRT      = restrt.o sread.o
YTOS        = ksym.o ktos.o rbesjy.o rbessk.o ytok.o ytokg.o
#
# set up dependencies needed for each program =================================
#
FLD_SUBS = bdctrl.o bdpset.o bd.prpini.o evmtch.o iosbin-dummy.o prnsum.o \
           $(CORE_ALL)
BD_SUBS  = cntrct.o $(FLD_SUBS)
MOL_SUBS = calca.o convrg.o eavg.o epsum.o erange.o findrx.o getlev.o \
           header.o ldrdwt.o locpol.o mol.prpini.o nexte.o \
           psiasy.o rdsigu.o scctrl.o scpset.o scwave.o sproc.o \
           $(CORE_ALL) $(IOS) $(PRBR) $(RESTRT) $(PROPS_OTHER) $(YTOS)
REPLACEABLE = bckgrd-dummy.o pertrb.o reiniv-dummy.o scapot-default.o
#
# dependency lists containing potential routines  =============================
#
POTUTILS    = gasleg.o gaushp.o gaussp.o herm.o hrecur.o plm.o yrr.o zbes.o
POTENL-GP   = potenl.o $(POTUTILS)
#
# set up 'normal' set of object files dependencies ============================
#
CORE_MOL = $(MODULES) mol.driver.o $(CLOCK) $(DIAG) $(LA_UTILS) $(MOL_SUBS)
CORE_BND = $(MODULES)  bd.driver.o $(CLOCK) $(DIAG) $(LA_UTILS) $(BD_SUBS)
CORE_FLD = $(MODULES) fld.driver.o $(CLOCK) $(DIAG) $(LA_UTILS) $(FLD_SUBS)
#
# dependency lists for system-specific potentials =============================
#
POT-BASIC   = $(POTENL-GP) vstar-dummy.o       vrtp-dummy.o
POT-H2H2    = $(POTENL-GP) vstar-H2_H2.o       vrtp-H2_H2.o
POT-ARHX    = $(POTENL-GP) vstar-dummy.o       vrtp-Rg_HX-eta.o \
                              extpot-Rg_HX.o sumleg.o pleg.o pasleg.o \
                              assleg.o damp.o
POT-ARCH4   = $(POTENL-GP) vstar-Rg_CH4_Buck.o vrtp-extpot_2ang.o \
                              extpot-Ar_CH4.o
POT-ARCO2   = $(POTENL-GP) vstar-dummy.o       vrtp-extpot_1ang.o \
                              extpot-Ar_CO2.o sumleg.o
POT-MGNH    = $(POTENL-GP) vstar-Mg_NH.o       vrtp-dummy.o
POT-RGH2    = potenl-Rg_H2.o $(POTUTILS)
#
# note that modules must be earlier in dependency lists than object codes
# that use them because they must be compiled first
POT-RB2     = pot_data_Tiemann-Rb2-2010.o \
                 $(POTENL-GP) vstar-Tiemann.o     vrtp-dummy.o
#
# dependency lists for plug-in basis-set suites ===============================
#
BASE9-UNUSED = base9-skeleton.o \
               efv9-dummy.o degen9-nondegenerate.o thrsh9-dummy.o
BASE9-ALK2   = base9-alk_alk_ucpld.o \
               efv9-dummy.o degen9-nondegenerate.o
BASE9-MGNH   = base9-1S_3Sigma_cpld.o \
               efv9-dummy.o degen9-nondegenerate.o potin9-example.o
###############################################################################
#
# list of targets (which are all executables) =================================
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
USER-PROGS = 
###############################################################################
#
# rules for making executables ================================================
#
.PHONY: $(PROGS) $(USER-PROGS)
#$(PROGS) : %: $(EXECDIR)/%

$(EXECDIR)/molscat-basic     : $(addprefix $(OBJDIR)/,main1M.o \
                                                      $(CORE_MOL) \
                                                      $(REPLACEABLE) \
                                                      $(BASE9-UNUSED) \
                                                      $(POT-BASIC))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/molscat-H2_H2     : $(addprefix $(OBJDIR)/,main1M.o \
                                                      $(CORE_MOL) \
                                                      $(REPLACEABLE) \
                                                      $(BASE9-UNUSED) \
                                                      $(POT-H2H2))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/molscat-Rg_HX     : $(addprefix $(OBJDIR)/,main1M.o \
                                                      $(CORE_MOL) \
                                                      $(REPLACEABLE) \
                                                      $(BASE9-UNUSED) \
                                                      $(POT-ARHX))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/molscat-Rg_CO2    : $(addprefix $(OBJDIR)/,main1M.o \
                                                      $(CORE_MOL) \
                                                      $(REPLACEABLE) \
                                                      $(BASE9-UNUSED) \
                                                      $(POT-ARCO2))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/molscat-Rg_H2     : $(addprefix $(OBJDIR)/,main1M.o \
                                                      $(CORE_MOL) \
                                                      $(REPLACEABLE) \
                                                      $(BASE9-UNUSED) \
                                                      $(POT-RGH2))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/molscat-Ar_CH4    : $(addprefix $(OBJDIR)/,main1M.o \
                                                      $(CORE_MOL) \
                                                      $(REPLACEABLE) \
                                                      $(BASE9-UNUSED) \
                                                      $(POT-ARCH4))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/molscat-Mg_NH     : $(addprefix $(OBJDIR)/,main1M.o \
                                                      $(CORE_MOL) \
                                                      $(REPLACEABLE) \
                                                      $(BASE9-MGNH) \
                                                      $(POT-MGNH))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/molscat-Rb2       : $(addprefix $(OBJDIR)/,main1M.o \
                                                      $(CORE_MOL) \
                                                      $(REPLACEABLE) \
                                                      $(BASE9-ALK2) \
                                                      $(POT-RB2))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/bound-basic       : $(addprefix $(OBJDIR)/,main1M.o \
                                                      $(CORE_BND) \
                                                      $(REPLACEABLE) \
                                                      $(BASE9-UNUSED) \
                                                      $(POT-BASIC))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/bound-Rg_HX       : $(addprefix $(OBJDIR)/,main1M.o \
                                                      $(CORE_BND) \
                                                      $(REPLACEABLE) \
                                                      $(BASE9-UNUSED) \
                                                      $(POT-ARHX))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/bound-Rb2         : $(addprefix $(OBJDIR)/,main1M.o \
                                                      $(CORE_BND) \
                                                      $(REPLACEABLE) \
                                                      $(BASE9-ALK2) \
                                                      $(POT-RB2))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/bound-H2_H2       : $(addprefix $(OBJDIR)/,main1M.o \
                                                      $(CORE_BND) \
                                                      $(REPLACEABLE) \
                                                      $(BASE9-UNUSED) \
                                                      $(POT-H2H2))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/bound-Rg_H2       : $(addprefix $(OBJDIR)/,main1M.o \
                                                      $(CORE_BND) \
                                                      $(REPLACEABLE) \
                                                      $(BASE9-UNUSED) \
                                                      $(POT-RGH2))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/bound-Ar_CH4      : $(addprefix $(OBJDIR)/,main1M.o \
                                                      $(CORE_BND) \
                                                      $(REPLACEABLE) \
                                                      $(BASE9-UNUSED) \
                                                      $(POT-ARCH4))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/bound-Rg_CO2      : $(addprefix $(OBJDIR)/,main1M.o \
                                                      $(CORE_BND) \
                                                      $(REPLACEABLE) \
                                                      $(BASE9-UNUSED) \
                                                      $(POT-ARCO2))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/bound-Mg_NH       : $(addprefix $(OBJDIR)/,main1M.o \
                                                      $(CORE_BND) \
                                                      $(REPLACEABLE) \
                                                      $(BASE9-MGNH) \
                                                      $(POT-MGNH))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/field-basic       : $(addprefix $(OBJDIR)/,main1M.o \
                                                      $(CORE_FLD) \
                                                      $(REPLACEABLE) \
                                                      $(BASE9-UNUSED) \
                                                      $(POT-BASIC))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/field-Rb2         : $(addprefix $(OBJDIR)/,main1M.o \
                                                      $(CORE_FLD) \
                                                      $(REPLACEABLE) \
                                                      $(BASE9-ALK2) \
                                                      $(POT-RB2))
	$(LINK.f) $^ $(LIBS) -o $@
$(EXECDIR)/field-Mg_NH       : $(addprefix $(OBJDIR)/,main1M.o \
                                                      $(CORE_FLD) \
                                                      $(REPLACEABLE) \
                                                      $(BASE9-MGNH) \
                                                      $(POT-MGNH))
	$(LINK.f) $^ $(LIBS) -o $@

SHELL = /bin/sh

############################## Section 1 #####################################
# Set up all the variables that govern compilation

ifeq ($(MOLSCAT_COMPILER),)
  Compiler = gfortran
else
  Compiler = $(MOLSCAT_COMPILER)
endif

LIBS = -llapack -lblas
ifeq ($(MOLSCAT_COMPILER),gfortran)
  FFLAGS = -g --std=legacy
endif

ifneq ($(MOLSCAT_LIBS),)
  LIBS = $(MOLSCAT_LIBS)
endif
ifneq ($(MOLSCAT_FLAGS),)
  FFLAGS = $(FFLAGS) $(MOLSCAT_FLAGS)
endif

ifeq ($(MOLSCAT_OBJDIR),)
  OBJDIR = object-$(Compiler)
else
  OBJDIR = $(MOLSCAT_OBJDIR)
endif

ifeq ($(MOLSCAT_EXECDIR),)
  EXECDIR = bin
else
  EXECDIR = $(MOLSCAT_EXECDIR)
endif

LINK.f = $(Compiler)
COMPILE.f = $(Compiler) -c

# If you forget to type an executable name, this is what will be made.  You
# can change it by setting the environment variable MOLSCAT_EXECUTABLE
#
ifeq ($(MOLSCAT_EXECUTABLE),)
  .DEFAULT_GOAL := molscat-basic
else
  .DEFAULT_GOAL := $(MOLSCAT_EXECUTABLE)
endif

# Clear out the list of suffixes and then add the suffixes that define
# pattern-matching default rules
.SUFFIXES:
.SUFFIXES: .f .f90 .o

# Setting a list of object files for each of the common modules.
# These object files are dependent on the .mod file and the object file
# for that module.
#
# Object codes that reference the angles module
angles_DEPS = iosbin.o iosout.o potenl.o vrtp-extpot_1ang.o                    \
      vrtp-extpot_2ang.o vrtp-H2_H2.o vrtp-Rg_HX-eta.o

$(addprefix $(OBJDIR)/,$(angles_DEPS)): angles.mod $(OBJDIR)/angles_module.o

# Object codes that reference the basis_data module
basis_data_DEPS = asrot.o base9-1S_3Sigma_cpld.o                               \
      base9-alk_alk_ucpld.o base9-skeleton.o base.o bd.driver.o                \
      fld.driver.o getlev.o iosbin.o iosout.o iospb.o j6to4.o mol.driver.o     \
      prbr.o rasymu.o restrt.o set4.o set6c.o set6.o setbas.o sproc.o          \
      surbas.o thresh.o

$(addprefix $(OBJDIR)/,$(basis_data_DEPS)): basis_data.mod $(OBJDIR)/basis_data_module.o

# Object codes that reference the efvs module
efvs_DEPS = base9-1S_3Sigma_cpld.o base9-alk_alk_ucpld.o                       \
      base9-skeleton.o base.o bckgrd-dummy.o bd.driver.o efv9-dummy.o          \
      fld.driver.o header.o ldrdwt.o locpol.o mol.driver.o prnsum.o            \
      rdsigu.o restrt.o scapot-default.o setefv.o sproc.o wvhead.o

$(addprefix $(OBJDIR)/,$(efvs_DEPS)): efvs.mod $(OBJDIR)/efvs_module.o

# Object codes that reference the pair_state module
pair_state_DEPS = base9-skeleton.o base.o bd.driver.o fld.driver.o             \
      j6to4.o mol.driver.o rasymu.o set4.o set6c.o set6.o set6i.o setbas.o     \
      sig6.o surbas.o

$(addprefix $(OBJDIR)/,$(pair_state_DEPS)): pair_state.mod $(OBJDIR)/pair_state_module.o

# Object codes that reference the physical_constants module
physical_constants_DEPS = base9-1S_3Sigma_cpld.o                               \
      base9-alk_alk_ucpld.o bd.driver.o chckru.o eavg.o ecnv.o                 \
      extpot-Ar_CH4.o extpot-Ar_CO2.o extpot-Rg_HX.o fld.driver.o              \
      Hannover_module.o Hdata-Rb2-2010.o mol.driver.o surbas.o                 \
      vstar-Hannover.o

$(addprefix $(OBJDIR)/,$(physical_constants_DEPS)): physical_constants.mod $(OBJDIR)/physical_constants_module.o

# Object codes that reference the potential module
potential_DEPS = aiprop.o base9-1S_3Sigma_cpld.o                               \
      base9-alk_alk_ucpld.o base9-skeleton.o base.o bd.driver.o chckru.o       \
      cheint.o cntrct.o dermat.o dvprop.o findrm.o fld.driver.o getlev.o       \
      gnham.o hammat.o header.o ldprop.o maprop.o mdprop.o mgprop.o            \
      mol.driver.o multop.o odprop.o potenl.o potenl-Rg_H2.o restrt.o          \
      scwave.o setefv.o sproc.o sumlam.o thresh.o vstar-Hannover.o vvprop.o    \
      wkb.o yinit.o ytrans.o

$(addprefix $(OBJDIR)/,$(potential_DEPS)): potential.mod $(OBJDIR)/potential_module.o

# Object codes that reference the sizes module
sizes_DEPS = angles_module.o basis_data_module.o bd.driver.o eavg.o            \
      fld.driver.o iosbin.o iosout.o j6to4.o locpol.o mol.driver.o             \
      potential_module.o prbr.o restrt.o setefv.o

$(addprefix $(OBJDIR)/,$(sizes_DEPS)): sizes.mod $(OBJDIR)/sizes_module.o

# Setting a list of object files for other modules that are referenced
# by source code files different from the source code file that creates
# the module.
# If the module is created in more than one place, we do not include
# the object files in the list of prerequisites.
#
# Object codes that reference the pot_data_hannover module
pot_data_hannover_DEPS = Hdata-Rb2-2010.o vstar-Hannover.o

$(addprefix $(OBJDIR)/,$(pot_data_hannover_DEPS)): pot_data_hannover.mod $(OBJDIR)/Hannover_module.o

# set up general variables ===================================================
 COMMON_MODULES = sizes_module.o angles_module.o basis_data_module.o \
                  efvs_module.o pair_state_module.o potential_module.o \
                  physical_constants_module.o

 GCLOCK        = gclock.o gdate.o
 DIAG          = diagvc.o diagvl.o
 LA_UTILS      = dgemul.o syminl-alloc.o dsyneg.o syminv.o
 UTILITIES     = dsyfil.o matprn.o
 ANG_MOM       = j3j000.o j6j.o j9j.o parsgn.o sixj.o threej.o thrj.o xninej.o

 YTRANS        = chnsrt.o ev_ord.o gtwvec.o multop.o ytrans.o
 PROPS_LD      = aiprop.o airymp.o corr.o ldprop.o maprop.o maxmgv.o mdprop.o  \
                 mgprop.o outmat.o potent.o scairy.o spropn.o trnsfm.o         \
                 trnsp.o waveig.o yinit.o
 BASE_FILES    = asrot.o check6.o couple.o cpl21.o cpl22.o cpl23.o cpl25.o     \
                 cpl3.o cpl4.o cplout.o dmsym.o esymtp.o find.o fsymtp.o       \
                 gsymtp.o idpart.o ipasym.o j6to4.o mcgcpl.o order.o qsymtp.o  \
                 rasymu.o rsymtp.o set4.o set6.o set6c.o setbas.o surbas.o
#BASE          = $(BASE_FILES) base.o
#INOLLS_BASE   = $(BASE_FILES) base.inolls.o
 CORE          = brent.o calck.o chckru.o cheint.o drcalc.o drset.o ecnv.o     \
                 ecnvx.o findrm.o gnham.o idpchk.o main.o progms.o propst.o    \
                 prprop.o setefv.o stsrch.o thresh.o thrlst.o wvbksb.o wvhead.o
 W_MATRIX      = dermat.o hammat.o sumlam.o

 MOL           = calca.o convrg.o eavg.o epsum.o erange.o findrx.o getlev.o    \
                 ldrdwt.o locpol.o mol.prpini.o nexte.o header.o rdsigu.o      \
                 psiasy.o scctrl.o scpset.o scwave.o sproc.o
 PROPS_OTHER   = delrd.o dvprop.o pert1.o pert2.o odprop.o rmprop.o sgnchk.o   \
                 stabil.o vvprop.o wkb.o
 IOS           = dasize.o get102.o gaussp.o gasleg.o iosbin.o iosclc.o         \
                 iosdrv.o iosout.o iospb.o isutp.o plm.o set6i.o sig6.o yrr.o  \
                 zbes.o
 PRBR          = prbr.o prbr3.o
 RESTRT        = restrt.o sread.o
 YTOS          = ksym.o ktos.o rbesjy.o rbessk.o ytok.o ytokg.o

 FLD           = bdctrl.o bdpset.o iosbin-dummy.o bd.prpini.o evmtch.o prnsum.o

 BD            = cntrct.o $(FLD)


 ALL_CODE      = $(COMMON_MODULES) $(GCLOCK) $(DIAG) $(LA_UTILS) $(UTILITIES)  \
                 $(ANG_MOM) $(YTRANS) $(PROPS_LD) $(CORE) $(W_MATRIX)  \
                 $(MOL_ABSORB)
 CORE_ALL      = $(ALL_CODE) $(BASE_FILES)

# set up general variables needed for each program ===========================
 MOL_SUBS      = $(CORE_ALL) $(MOL) $(PROPS_OTHER) $(IOS) $(PRBR) $(RESTRT)    \
                 $(YTOS)
 FLD_SUBS      = $(CORE_ALL) $(FLD)
 BD_SUBS       = $(CORE_ALL) $(BD)
#
# set up dependencies needed for each program ===============================
#
 CORE_MOL      = $(MOL_SUBS) base.o mol.driver.o
 CORE_BD       = $(BD_SUBS)  base.o bd.driver.o
 CORE_FLD      = $(FLD_SUBS) base.o fld.driver.o
 INOLLS_MOL    = $(MOL_SUBS) i_nolls_module.o base.inolls.o mol.driver.inolls.o
 INOLLS_BD     = $(BD_SUBS)  i_nolls_module.o base.inolls.o bd.driver.inolls.o
 INOLLS_FLD    = $(FLD_SUBS) i_nolls_module.o base.inolls.o fld.driver.inolls.o

# Normal version of potential scaling
 SCALING       = reiniv-dummy.o scapot-default.o
 REPLACEABLE   = bckgrd-dummy.o pertrb.o $(SCALING)

# basic set of potenl and base9 routine
 POTUTILS      = gasleg.o gaushp.o gaussp.o herm.o hrecur.o plm.o yrr.o zbes.o
 POTENL-GP     = potenl.o $(POTUTILS)
 MAKEV         = vrtp-dummy.o vstar-dummy.o
#
 POT-BASIC     = $(POTENL-GP) $(MAKEV)
#
 EFV9DUM       = efv9-dummy.o
 THRSH9DUM     = thrsh9-dummy.o
 POTIN9DUM     = potin9-example.o
 DEGEN9DUM     = degen9-nondegenerate.o
 B9UTILS       = $(EFV9DUM) $(THRSH9DUM) $(POTIN9DUM) $(DEGEN9DUM)
 B9UTILS-EVF9     =            $(THRSH9DUM) $(POTIN9DUM) $(DEGEN9DUM)
 B9UTILS-THRSH9   = $(EFV9DUM)              $(POTIN9DUM) $(DEGEN9DUM)
 B9UTILS-POTIN9   = $(EFV9DUM) $(THRSH9DUM)              $(DEGEN9DUM)
 B9UTILS-DEGEN9   = $(EFV9DUM) $(THRSH9DUM) $(POTIN9DUM)
 EFV+THRSH_DUM    = $(EFV9DUM) $(THRSH9DUM)
 EFV+DEGEN_DUM    = $(EFV9DUM)                           $(DEGEN9DUM)
 THRSH+DEGEN_DUM  =            $(THRSH9DUM)              $(DEGEN9DUM)

 BASE9-UNUSED  = base9-skeleton.o $(B9UTILS-POTIN9)

#
# dependency lists for system-specific potentials =============================
#
 SUMLEG    = sumleg.o
 PLEG      = pleg.o
 PASLEG    = pasleg.o assleg.o
 DAMP      = damp.o
 
 V-Rg+HX   = $(VSTAR-DUM) vrtp-Rg_HX-eta.o
 
 VSTAR-DUM = vstar-dummy.o
 VRTP-DUM  = vrtp-dummy.o
 VRTP-1ang = vrtp-extpot_1ang.o
 VRTP-2ang = vrtp-extpot_2ang.o

# note that modules must be earlier in dependency lists than object codes
# that use them because they must be compiled first
#
 V-1ANG    = $(VRTP-1ang) $(VSTAR-DUM)
 V-2ANG    = $(VRTP-2ang) $(VSTAR-DUM)
 V-Hann      = $(VRTP-DUM) Hannover_module.o vstar-Hannover.o
 V-Mg+NH   = vstar-Mg_NH.o $(VRTP-DUM)
 
 DATA-RB2-Hann-2010 = Hdata-Rb2-2010.o
 
 EXTPOT-Rg+HX     = extpot-Rg_HX.o $(DAMP) $(SUMLEG)
 
 POT-H2H2  = $(POTENL-GP) vstar-H2_H2.o       vrtp-H2_H2.o
 POT-ARHX  = $(POTENL-GP) $(V-Rg+HX) $(EXTPOT-Rg+HX) $(PLEG) $(PASLEG)
 POT-ARCH4 = $(POTENL-GP) vstar-Rg_CH4_Buck.o $(VRTP-2ang) extpot-Ar_CH4.o
 POT-ARCO2 = $(POTENL-GP) $(V-1ANG) extpot-Ar_CO2.o $(SUMLEG)
 POT-MGNH  = $(POTENL-GP) $(V-Mg+NH)
 POT-RGH2  = potenl-Rg_H2.o $(POTUTILS)
 
 POT-RB2_Hann-2010 = $(POTENL-GP) $(V-Hann) $(DATA-RB2-Hann-2010)

#
# dependency lists for plug-in basis-set suites ===============================
#
 BASE9-UNUSED        = base9-skeleton.o $(B9UTILS-POTIN9)
 BASE9-ALK-ALK-UCPLD = base9-alk_alk_ucpld.o $(EFV+DEGEN_DUM)
 BASE9-1S+3Sigma     = base9-1S_3Sigma_cpld.o $(B9UTILS-THRSH9)

#
# lists of object files for each executable
#
 molscat-basic_OBJS  = $(CORE_MOL) $(REPLACEABLE) $(BASE9-UNUSED) $(POT-BASIC)
 bound-basic_OBJS    = $(CORE_BD)  $(REPLACEABLE) $(BASE9-UNUSED) $(POT-BASIC)
 field-basic_OBJS    = $(CORE_FLD) $(REPLACEABLE) $(BASE9-UNUSED) $(POT-BASIC)
 molscat-Ar_CH4_OBJS = $(CORE_MOL) $(REPLACEABLE) $(BASE9-UNUSED) $(POT-ARCH4)
 bound-Ar_CH4_OBJS   = $(CORE_BD)  $(REPLACEABLE) $(BASE9-UNUSED) $(POT-ARCH4)
 molscat-Rg_CO2_OBJS = $(CORE_MOL) $(REPLACEABLE) $(BASE9-UNUSED) $(POT-ARCO2)
 bound-Rg_CO2_OBJS   = $(CORE_BD)  $(REPLACEABLE) $(BASE9-UNUSED) $(POT-ARCO2)
 molscat-Rg_HX_OBJS  = $(CORE_MOL) $(REPLACEABLE) $(BASE9-UNUSED) $(POT-ARHX)
 bound-Rg_HX_OBJS    = $(CORE_BD)  $(REPLACEABLE) $(BASE9-UNUSED) $(POT-ARHX)
 molscat-H2_H2_OBJS  = $(CORE_MOL) $(REPLACEABLE) $(BASE9-UNUSED) $(POT-H2H2)
 bound-H2_H2_OBJS    = $(CORE_BD)  $(REPLACEABLE) $(BASE9-UNUSED) $(POT-H2H2)
 molscat-Rg_H2_OBJS  = $(CORE_MOL) $(REPLACEABLE) $(BASE9-UNUSED) $(POT-RGH2)
 bound-Rg_H2_OBJS    = $(CORE_BD)  $(REPLACEABLE) $(BASE9-UNUSED) $(POT-RGH2)
 molscat-Mg_NH_OBJS  = $(CORE_MOL) $(REPLACEABLE) $(BASE9-1S+3Sigma) $(POT-MGNH)
 bound-Mg_NH_OBJS    = $(CORE_BD)  $(REPLACEABLE) $(BASE9-1S+3Sigma) $(POT-MGNH)
 field-Mg_NH_OBJS    = $(CORE_FLD) $(REPLACEABLE) $(BASE9-1S+3Sigma) $(POT-MGNH)
 molscat-Rb2_OBJS    = $(CORE_MOL) $(REPLACEABLE) $(BASE9-ALK-ALK-UCPLD) \
                       $(POT-RB2_Hann-2010)
 bound-Rb2_OBJS      = $(CORE_BD)  $(REPLACEABLE) $(BASE9-ALK-ALK-UCPLD) \
                       $(POT-RB2_Hann-2010)
 field-Rb2_OBJS      = $(CORE_FLD) $(REPLACEABLE) $(BASE9-ALK-ALK-UCPLD) \
                       $(POT-RB2_Hann-2010)
#
# list of executables
#
 PUBLIC-PROGS = molscat-basic bound-basic field-basic \
                molscat-Ar_CH4 bound-Ar_CH4 \
                molscat-Rg_CO2 bound-Rg_CO2 \
                molscat-Rg_HX bound-Rg_HX \
                molscat-H2_H2 bound-H2_H2 \
                molscat-Rg_H2 bound-Rg_H2 \
                molscat-Mg_NH bound-Mg_NH field-Mg_NH \
                molscat-Rb2 bound-Rb2 field-Rb2
 
 PROGS += $(PUBLIC-PROGS)

 USER-PROGS = 

 PROGS += $(USER-PROGS)

#
# Set up the dependency of each executable on its list of object files
#
define PROGS_template =
  $(1): $$(addprefix $$(OBJDIR)/,$$($(1)_OBJS))
  vpath $(1) $(EXECDIR)
endef

#
# Cycle over the executables setting up their dependencies
#
$(foreach prog,$(PROGS),$(eval $(call PROGS_template,$(prog))))


############################## Section 5 #####################################
#
# Set up rules
#

.PHONY: all

all: $(PROGS)

#
# Rule to make any executable
#
$(PROGS): | $(EXECDIR)
	$(LINK.f) $^ $(LIBS) -o $(EXECDIR)/$@
	@rm -f $(addprefix $(OBJDIR)/,$($(@F)-remove)) 2>&1 > /dev/null

#
# Rule to make the object directory
#
$(OBJDIR):
	@mkdir $(OBJDIR)

#
# Rule to make the executable directory
#
$(EXECDIR):
ifneq ($(EXECDIR),$(OBJDIR))
	@mkdir $(EXECDIR)
endif

#
# Default rules based on suffixes
#
$(OBJDIR)/%.o: %.f | $(OBJDIR)
	$(COMPILE.f) $(FFLAGS) $< -o $@

$(OBJDIR)/%.o: %.f90 | $(OBJDIR)
	$(COMPILE.f) $(FFLAGS) $< -o $@


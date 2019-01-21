# MOLSCAT

MOLSCAT is a general-purpose package for performing non-reactive quantum scattering calculations for atomic and molecular collisions using coupled-channel methods. Simple atom-molecule and molecule-molecule collision types are coded internally and additional ones may be handled with plug-in routines. Plug-in routines may include external magnetic, electric or photon fields (and combinations of them). Simple interaction potentials may be specified in input data and more complicated ones may be handled with plug-in routines.

BOUND is a general-purpose package for performing calculations of bound-state energies in weakly bound atomic and molecular systems using coupled-channel methods. It solves the same sets of coupled equations as MOLSCAT, and can use the same plug-in routines if desired, but applies bound-state boundary conditions.

FIELD is a development of BOUND that locates values of external fields at which a bound state exists with a specified energy.  One important use is to locate the positions of magnetically tunable Feshbach resonance positions in ultracold collisions.

The three programs have different applications, but they use closely related methods and share many subroutines. They are therefore released with a single code base and are documented together.

The programs have built-in capabilities to generate the coupled equations for
- Atom + linear rigid rotor;
- Atom + vibrating diatom;
- Linear rigid rotor + linear rigid rotor;
- Atom + symmetric top;
- Atom + asymmetric top;
- Asymmetric top + linear molecule;
- Atom + rigid corrugated surface: diffractive (elastic) scattering.

For these cases, the programs implement both close-coupling calculations, with no dynamical approximations, and a variety of approximate methods including the coupled states and helicity decoupling approximations. MOLSCAT can loop over total angular momentum (partial wave) to calculate elastic and inelastic integral cross sections and spectroscopic line-shape cross sections. Post-processors are available that read MOLSCAT S-matrix files and calculate differential cross sections, transport, relaxation and Senftleben-Beenakker cross sections, and fit the parameters of scattering resonances.

The programs provide an interface for a plug-in basis-set suite to set up other sets of coupled equations. This capability has been widely used in bound-state and scattering calculations related to the spectroscopy of Van der Waals complexes and in ultracold atomic and molecular collisions. Basis-set suites can be programmed to handle a wide variety of interacting species and to take account of one or more external fields, such as electric, magnetic and photon fields. Two such suites are included in this distribution, for
- closed-shell atom + triplet-Sigma diatom in a magnetic field;
- Alkali-metal atom + alkali-metal atom in a magnetic field, including hyperfine structure.

For low-energy scattering, MOLSCAT can calculate scattering lengths and effective ranges and can locate and characterize scattering resonances as a function of an external variable such as the magnetic field.

Further information on the programs is available at
- MOLSCAT: https://arxiv.org/abs/1811.09584
- BOUND and FIELD: https://arxiv.org/abs/1811.09111

Full program documentation is in the file molscat_bound_field.pdf

The programs are free software, distributed under the terms of the GNU General Public License, Version 3.

The programs are written in FORTRAN 90, though with many features from older versions of FORTRAN. They have been tested with common FORTRAN compilers including gfortran and ifort. 

--------------------------------------------------------------------------------

This release contains the following files:

1. Fortran source code files ending in .f or .f90, which contain the code needed to construct basic executables for MOLSCAT, BOUND or FIELD, together with some additional source code files needed to produce variants of MOLSCAT, BOUND or FIELD with special-purpose basis-set suites and/or potential routines, used for the example calculations outlined in the documentation.  These are contained in the directory source_code
2. A makefile (GNUmakefile) which may be used to make the MOLSCAT, BOUND or FIELD executables used for the example calculations.  This is also contained in the directory source_code
3. All the input data files  used for example calculations in chapter 3 ('Using the programs: a basic guide', sections 3.7, 3.8 and 3.9) and for the example calculations in chapter 13 ('Example input and output files') of the documentation, together with associated output files.  The input files are contained in the directory examples/input/ and the output files are contained in the directory examples/output/.
4. Data files used to construct the potential for Mg + NH in the executables molscat-Mg_NH, bound-Mg_NH and field-Mg_NH and for H2 + H2 in the executables bound-H2_H2 and molscat-H2_H2.  These are contained in the directory data

--------------------------------------------------------------------------------

Editing the makefile:

The programs are supplied with a makefile named GNUmakefile.  This is designed for use with GNU make (gmake), which is standard in most Linux distributions.

GNUmakefile first needs to be adapted for your own working environment.  MOLSCAT, BOUND and FIELD executables make use of BLAS and LAPACK library routines.  If optimised versions of these library routines are available, include the libraries containing them in the variable LIBS.  Otherwise, download the requisite routines from the Netlib repository and include the object code files in the variable LIBUTILS.

The BLAS routines used are:

daxpy     dcopy     dgemm     dgemv     dgesv     dscal     dswap     dsymm     idamax

and the LAPACK routines are:

dgesv     dlamch    dsyevr    dsyr2k    dsytrf    dsytri    ilaenv    lsame

Many of the LAPACK routines call other BLAS and LAPACK routines.

GNUmakefile sets the compiler (in the variable Compiler) to be gfortran.  If gfortran is not available, or you prefer to use another compiler, you will need to change this.

The optimisation level is set at 0 and traceback is switched on.  These are set in the variable COMPILE.f.  Once you have assured yourself that the code works you may wish to increase the optimisation level, and switch the traceback off.

The executables all make use of date and time routines GDATE, GTIME and GCLOCK. GDATE and GTIME use the f90 intrinsic subroutine date_and_time, and GCLOCK uses the f90 intrinsic subroutine cpu_time.  If these are not available, they can be replaced with local routines.  In the last resort, GDATE and GTIME should return blank character variables and GCLOCK a double precision value of 0.D0.

--------------------------------------------------------------------------------

Making the executables:

GNUmakefile contains rules to make executables for several variants of MOLSCAT, BOUND or FIELD, with filenames prefixed by molscat-, bound- or field- as appropriate.

They may be made with Linux commands such as
make molscat-H2_H2
or analogous commands for other variants of the programs.  The list of executables for which GNUmakefile contains rules is contained in the variable PROGS.

GNUmakefile places the object files in the directory named in the variable OBJDIR and the executables in the directory named in the variable EXECDIR.  These are set to the current directory (.), but you may wish to change this, in which case you need to change the directories named in these variables and ensure that they already exist. You will also need to remove the hash character at the beginning of the line

#$(PROGS) : %: $(EXECDIR)/%

(this line is commented out in GNUmakefile because otherwise make will produce a warning about a circular dependency if the variable EXECDIR is set to the current directory.)

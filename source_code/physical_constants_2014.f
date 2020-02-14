      module physical_constants
C  Copyright (C) 2020 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      implicit none

c  These are the CODATA 2014 recommended values, inserted 13-10-16
      integer,parameter:: NIST_year=2014
      double precision,parameter::
     &bohr_in_SI                       = 0.52917721067d-10,
     &hartree_in_SI                    = 4.359744650d-18,
     &hartree_in_inv_cm                = 2.194746313702d5,
     &Hz_in_inv_cm                     = 3.335640951d-11,
     &eV_in_inv_cm                     = 8.065544005d3,
     &K_in_inv_cm                      = 0.69503457d0,
     &g_e                              = -2.00231930436182d0,
     &fine_structure_constant          = 7.2973525664d-3,
     &inverse_fine_structure_constant  = 137.035999139d0,
     &bohr_magneton_in_inv_cm_per_T    = 0.4668644814d0,
     &nuclear_magneton_in_inv_cm_per_T = 2.542623432d-4,
     &Planck_constant_in_SI            = 6.626070040d-34,
     &hbar_in_SI                       = 1.054571800d-34,
     &speed_of_light_in_SI             = 2.99792458d8,
     &pi_value                         = 3.141592653589793d0,
     &atomic_mass_constant             = 1.660539040d-27,
     &cal_in_SI                        = 4.184d0,
     &Joule_in_inv_cm                  = 5.034116651d22,
     &Avogadro_number                  = 6.022140857d23

      double precision,parameter::
     &Angstrom_in_SI = 1d-10,
     &kilo_in_SI     = 1d3,
     &Mega_in_SI     = 1d6,
     &Giga_in_SI     = 1d9,
     &Gauss_in_SI    = 1d4,
     &centi_in_SI    = 1d2,
     &erg_in_SI      = 1d-7

      double precision,parameter::
     &bohr_to_Angstrom           = bohr_in_SI/Angstrom_in_SI,
     &MHz_in_inv_cm              = Hz_in_inv_cm*Mega_in_SI,
     &GHz_in_inv_cm              = Hz_in_inv_cm*Giga_in_SI,
     &kJ_per_mol_in_inv_cm       = Joule_in_inv_cm*kilo_in_SI/
     &                             Avogadro_number,
     &kcal_per_mol_in_inv_cm     = kJ_per_mol_in_inv_cm*cal_in_SI,
     &Debye_in_SI                = Hz_in_inv_cm/1d19,
     &Debye_Volt_metre_in_inv_cm = Debye_in_SI*Joule_in_inv_cm,
     &erg_in_inv_cm              = Joule_in_inv_cm*erg_in_SI,
     &bohr_magneton        = bohr_magneton_in_inv_cm_per_T/Gauss_in_SI,
     &speed_of_light_in_cm = speed_of_light_in_SI*centi_in_SI,
     &nuclear_magneton = nuclear_magneton_in_inv_cm_per_T/Gauss_in_SI,
     &bfct             = hbar_in_SI/(4d0*pi_value*speed_of_light_in_cm)/
     &                   (Angstrom_in_SI)**2/amu_in_SI

      end module physical_constants

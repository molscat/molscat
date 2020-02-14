      module physical_constants
C  Copyright (C) 2020 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      implicit none

c  These values are from the BIPM draft brochure.  Values inserted 27-03-2019
      integer,parameter:: NIST_year=2018
      double precision, parameter ::
     & speed_of_light_in_SI                    = 2.99792458d8,     ! m s^-1
     & Planck_constant_in_SI                   = 6.62607015d-34,   ! J s
     & Boltzmann_constant_in_SI                = 1.380649d-23,     ! J K^-1
     & Avogadro_constant_in_SI                 = 6.02214076d23,    ! mol^-1
     & Cs_hyperfine_transition_frequency_in_SI = 9.192631770d9,    ! Hz
     & elementary_charge_in_SI                 = 1.6021766340d-19, ! C
c  and these are other defined values
     & pi_value                                = 3.141592653589793d0,
     & hbar_in_SI                              = Planck_constant_in_SI
     &                                           /2d0/pi_value,
     & cal_in_SI                               = 4.184d0     ! J

      double precision, parameter ::
     & inverse_fine_structure_constant = 137.035999084d0,    ! 2018 CODATA value
     & fine_structure_constant         = 7.2973525693d-3,    ! 2018 CODATA value
c  electron spin g-factor (2018 CODATA value)
     & g_e                             =-2.00231930436256d0,
c  atomic mass constant (2018 CODATA value)
     & atomic_mass_constant            = 1.66053906660d-27,
c  proton mass (2018 CODATA value)
     & proton_mass                     = 1.67262192369d-27,
c
c  electron mass (2018 CODATA value)
     & electron_mass                   = 9.1093837015d-31

c  These values are calculated from the above (the values recorded here are from 2018 CODATA tables)
      double precision,parameter::
c  bohr radius is defined as hbar/(electron_mass*speed_of_light*fine_structure_constant)
     & bohr_in_SI                       = 0.529177210903d-10,        ! m
c  hartree is defined as electron_mass*(speed_of_light*fine_structure_constant)**2
     & hartree_in_SI                    = 4.3597447222071d-18,       ! J
c  bohr magneton is defined as elementary_charge*hbar/(2d0*electron_mass)
c  (nuclear magneton is defined as bohr_magneton*electron_mass/proton_mass)
     & bohr_magneton_in_SI              = 927.40100783d-26,          ! J/T
     & electronvolt_in_SI               = elementary_charge_in_SI,   ! V
c  conversion factor from Joule to inverse cm
     & Joule_in_inv_cm                  = 1d-2/Planck_constant_in_SI ! cm^-1/J
     &                                        /speed_of_light_in_SI,
c  and these are energy equivalents (or prefactors for the magneton values) in inverse cm
     & hartree_in_inv_cm                = 2.1947463136320d5, ! hartree_in_SI*Joule_in_inv_cm
     & Hz_in_inv_cm                     = 3.335640951d-11,   ! 1d-2/speed_of_light_in_SI
     & eV_in_inv_cm                     = 8.065543937d3,     ! electronvolt_in_SI*Joule_in_inv_cm
     & K_in_inv_cm                      = 0.6950348004d0,    ! Boltzmann_constant_in_SI*Joule_in_inv_cm
     & bohr_magneton_in_inv_cm_per_T    = 0.46686447783d0,   ! bohr_magneton_in_SI*Joule_in_inv_cm
     & nuclear_magneton_in_inv_cm_per_T = 2.54262341353d-4   ! bohr_magneton_in_inv_cm_per_T*electron_mass/proton_mass

      double precision,parameter::
     & Angstrom_in_SI = 1d-10, ! m per Angstrom
     & kilo_in_SI     = 1d3,
     & Mega_in_SI     = 1d6,
     & Giga_in_SI     = 1d9,
     & Gauss_in_SI    = 1d4,
     & centi_in_SI    = 1d2,
     & erg_in_SI      = 1d-7

      double precision,parameter::
     & bohr_to_Angstrom           = bohr_in_SI/Angstrom_in_SI,          ! Angstrom
     & MHz_in_inv_cm              = Hz_in_inv_cm*Mega_in_SI,            ! cm^-1 per MHz
     & GHz_in_inv_cm              = Hz_in_inv_cm*Giga_in_SI,            ! cm^-1 per GHz
     & kJ_per_mol_in_inv_cm       = Joule_in_inv_cm*kilo_in_SI          ! cm^-1 per kJ/mol
     &                              /Avogadro_constant_in_SI,
     & kcal_per_mol_in_inv_cm     = kJ_per_mol_in_inv_cm*cal_in_SI,     ! cm^-1 per kcal/mol
     & Debye_in_SI                = Hz_in_inv_cm/1d19,                  ! D
     & Debye_Volt_metre_in_inv_cm = Debye_in_SI*Joule_in_inv_cm,        ! cm^-1 per D Vm
     & erg_in_inv_cm              = Joule_in_inv_cm*erg_in_SI,          ! cm^-1 per erg
     & bohr_magneton              = bohr_magneton_in_inv_cm_per_T       ! cm^-1/G
     &                              /Gauss_in_SI,
     & speed_of_light_in_cm       = speed_of_light_in_SI*centi_in_SI,   ! cm/s
     & nuclear_magneton           = nuclear_magneton_in_inv_cm_per_T    ! cm^-1/G
     &                              /Gauss_in_SI,
     & bfct                       = hbar_in_SI                          ! cm^-1
     &                              /(4d0*pi_value*speed_of_light_in_cm)
     &                              /(Angstrom_in_SI)**2
     &                              /atomic_mass_constant

      end module physical_constants

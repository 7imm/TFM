!-------------------------------------------------------------------------------
!    _\/  \/_      ______________________  ___
!     _\/\/_       ___  __/__  ____/__   |/  /
! _\_\_\/\/_/_/_   __  /  __  /_   __  /|_/ /
!  / /_/\/\_\ \    _  /   _  __/   _  /  / /
!     _/\/\_       /_/    /_/      /_/  /_/
!     /\  /\
!
! Copyright 2024 Timm Schultz
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!-------------------------------------------------------------------------------
MODULE tfm_density_tools
!-------------------------------------------------------------------------------
! Module: tfm_density_tools
!
! This module contains procedures commonly used by density related
! functions and subroutines.
!
! Dependencies: tfm_essentials, tfm_constants
!
! Interfaces:
!   tfm_density_arrhenius: Interface for functions
!                          "tfm_density_arrhenius_i" and
!                          "tfm_density_arrhenius_n".
!
! Functions:
!   tfm_density_arrhenius_i: Arrhenius equation for dimension 1.
!   tfm_density_arrhenius_n: Arrhenius equation for dimension n.
!
! Subroutines:
!   tfm_density_lin_interp: Simple linear interpolation.
!   tfm_density_bagmean: Bagmean from firn profile.
!-------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! modules used
  !-----------------------------------------------------------------------------
  USE, INTRINSIC :: iso_fortran_env, ONLY : dp => real64

  USE tfm_constants, ONLY : &
    MELT_TEMP,              &
    ICE_DENSITY,            &
    GAS_CONST

  !-----------------------------------------------------------------------------
  ! declarations
  !-----------------------------------------------------------------------------
  IMPLICIT NONE (TYPE, EXTERNAL)
  PRIVATE

  !-----------------------------------------------------------------------------
  ! public routines
  !-----------------------------------------------------------------------------
  PUBLIC ::                &
    tfm_density_arrhenius, &
    tfm_density_bagmean,   &
    tfm_density_lin_interp


  !-----------------------------------------------------------------------------
  ! interfaces
  !-----------------------------------------------------------------------------
  INTERFACE tfm_density_arrhenius
    MODULE PROCEDURE           &
      tfm_density_arrhenius_i, &
      tfm_density_arrhenius_n
  END INTERFACE tfm_density_arrhenius

  !-----------------------------------------------------------------------------
  ! routines of this module
  !-----------------------------------------------------------------------------
  CONTAINS


  FUNCTION tfm_density_arrhenius_i( &
  &  i,                             &
  &  factor,                        &
  &  activation_energy,             &
  &  temperature                    &
  ) RESULT(arrhenius)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_arrhenius_i
    !
    ! The function computes an arrhenius equation for a given factor,
    ! activation energy and temperature. The dimension of the input and
    ! output is always one.
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   i: Dummy variable needed for the interface (tfm_density_arrhenius).
    !   factor: Pre factor of the arrhenius equation.
    !   activation_energy: Activation energy of the arrhenius equation.
    !   tempreature: Temperature (K).
    !
    ! Result:
    !   arrhenius: Arrhenius Factor (dimension 1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: i

    REAL(dp), INTENT(IN) :: &
      factor,               &
      activation_energy,    &
      temperature

    INTEGER :: m

    REAL(dp) :: arrhenius

    !---------------------------------------------------------------------------

    m = 1 * i
    arrhenius = factor * exp(-activation_energy / (GAS_CONST * temperature))
  END FUNCTION tfm_density_arrhenius_i


  FUNCTION tfm_density_arrhenius_n( &
  &  n,                             &
  &  factor,                        &
  &  activation_energy,             &
  &  temperature                    &
  ) RESULT(arrhenius)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !--------------------------------------------------------------------------
    ! Function: tfm_density_arrhenius_n
    !
    ! The function computes an arrhenius equation for a given factor,
    ! activation energy and temperature. The dimension of the input and
    ! output is defined by input.
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   n: Dimension of temperature input and arrhenius factor output.
    !   factor: Pre factor of the arrhenius equation.
    !   activation_energy: Activation energy of the arrhenius equation.
    !   temperature: Temperature (K).
    !
    ! Result:
    !   arrhenius: Arrhenius factor (dimension n).
    !--------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: n

    REAL(dp), INTENT(IN) :: &
      factor,               &
      activation_energy

    REAL(dp), DIMENSION(n), INTENT(IN) :: temperature

    REAL(dp), DIMENSION(n) :: arrhenius

    !--------------------------------------------------------------------------

    arrhenius = factor * exp(-activation_energy / (GAS_CONST * temperature))
  END FUNCTION tfm_density_arrhenius_n


  SUBROUTINE tfm_density_lin_interp( &
  &  z0,                             &
  &  z1,                             &
  &  v0,                             &
  &  v1,                             &
  &  dz,                             &
  &  v                               &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Subroutine: tfm_density_lin_interp
    !
    ! Simple routine for linear interpolation between given values.
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   z0: First given x-value.
    !   z1: Second given x-value.
    !   v0: First fiven y-value.
    !   v1: Second given y-value.
    !   dz: Distance to value z0 at which to interpolate at.
    !   v - on input: Variable to store the result.
    !
    ! Result:
    !   v - on output: Value from linear interpolation at z0 + dz.
    !---------------------------------------------------------------------------

    REAL(dp), INTENT(IN) :: &
      z0,                   &
      z1,                   &
      v0,                   &
      v1,                   &
      dz

    REAL(dp), INTENT(INOUT) :: v

    !---------------------------------------------------------------------------

    v = v0 + ((v1 - v0) / (z1 - z0)) * dz
  END SUBROUTINE tfm_density_lin_interp


  SUBROUTINE tfm_density_bagmean( &
  &  nz,                          &
  &  depth,                       &
  &  density,                     &
  &  dz,                          &
  &  mz,                          &
  &  bagmean                      &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Subroutine: tfm_density_bagmean
    !
    ! The routine computes bagmean values of a given density profile.
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of arrays "depth" and "density".
    !   depth: Array containing the depth values of the density profile (m).
    !   density: Array containing the density values of the density
    !            profile (kg m**-3).
    !   dz: Size of the "window" the bagmean values are computed from (m).
    !   mz: Length of the bagmean density profile (total length of
    !       the profile is mz*dz)
    !   bagmean - on input: Variable to store the bagmean profile.
    !
    ! Result:
    !   bagmean - on output: Array of dimension "mz" contaning the depth
    !                        (1,:) values (m) of the bagmean profile and
    !                        the corresponding density (2,:)
    !                        values (kg m**-3).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: &
      nz,                  &
      mz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      density

    REAL(dp), INTENT(IN) :: dz

    REAL(dp), DIMENSION(2,mz), INTENT(INOUT) :: bagmean

    INTEGER :: &
      n,       &
      m

    REAL(dp) :: d

    !---------------------------------------------------------------------------

    bagmean = 0.0_dp

    n = nz
    DO m = 1, mz, 1

      d = 0.0_dp
      bagmean(1,mz-m+1) = depth(nz) - (m * dz) + (0.5_dp * dz)

      DO WHILE ( depth(n) > (depth(nz) - (m * dz)) )
        bagmean(2,mz-m+1) = bagmean(2,mz-m+1) + density(n)
        n = n - 1
        d = d + 1.0_dp
      END DO

      bagmean(2,mz-m+1) = bagmean(2,mz-m+1) / d
    END DO
  END SUBROUTINE tfm_density_bagmean
END MODULE tfm_density_tools



MODULE tfm_density_stress
!-------------------------------------------------------------------------------
! Module: tfm_density_stress
!
! The module contains common stress related procedures.
!
! Dependencies: tfm_essentials, tfm_constants
!
! Functions:
!  tfm_density_boyleMariotte: Boyle-Mariotte factor from absolute
!                             density.
!  tfm_rel_density_boyleMariotte: Boyle-Mariotte factor from
!                                 relative density
!
! Subroutines:
!  tfm_density_computeStress: stres from depth and density
!-------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! modules used
  !-----------------------------------------------------------------------------
  USE, INTRINSIC :: iso_fortran_env, ONLY : dp => real64

  USE tfm_constants, ONLY : &
    CLOSEOFF_DENSITY,       &
    ICE_DENSITY,            &
    ACC_GRAVITY

  !-----------------------------------------------------------------------------
  ! declarations
  !-----------------------------------------------------------------------------
  IMPLICIT NONE (TYPE, EXTERNAL)
  PRIVATE

  !-----------------------------------------------------------------------------
  ! public routines
  !-----------------------------------------------------------------------------
  PUBLIC ::                    &
    tfm_density_computestress, &
    tfm_density_rel_boylemariotte

  !-----------------------------------------------------------------------------
  ! routines of this module
  !-----------------------------------------------------------------------------
  CONTAINS


  SUBROUTINE tfm_density_computeStress( &
  &  nz,                                &
  &  depth,                             &
  &  density,                           &
  &  stress                             &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Subroutine tfm_density_computeStress
    !
    ! Routine to compute the stress along a given firn profile, by
    ! integration of the density.
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "depth", "density", and "stress".
    !   depth: Depth along the firn profile (m).
    !   density: Density along the firn profile (kg m**-3).
    !   stress - on input: Variable to store the result.
    !
    ! Result:
    !   stress - on output: Stress computed from "depth" and
    !                       "density" (Pa).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: &
      nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      density

    REAL(dp), DIMENSION(nz), INTENT(INOUT) :: &
      stress

    INTEGER                 :: n
    REAL(dp), DIMENSION(nz) :: dz

    !---------------------------------------------------------------------------

    dz(nz) = 0.0_dp
    dz(1:nz-1) =  depth(2:nz) - depth(1:nz-1)

    stress = (dz * density * ACC_GRAVITY)
    DO n = nz - 1, 1, -1
      stress(n) = stress(n) + stress(n+1)
    END DO
  END SUBROUTINE tfm_density_computeStress


  FUNCTION tfm_density_boyleMariotte( &
  &  density                          &
  ) RESULT(boyle_mariotte)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_boyleMariotte
    !
    ! Boyle-Mariotte equation to compute the pore pressure within bubbly
    ! ice from given absoulte density.
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   density: Absolute density (kg m**-3).
    !
    ! Result:
    !   boyle_mariotte: Boyle-Mariotte factor. The pore pressure is the
    !                   product of this factor and the pressure at pore
    !                   close off (1).
    !---------------------------------------------------------------------------

    REAL(dp), INTENT(IN) :: density
    REAL(dp)             :: boyle_mariotte

    !---------------------------------------------------------------------------

    boyle_mariotte = (                                &
    &  (density * (ICE_DENSITY - CLOSEOFF_DENSITY))   &
    &  / (CLOSEOFF_DENSITY * (ICE_DENSITY - density)) &
    )
  END FUNCTION tfm_density_boyleMariotte


  FUNCTION tfm_density_rel_boyleMariotte( &
  &  rel_density                          &
  ) RESULT(rel_boyle_mariotte)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_rel_boyleMariotte
    !
    ! Boyle-Mariotte equation to compute the pore pressure within bubbly
    ! ice from given relative density.
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   density: Relative density (1).
    !
    ! Result:
    !   boyle_mariotte: Boyle-Mariotte factor. The pore pressure is the
    !   product of this factor and the pressure at pore close off (1).
    !---------------------------------------------------------------------------

    REAL(dp), INTENT(IN) :: rel_density
    REAL(dp)             :: rel_boyle_mariotte

    !---------------------------------------------------------------------------

    rel_boyle_mariotte = (                                           &
    &  rel_density * (1.0_dp - (CLOSEOFF_DENSITY / ICE_DENSITY))     &
    &  / ((CLOSEOFF_DENSITY / ICE_DENSITY) * (1.0_dp - rel_density)) &
    )
  END FUNCTION tfm_density_rel_boyleMariotte
END MODULE tfm_density_stress


MODULE tfm_density_processes
! ------------------------------------------------------------------------------
! Module: tfm_density_processes
!
! The module contains functions for process based firn densification
! modelling.
!
! Dependencies: tfm_constants, tfm_density_tools
!
! Functions:
!   tfm_density_NabarroHerringCreep
!   tfm_density_NabarroHerringCreepMod
!   tfm_density_CobleCreep
!   tfm_density_DislocationCreep
!   tfm_density_DislocationCreepMod
!   tfm_density_GrainBoundarySliding
! ------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! modules used
  !-----------------------------------------------------------------------------
  USE, INTRINSIC :: iso_fortran_env, ONLY : dp => real64

  USE tfm_constants, ONLY :               &
    ACTIVATION_ENERGY_BOUNDARY_DIFFUSION, &
    BOLTZMANN,                            &
    BURGERS_VECTOR,                       &
    ICE_DENSITY,                          &
    PRE_FACTOR_BOUNDARY_DIFFUSION,        &
    VOLUME_H2o,                           &
    ICE_N,                                &
    PRE_FACTOR_LATTICE_DIFFUSION,         &
    ACTIVATION_ENERGY_LATTICE_DIFFUSION

  USE tfm_density_tools, ONLY : &
    tfm_density_arrhenius

  !-----------------------------------------------------------------------------
  ! declarations
  !-----------------------------------------------------------------------------
  IMPLICIT NONE (TYPE, EXTERNAL)
  PRIVATE

  !-----------------------------------------------------------------------------
  ! public routines
  !-----------------------------------------------------------------------------
  PUBLIC ::                             &
    tfm_density_CobleCreep,             &
    tfm_density_DislocationCreep,       &
    tfm_density_DislocationCreepMod,    &
    tfm_density_NabarroHerringCreep,    &
    tfm_density_NabarroHerringCreepMod, &
    tfm_density_GrainBoundarySliding

  !-----------------------------------------------------------------------------
  ! routines of this module
  !-----------------------------------------------------------------------------
  CONTAINS


  FUNCTION tfm_density_NabarroHerringCreep( &
  &  temperature,                           &
  &  grain_radius,                          &
  &  driving_force                          &
  ) RESULT(strain_rate)
    USE tfm_density_tools, ONLY : tfm_density_arrhenius
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_NabarroHerringCreep
    !
    ! Maeno, N. and Ebinuma, T. Pressure Sintering of Ice and Its
    ! Implication to the Densification of Snow at Polar Glaciers and Ice
    ! Sheets. J. Phys. Chem., 87, 4103-4110, (1983).
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   temperature: Temperaure (K).
    !   grain_radius: Grain radius (m).
    !   driving force: Sintering driving force (Pa).
    !
    ! Result:
    !   strain_rate: Resulting strain rate (s**-1).
    !---------------------------------------------------------------------------

    REAL(dp), INTENT(IN) :: &
      temperature,          &
      grain_radius,         &
      driving_force

    REAL(dp) :: strain_rate
    REAL(dp) :: rate_factor_lattice_diffusion

    !---------------------------------------------------------------------------

    rate_factor_lattice_diffusion = tfm_density_arrhenius(    &
    &  1,                                                     &
    &  factor=PRE_FACTOR_LATTICE_DIFFUSION,                   &
    &  activation_energy=ACTIVATION_ENERGY_LATTICE_DIFFUSION, &
    &  temperature=temperature                                &
    )

    strain_rate = (                                                   &
    &  -((10.0_dp * VOLUME_H2O) / (3.0_dp * BOLTZMANN * temperature)) &
    &  * (1.0_dp / (grain_radius**2.0_dp))                            &
    &  * rate_factor_lattice_diffusion                                &
    &  * driving_force                                                &
    )
  END FUNCTION tfm_density_NabarroHerringCreep


  FUNCTION tfm_density_NabarroHerringCreepMod( &
  &  temperature,                              &
  &  density,                                  &
  &  grain_radius,                             &
  &  pore_radius,                              &
  &  driving_force                             &
  ) RESULT(strain_rate)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_NabarroHerringCreepMod
    !
    ! Maeno, N. and Ebinuma, T. Pressure Sintering of Ice and Its
    ! Implication to the Densification of Snow at Polar Glaciers and Ice
    ! Sheets. J. Phys. Chem., 87, 4103-4110, (1983).
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   temperature: Temperature (K).
    !   density: Density (kg m**-3).
    !   grain_radius: Grain radius (m).
    !   pore_radius: Pore radius (m).
    !   driving_force: Sintering dirving force (Pa).
    !
    ! Result:
    !   strain_rate: Resulting strain rate (s**-1).
    !---------------------------------------------------------------------------

    REAL(dp), INTENT(IN) :: &
      temperature,          &
      density,              &
      grain_radius,         &
      pore_radius,          &
      driving_force

    REAL(dp) :: strain_rate
    REAL(dp) :: rate_factor_lattice_diffusion

    !---------------------------------------------------------------------------

    rate_factor_lattice_diffusion = tfm_density_arrhenius(    &
    &  1,                                                     &
    &  factor=PRE_FACTOR_LATTICE_DIFFUSION,                   &
    &  activation_energy=ACTIVATION_ENERGY_LATTICE_DIFFUSION, &
    &  temperature=temperature                                &
    )

    strain_rate = (                                                  &
    &  -((3.0_dp * VOLUME_H2O) / (BOLTZMANN * temperature))          &
    &  * (1.0_dp / (grain_radius**2.0_dp))                           &
    &  * (ICE_DENSITY / density)                                     &
    &  * ((1.0_dp / grain_radius) * (                                &
    &    (pore_radius * grain_radius) / (grain_radius - pore_radius) &
    &  ))                                                            &
    &  * rate_factor_lattice_diffusion                               &
    &  * driving_force                                               &
    )
  END FUNCTION tfm_density_NabarroHerringCreepMod


  FUNCTION tfm_density_CobleCreep( &
  &  temperature,                  &
  &  grain_radius,                 &
  &  driving_force                 &
  ) RESULT(strain_rate)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_CobleCreep
    !
    ! Maeno, N. and Ebinuma, T. Pressure Sintering of Ice and Its
    ! Implication to the Densification of Snow at Polar Glaciers and Ice
    ! Sheets. J. Phys. Chem., 87, 4103-4110, (1983).
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   temperature: Temeprature (K).
    !   grain_radius: Grain raidus (m).
    !   driving_force: Sintering driving force (Pa).
    !
    ! Result:
    !   strain_rate: Resulting strain rate (s**-1).
    !---------------------------------------------------------------------------

    REAL(dp), INTENT(IN) :: &
      temperature,          &
      grain_radius,         &
      driving_force

    REAL(dp) :: strain_rate
    REAL(dp) :: rate_factor_boundary_diffusion

    !---------------------------------------------------------------------------

    rate_factor_boundary_diffusion = tfm_density_arrhenius(    &
    &  1,                                                      &
    &  factor=PRE_FACTOR_BOUNDARY_DIFFUSION,                   &
    &  activation_energy=ACTIVATION_ENERGY_BOUNDARY_DIFFUSION, &
    &  temperature=temperature                                 &
    )

    strain_rate = (                                       &
    &  -(                                                 &
    &    (37.0_dp * VOLUME_H2O * 2.0_dp * BURGERS_VECTOR) &
    &    / (2.0_dp * BOLTZMANN * temperature)             &
    &  )                                                  &
    &  * (1.0_dp / (grain_radius**3.0_dp))                &
    &  * rate_factor_boundary_diffusion                   &
    &  * driving_force                                    &
    )
  END FUNCTION tfm_density_CobleCreep


  FUNCTION tfm_density_DislocationCreep( &
  &  temperature,                        &
  &  density,                            &
  &  driving_force                       &
  ) RESULT(strain_rate)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_DislocationCreep
    !
    ! Maeno, N. and Ebinuma, T. Pressure Sintering of Ice and Its
    ! Implication to the Densification of Snow at Polar Glaciers and Ice
    ! Sheets. J. Phys. Chem., 87, 4103-4110, (1983).
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   temperature: Temeprature (K).
    !   density: Density (kg m**-3).
    !   driving_force: Sintering driving force (Pa).
    !
    ! Result:
    !   strain_rate: Resulting strain rate (s**-1).
    !---------------------------------------------------------------------------

    REAL(dp), INTENT(IN) :: &
      temperature,          &
      density,              &
      driving_force

    REAL(dp) :: strain_rate
    REAL(dp) :: rate_factor_dislocation_creep

    !---------------------------------------------------------------------------

    rate_factor_dislocation_creep = tfm_density_arrhenius( &
    &  1,                                                  &
    &  factor=3.22E-11_dp,                                 &
    &  activation_energy=74500.0_dp,                       &
    &  temperature=temperature                             &
    )

    strain_rate = (                                                                   &
    &  (-2.0_dp * rate_factor_dislocation_creep)                                      &
    &  * (                                                                            &
    &    (1.0_dp - (density / ICE_DENSITY))                                           &
    &    / ((1.0_dp - ((1.0_dp - (density / ICE_DENSITY))**(1.0_dp / ICE_N)))**ICE_N) &
    &  )                                                                              &
    &  * (((2.0_dp / ICE_N) * driving_force)**ICE_N)                                  &
    )
  END FUNCTION tfm_density_DislocationCreep


  FUNCTION tfm_density_DislocationCreepMod( &
  &  temperature,                           &
  &  density,                               &
  &  driving_force                          &
  ) RESULT(strain_rate)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_DislocationCreep
    !
    ! Maeno, N. and Ebinuma, T. Pressure Sintering of Ice and Its
    ! Implication to the Densification of Snow at Polar Glaciers and Ice
    ! Sheets. J. Phys. Chem., 87, 4103-4110, (1983).
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   temperature: Temeprature (K).
    !   density: Density (kg m**-3).
    !   driving_force: Sintering driving force (Pa).
    !
    ! Result:
    !   strain_rate: Resulting strain rate (s**-1).
    !---------------------------------------------------------------------------

    REAL(dp), INTENT(IN) :: &
      temperature,          &
      density,              &
      driving_force

    REAL(dp) :: strain_rate
    REAL(dp) :: rate_factor_dislocation_creep

    !---------------------------------------------------------------------------

    rate_factor_dislocation_creep = tfm_density_arrhenius( &
    &  1,                                                  &
    &  factor=3.22E-11_dp,                                 &
    &  activation_energy=74500.0_dp,                       &
    &  temperature=temperature                             &
    )

    strain_rate = (                                                                   &
    &  (-3.0_dp / 2.0_dp)                                                             &
    &  * rate_factor_dislocation_creep                                                &
    &  * (                                                                            &
    &    (1.0_dp - (density / ICE_DENSITY))                                           &
    &    / ((1.0_dp - ((1.0_dp - (density / ICE_DENSITY))**(1.0_dp / ICE_N)))**ICE_N) &
    &  )                                                                              &
    &  * (((3.0_dp / (2.0_dp * ICE_N)) * driving_force)**ICE_N)                       &
    )
  END FUNCTION tfm_density_DislocationCreepMod


  FUNCTION tfm_density_GrainBoundarySliding( &
  &  temperature,                            &
  &  density,                                &
  &  grain_radius,                           &
  &  neck_radius,                            &
  &  amplitude_gbo,                          &
  &  driving_force                           &
  ) RESULT(strain_rate)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_GrainBoundarySliding
    !
    ! Alley, R. B. Firn Densification by Grain Boundary Sliding: A First
    ! Model. J. Phys. Colloques, 48, C1, C1-249-C1-256, (1987), https://
    ! doi.org/10.1051/jphyscol:1987135
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   temperature: Temperature (K).
    !   density: Density (kg m**-3).
    !   grain_radius: Grain radius (m).
    !   neck_radius: Radius of necks between grains (m).
    !   amplitude_gbo: Amplitude of grain boundary obstructions (m).
    !   driving_force: Sintering driving force (Pa).
    !
    ! Reuslt:
    !   strain_rate: Resulting strain rate (s**-1).
    !---------------------------------------------------------------------------

    REAL(dp), INTENT(IN) :: &
      temperature,          &
      density,              &
      grain_radius,         &
      neck_radius,          &
      amplitude_gbo,        &
      driving_force

    REAL(dp) :: strain_rate
    REAL(dp) :: rate_factor_boundary_diffusion

    !---------------------------------------------------------------------------

    rate_factor_boundary_diffusion = tfm_density_arrhenius(    &
    &  1,                                                      &
    &  factor=PRE_FACTOR_BOUNDARY_DIFFUSION,                   &
    &  activation_energy=ACTIVATION_ENERGY_BOUNDARY_DIFFUSION, &
    &  temperature=temperature                                 &
    )

    strain_rate = (                                               &
    &  (-2.0_dp / 15.0_dp)                                        &
    &  * (2.0_dp * BURGERS_VECTOR)                                &
    &  * (                                                        &
    &    (8.0_dp * rate_factor_boundary_diffusion * VOLUME_H2O)   &
    &    / (BOLTZMANN * temperature * (amplitude_gbo**2.0_dp))    &
    &  )                                                          &
    &  * (grain_radius / (neck_radius**2.0_dp))                   &
    &  * ((ICE_DENSITY / density)**3.0_dp)                        &
    &  * (1.0_dp - ((5.0_dp * density) / (3.0_dp * ICE_DENSITY))) &
    &  * driving_force                                            &
    )
  END FUNCTION tfm_density_GrainBoundarySliding
END MODULE Tfm_density_processes



MODULE tfm_density_herronLangway
! ------------------------------------------------------------------------------
! Module: tfm_density_herronLangway
!
! This module contains generic routines to solve firn densification
! models of the Herron & Langway type. These routines can be used to
! solve various firn densification models based on this model type.
!
! Herron, M. M. and Langway, C. C. Firn Densification: An Empirical
! Model. Journal of Glaciology, 25 (93), 373-385, (1980). https://
! doi.org/10.3189/S0022143000015239
!
! Dependencies: tfm_essentials, tfm_constants
!
! Subroutines:
!   tfm_density_HLtype: Generic function for Herron & Langway type
!                       models.
! ------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! modules used
  !-----------------------------------------------------------------------------
  USE, INTRINSIC :: iso_fortran_env, ONLY : dp => real64

  USE tfm_essentials, ONLY : &
    tfm_essentials_mean_acc

  USE tfm_constants, ONLY : &
    ACC_GRAVITY,            &
    GAS_CONST,              &
    ICE_DENSITY,            &
    SECONDS_YEAR,           &
    WATER_DENSITY

  !-----------------------------------------------------------------------------
  ! declarations
  !-----------------------------------------------------------------------------
  IMPLICIT NONE (TYPE, EXTERNAL)
  PRIVATE

  !-----------------------------------------------------------------------------
  ! public routines
  !-----------------------------------------------------------------------------
  PUBLIC :: tfm_density_hltype

  !-----------------------------------------------------------------------------
  ! routines of this module
  !-----------------------------------------------------------------------------
  CONTAINS


  SUBROUTINE tfm_density_HLtype( &
  &  nz,                         &
  &  dt,                         &
  &  stage1_params,              &
  &  stage2_params,              &
  &  depth,                      &
  &  temperature,                &
  &  density,                    &
  &  age,                        &
  &  d_density                   &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Subroutine: tfm_density_HLtype
    !
    ! The routine computes the density change resulting from a "generic"
    ! densification model of the Herron & Langway type.
    !
    ! Herron, M. M. and Langway, C. C. Firn Densification: An Empirical
    ! Model. Journal of Glaciology, 25 (93), 373-385, (1980). https://
    ! doi.org/10.3189/S0022143000015239
    !
    ! Many firn densification models rely on this type of model. The idea
    ! is to formulate a generic function following the general form of the
    ! material model. It is then solved by passing three arguments for the
    ! first and the second stage of firn densification. All models following
    ! the concept of Herron & Langway can be broken down to these six
    ! parameters.
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension variables "depth", "temperature", "density", "age",
    !       "d_density".
    !   dt: Time step (s).
    !   stage1_params: Parameters for the first densification stage
    !                  (dimension (3)).
    !   stage2_params: Parameters for the seconds densification stage
    !                  (dimension (3)).
    !   depth: Depth of the firn profile (m).
    !   temperature: Temperature of the firn profile (K).
    !   density: Density of the firn profile (kg m**-3).
    !   age: Age of the firn profile (s).
    !   d_density - on input: Variable for the result.
    !
    ! Result:
    !   d_density - on output: Change of the density (kg m**-3).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN)  :: nz
    REAL(dp), INTENT(IN) :: dt

    REAL(dp), DIMENSION(3), INTENT(IN)  :: &
      stage1_params,                       &
      stage2_params

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      temperature,                         &
      density,                             &
      age

    REAL(dp), DIMENSION(nz), INTENT(INOUT) :: d_density

    REAL(dp), DIMENSION(nz) :: &
      mean_acc,                &
      c

    integer :: n

    !---------------------------------------------------------------------------

    ! computation of the mean accumulation rate
    ! over the lifetime of the firn parcel (kg a-1 m-2)
    CALL tfm_essentials_mean_acc( &
    &  nz=nz,                     &
    &  depth=depth,               &
    &  density=density,           &
    &  age=age,                   &
    &  mean_acc=mean_acc          &
    )
    mean_acc = (mean_acc * WATER_DENSITY)

    ! boundary between first and seconds stage
    DO n = nz, 1, -1

      ! first stage
      IF ( ( density(n) > 0.0_dp ) .AND. ( density(n) <= 550.0_dp ) ) THEN

        c(n) = (                                                   &
        &  stage1_params(1)                                        &
        &  * (mean_acc(n)**stage1_params(2))                       &
        &  * ACC_GRAVITY                                           &
        &  * exp(-stage1_params(3) / (GAS_CONST * temperature(n))) &
        )

      ! seconds stage
      ELSE IF ( ( density(n) > 550.0_dp ) .AND. ( density(n) <= ICE_DENSITY ) ) THEN

        c(n) = (                                                   &
        &  stage2_params(1)                                        &
        &  * (mean_acc(n)**stage2_params(2))                       &
        &  * ACC_GRAVITY                                           &
        &  * exp(-stage2_params(3) / (GAS_CONST * temperature(n))) &
        )

      ELSE

        PRINT *, 'module: tfm_density, subroutine: tfm_density_HLtype'
        PRINT *, 'The density seems to show an irregular value.'
        PRINT *, 'Something went wrong! Stopping right here!'
        STOP

      END IF
    END DO

    ! density change
    d_density = ((dt / SECONDS_YEAR) * (c * (ICE_DENSITY - density)))
  END SUBROUTINE tfm_density_HLtype
END MODULE TFM_density_herronLangway



MODULE tfm_density_fischmeisterArzt
! ---------------------------------------------------------------------
! Module: tfm_density_fischmeisterArzt
!
! The module contains routines to compute the coordination number and
! contact area of particles during sintering based on the density
! following the work of Fischmeister & Arzt. For different uses
! different functions handling different dimensions are connected via
! suitable interfaces.
!
! Fischmeister, H. F. and Arzt, E. Densification of powders by particle
! deformation. Powder Metallurgy, 26 (2), 82-88, (1983).
! https://doi.org/10.1179/pom.1983.26.2.82
!
! Dependencies: tfm_essentials, tfm_constants
!
! Interfaces:
!   tfm_density_arztCoordination: Interface for coordination number
!                                 ("tfm_density_arztCoordination_i",
!                                 "tfm_arztCoordination_n").
!   tfm_density_arztContactarea: Interface for contact area
!                                ("tfm_density_arztContactarea_i",
!                                "tfm_arztContactarea_n").
!
! Parameters:
!   ZZERO: Coordination number at densest packing (1).
!   CARZT: Factor (1)
!
! Functions:
!   tfm_density_arztCoordination_i: Coordination number, dimension (1).
!   tfm_density_arztCoordination_n: Coordination number, dimenison (n).
!   tfm_density_arztContactarea_i: Contact area, dimension(1).
!   tfm_density_arztContactarea_n: Contact area, dimension(n).
! ---------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! modules used
  !-----------------------------------------------------------------------------
  USE, INTRINSIC :: iso_fortran_env, ONLY : dp => real64

  USE tfm_constants, ONLY : PI

  !-----------------------------------------------------------------------------
  ! declarations
  !-----------------------------------------------------------------------------
  IMPLICIT NONE (TYPE, EXTERNAL)
  PRIVATE

  !-----------------------------------------------------------------------------
  ! public routines
  !-----------------------------------------------------------------------------
  PUBLIC ::                      &
    tfm_density_arztContactarea, &
    tfm_density_arztCoordination

  !-----------------------------------------------------------------------------
  ! interfaces
  !-----------------------------------------------------------------------------
  INTERFACE tfm_density_arztCoordination
    MODULE PROCEDURE                  &
      tfm_density_arztCoordination_i, &
      tfm_density_arztCoordination_n
  END INTERFACE tfm_density_arztCoordination

  INTERFACE tfm_density_arztContactarea
    MODULE PROCEDURE                 &
      tfm_density_arztContactarea_i, &
      tfm_density_arztContactarea_n
  END INTERFACE tfm_density_arztContactarea

  !-----------------------------------------------------------------------------
  ! parameters
  !-----------------------------------------------------------------------------
  REAL(dp), PARAMETER :: ZZERO = 7.3_dp
  !REAl(dp), parameter :: ZZERO = 4.81_dp
  REAL(dp), PARAMETER :: CARZT = 15.5_dp


  !-----------------------------------------------------------------------------
  ! routines of this module
  !-----------------------------------------------------------------------------
  CONTAINS


  FUNCTION tfm_density_arztCoordination_i( &
  &  n,                                    &
  &  rel_density,                          &
  &  d_zero                                &
  ) RESULT(coordination_number)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function tfm_density_arztCoordination_i
    !
    ! Computation of the coordination number following Fischmeister & Arzt
    ! (1983) for dimension (1).
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   n: Dummy argument for the interface (always = 1).
    !   rel_density: Relative density (1).
    !   d_zero: Relative density at densest packing (1).
    !
    ! Result:
    !   coordination_number: Coordination number at given density (1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: n

    REAL(dp), INTENT(IN) :: &
      rel_density,          &
      d_zero

    REAL(dp) :: coordination_number
    INTEGER  :: m

    !---------------------------------------------------------------------------

    m = 1 * n

    coordination_number = (                                             &
    &  ZZERO                                                            &
    &  + CARZT * (((rel_density / d_zero)**(1.0_dp / 3.0_dp)) - 1.0_dp) &
    )
  END FUNCTION tfm_density_arztCoordination_i


  FUNCTION tfm_density_arztCoordination_n( &
  &  n,                                    &
  &  rel_density,                          &
  &  d_zero                                &
  ) RESULT(coordination_number)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function tfm_density_arztCoordination_n
    !
    ! Computation of the coordination number following Fischmeister & Arzt
    ! (1983) for dimension (n).
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   n: Dimension of variables "rel_density", "coordination_number".
    !   rel_density: Relative density (1).
    !   d_zero: Relative density at densest packing (1).
    !
    ! Result:
    !   coordination_number: Coordination number at given density (1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN)                :: n
    REAL(dp), DIMENSION(n), INTENT(IN) :: rel_density
    REAL(dp), INTENT(IN)               :: d_zero
    REAL(dp), DIMENSION(n)             :: coordination_number

    !---------------------------------------------------------------------------

    coordination_number = (                                                   &
    &  ZZERO + CARZT * (((rel_density / d_zero)**(1.0_dp / 3.0_dp)) - 1.0_dp) &
    )
  END FUNCTION tfm_density_arztCoordination_n


  FUNCTION tfm_density_arztContactarea_i( &
  &  n,                                   &
  &  rel_density,                         &
  &  d_zero                               &
  ) RESULT(contactarea)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function tfm_density_arztContactarea_i
    !
    ! Computation of the contact area between grains  following
    ! Fischmeister & Arzt (1983) for dimension (1).
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   n: Dummy argument for the interface (always = 1).
    !   rel_density: Relative density (1).
    !   d_zero: Relative density at densest packing (1).
    !
    ! Result:
    !   contact_area: Contact area at given density (1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: n

    REAL(dp), INTENT(IN) :: &
      rel_density,          &
      d_zero

    REAL(dp) ::            &
      coordination_number, &
      r_i,                 &
      r_ii,                &
      contactarea

    INTEGER :: m

    !---------------------------------------------------------------------------

    m = 1 * n

    coordination_number = tfm_density_arztCoordination(n, rel_density, d_zero)

    r_i = (rel_density / d_zero)**(1.0_dp / 3.0_dp)
    r_ii = r_i + (                                                                 &
    &  (                                                                           &
    &    ((4.0_dp * ZZERO) * ((r_i - 1.0_dp)**2.0_dp) * ((2.0_dp * r_i) + 1.0_dp)) &
    &    + (CARZT * ((r_i - 1.0_dp)**3.0_dp) * ((3.0_dp * r_i) + 1.0_dp))          &
    &  )                                                                           &
    &  / (                                                                         &
    &    (12.0_dp * r_i)                                                           &
    &    * (                                                                       &
    &      (4.0_dp * r_i)                                                          &
    &      - (2.0_dp * ZZERO * (r_i - 1.0_dp))                                     &
    &      - (CARZT * ((r_i - 1.0_dp)**2.0_dp))                                    &
    &     )                                                                        &
    &  )                                                                           &
    )

    contactarea = (                                              &
    &  (PI / (3.0_dp * coordination_number * (r_i**2.0_dp)))     &
    &  * (                                                       &
    &    (3.0_dp * ((r_ii**2.0_dp) - 1.0_dp) * ZZERO)            &
    &    + ((r_ii**2.0_dp) * CARZT * ((2.0_dp * r_ii) - 3.0_dp)) &
    &    + (CARZT)                                               &
    &  )                                                         &
    )
  END FUNCTION tfm_density_arztContactarea_i


  FUNCTION tfm_density_arztContactarea_n( &
  &  n,                                   &
  &  rel_density,                         &
  &  d_zero                               &
  ) RESULT(contactarea)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function tfm_density_arztContactarea_n
    !
    ! Computation of the contact area between grains  following
    ! Fischmeister & Arzt (1983) for dimension (n).
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   n: Dimension of variables "rel_density", "coordination_number".
    !   rel_density: Relative density (1).
    !   d_zero: Relative density at densest packing (1).
    !
    ! Result:
    !   contact_area: Contact area at given density (1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN)                :: n
    REAL(dp), DIMENSION(n), INTENT(IN) :: rel_density
    REAL(dp), INTENT(IN)               :: d_zero

    REAL(dp), DIMENSION(n) :: &
      coordination_number,    &
      r_i,                    &
      r_ii,                   &
      contactarea

    !---------------------------------------------------------------------------

    coordination_number = tfm_density_arztCoordination( &
    &  n,                                               &
    &  rel_density=rel_density,                         &
    &  d_zero=d_zero                                    &
    )

    r_i = (rel_density / d_zero)**(1.0_dp / 3.0_dp)
    r_ii = r_i + (                                                                 &
    &  (                                                                           &
    &    ((4.0_dp * ZZERO) * ((r_i - 1.0_dp)**2.0_dp) * ((2.0_dp * r_i) + 1.0_dp)) &
    &    + (CARZT * ((r_i - 1.0_dp)**3.0_dp) * ((3.0_dp * r_i) + 1.0_dp))          &
    &  )                                                                           &
    &  / (                                                                         &
    &    (12.0_dp * r_i)                                                           &
    &    * (                                                                       &
    &      (4.0_dp * r_i)                                                          &
    &      - (2.0_dp * ZZERO * (r_i - 1.0_dp))                                     &
    &      - (CARZT * ((r_i - 1.0_dp)**2.0_dp)))                                   &
    &  )                                                                           &
    )

    contactarea = (                                              &
    &  (PI / (3.0_dp * coordination_number * (r_i**2.0_dp)))     &
    &  * (                                                       &
    &    (3.0_dp * ((r_ii**2.0_dp) - 1.0_dp) * ZZERO)            &
    &    + ((r_ii**2.0_dp) * CARZT * ((2.0_dp * r_ii) - 3.0_dp)) &
    &    + (CARZT)                                               &
    &  )                                                         &
    )
  END FUNCTION tfm_density_arztContactarea_n
END MODULE tfm_density_fischmeisterArzt


MODULE tfm_density_gagliardini
! ------------------------------------------------------------------------------
! Module: tfm_density_gagliardini
!
! The module contains routines to solve a firn densification model of
! the kind first described for firn by:
!
! Gagliardini, O. and Meyssonnier, J. Flow simulation of a firn-covered
! cold glacier. Annals of Glaciology, 24, 242-248 (1997).
! https://doi.org/10.3189/S0260305500012246
!
! The model is based on the work of Duva & Crow:
!
! Duva, J. M. and Crow, P. D. Analysis of consolidation of reinforced
! materials by power-law creep. Mechanics of Materials, 17 (1), 25-32,
! (1994). https://doi.org/10.1016/0167-6636(94)90011-6
!
! Dependencies: tfm_essentials, tfm_constants, tfm_denisty_tools
!
! Interfaces:
!   invariant_inter: Interface for the invariant function.
!   viscosity_inter: Interface for viscosity functions.
!
! Parameters:
!   STAGE_DIV: Relative density at transition stage II - stage III (1).
!
! Functions:
!   tfm_density_gagliardiniParamA0: Parameter A0 of the model.
!   tfm_density_gagliardiniParamB0: Parameter B0 of the model.
!   tfm_density_gagliardiniRate: Rate factor.
!   tfm_density_gagliradiniSolve: Function for solving the model
!                                 interatively.
! ------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! modules used
  !-----------------------------------------------------------------------------
  USE, INTRINSIC :: iso_fortran_env, ONLY : dp => real64

  USE tfm_constants, ONLY : &
    ICE_DENSITY,            &
    ICE_N

  USE tfm_density_tools, ONLY : &
    tfm_density_arrhenius

  !-----------------------------------------------------------------------------
  ! declarations
  !-----------------------------------------------------------------------------
  IMPLICIT NONE (TYPE, EXTERNAL)
  PRIVATE

  !-----------------------------------------------------------------------------
  ! public routines
  !-----------------------------------------------------------------------------
  PUBLIC ::                         &
    tfm_density_gagliardiniParamA0, &
    tfm_density_gagliardiniParamB0, &
    tfm_density_gagliardiniRate,    &
    tfm_density_gagliardiniSolve

  !-----------------------------------------------------------------------------
  ! interfaces
  !-----------------------------------------------------------------------------
  INTERFACE
    FUNCTION invariant_inter( &
    &  nz,                    &
    &  param_a,               &
    &  param_b,               &
    &  strain_rate_inp        &
    ) RESULT(invariant)
      USE, INTRINSIC :: iso_fortran_env, ONLY : dp => real64
      IMPLICIT NONE (TYPE, EXTERNAL)

      INTEGER, INTENT(IN)                 :: nz

      REAL(dp), DIMENSION(nz), INTENT(IN) :: &
        param_a,                             &
        param_b

      REAL(dp), DIMENSION(nz), INTENT(IN), OPTIONAL :: strain_rate_inp

      REAL(dp), DIMENSION(nz) :: invariant
    END FUNCTION invariant_inter
  END INTERFACE


  INTERFACE
    FUNCTION viscosity_inter( &
    &  nz,                    &
    &  param,                 &
    &  rate_factor,           &
    &  invariant              &
    ) RESULT(viscosity)
      USE, INTRINSIC :: iso_fortran_env, ONLY : dp => real64
      IMPLICIT NONE (TYPE, EXTERNAL)

      INTEGER, INTENT(IN)                 :: nz

      REAL(dp), DIMENSION(nz), INTENT(IN) :: &
        param,                               &
        rate_factor,                         &
        invariant

      REAL(dp), DIMENSION(nz) :: viscosity
    END FUNCTION viscosity_inter
  END INTERFACE

  !-----------------------------------------------------------------------------
  ! parameters
  !-----------------------------------------------------------------------------
  REAL(dp), PARAMETER :: STAGE_DIV = 0.81_dp

  !-----------------------------------------------------------------------------
  ! routines of this module
  !-----------------------------------------------------------------------------
  CONTAINS


  FUNCTION tfm_density_gagliardiniParamA0( &
  &  density                               &
  ) RESULT(param_a0)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_gagliardiniParmA0
    !
    ! Parameter a0 of the firn densification model as described by
    ! Gagliardini & Meyssonnier (1997), following Duva & Crow (1994).
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   density: Absolute density (kg m**-3).
    !
    ! Result:
    !   param_a0: Parameter a0.
    !---------------------------------------------------------------------------

    REAL(dp), INTENT(IN) :: density

    REAL(dp) :: &
      param_a0, &
      rel_density

    !---------------------------------------------------------------------------

    rel_density = (density / ICE_DENSITY)

    param_a0 = (                                               &
    &  (1.0_dp + ((2.0_dp / 3.0_dp) * (1.0_dp - rel_density))) &
    &  * (rel_density)**((-2.0_dp * ICE_N) / (ICE_N + 1.0_dp)) &
    )
  END FUNCTION tfm_density_gagliardiniParamA0


  FUNCTION tfm_density_gagliardiniParamB0( &
  &  density                               &
  ) RESULT(param_b0)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_gagliardiniParamB0
    !
    ! Parameter b0 of the firn densification model as described by
    ! Gagliardini & Meyssonnier (1997), following Duva & Crow (1994).
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   density: Absolute density (kg m**-3).
    !
    ! Result:
    !   param_b0: Parameter b0.
    !---------------------------------------------------------------------------

    REAL(dp), INTENT(IN) :: density

    REAL(dp) :: &
      param_b0, &
      rel_density

    !---------------------------------------------------------------------------

    rel_density = (density / ICE_DENSITY)

    param_b0 = (                                                           &
    &  (3.0_dp / 4.0_dp)                                                   &
    &  * (                                                                 &
    &    ((1.0_dp - rel_density)**(1.0_dp / ICE_N))                        &
    &    / (ICE_N * (1.0_dp - ((1.0_dp - rel_density)**(1.0_dp / ICE_N)))) &
    &  )**((2.0_dp * ICE_N) / (ICE_N + 1.0_dp))                            &
    )
  END FUNCTION tfm_density_gagliardiniParamB0


  FUNCTION tfm_density_gagliardiniRate( &
  &  nz,                                &
  &  temperature                        &
  ) RESULT(rate_factor)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_gagliardiniRate
    !
    ! Rate factor for ice flow following:
    ! Greve, R. and Blatter, H. Dynamics of Ice Sheets and Galciers.
    ! Springer, Berlin, (2009).
    !
    ! Author: Timm Schultz
    !
    ! Parameters (Greve & Blatter, 2009):
    !   TEMP_DIV: Temperature dividing high and low temperature (K).
    !   PRE_FACTOR_LOW: Pre factor of the arrhenius equation at low
    !                   temperatures.
    !   PRE_FACTOR_HIGH: Pre factor of the arrhenius equation at high
    !                    temperatures.
    !   ACTIVATION_ENERGY_LOW: Activation energy of the arrhenius equation
    !                          at low temperatures.
    !   ACTIVATION_ENERGY_HIGH: Activation energy of the arrhenius equation
    !                           at hight temperatures.
    !
    ! Arguments:
    !   nz: Dimension of variables "temperature" and "rate_factor".
    !   temperature: Temperature (K).
    !
    ! Result:
    !   rate_factor: Rate factor.
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), DIMENSION(nz), INTENT(IN) :: temperature

    INTEGER :: n
    REAL(dp), DIMENSION(nz) :: rate_factor

    ! parameters (values from Greve & Blatter, 2009)
    REAL(dp), PARAMETER ::                  &
      TEMP_DIV = 263.15_dp,                 &
      PRE_FACTOR_LOW  = 3.985e-13_dp,       &
      PRE_FACTOR_HIGH = 1.916e3_dp,         &
      ACTIVATION_ENERGY_LOW  =  60000.0_dp, &
      ACTIVATION_ENERGY_HIGH = 139000.0_dp

    !---------------------------------------------------------------------------

    DO n = 1, nz, 1
      IF ( temperature(n) <= TEMP_DIV ) THEN

        rate_factor(n) = tfm_density_arrhenius( &
        &  1,                                   &
        &  PRE_FACTOR_LOW,                      &
        &  ACTIVATION_ENERGY_LOW,               &
        &  temperature(n)                       &
        &)

      ELSE IF ( temperature(n) > TEMP_DIV ) THEN

        rate_factor(n) = tfm_density_arrhenius( &
        &  1,                                   &
        &  PRE_FACTOR_HIGH,                     &
        &  ACTIVATION_ENERGY_HIGH,              &
        &  temperature(n)                       &
        &)

      ELSE

        ! catch exception
        PRINT *, 'module: tfm_density_gagliardini'
        PRINT *, 'function: tfm_density_gagliardiniRate'
        PRINT *, ''
        PRINT *, 'It seems there are irregular temperature values!'
        PRINT *, ''
        PRINT *, 'Stopping right here!'
        STOP
      END IF
    END DO

    rate_factor = rate_factor**(-1.0_dp / ICE_N)
  END FUNCTION tfm_density_gagliardiniRate


  FUNCTION tfm_density_gagliardiniSolve( &
  &  nz,                                 &
  &  density,                            &
  &  stress,                             &
  &  dt,                                 &
  &  param_a,                            &
  &  param_b,                            &
  &  rate_factor,                        &
  &  invariant_func,                     &
  &  shear_visco_func,                   &
  &  bulk_visco_func                     &
  ) RESULT(d_density)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function tfm_density_gagliardiniSolve
    !
    ! This function solves a firn densification model of the kind described
    ! by Gagliardini & Meyssonnier (1997). Because this kind of model is
    ! non-linear it is solved iteratively.
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "density", "stress", "param_a",
    !       "param_b", "rate_factor".
    !   density: Absolute density (kg m**-3).
    !   stress: Stress due to overburden firn (Pa).
    !   dt: Time step (s).
    !   param_a: Parameter a along the firn profile.
    !   param_b: Parameter b along the firn profile.
    !   rate_factor: Ice flow rate factor.
    !   invariant_func: Function defining the strain invariant (procedure,
    !                   interface: invariant_func).
    !   shear_visco_func: Function defining the shear viscosity (procedure,
    !                     interface: viscosity_inter).
    !   bulk_visco_func: Function defining the bulk viscosity (procedure,
    !                    interface: viscosity_inter).
    !
    ! Result:
    !   d_density: Density change along the firn profile (kg m**-3).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      density,                             &
      stress,                              &
      param_a,                             &
      param_b,                             &
      rate_factor

    PROCEDURE(invariant_inter) :: invariant_func
    PROCEDURE(viscosity_inter) :: &
      shear_visco_func,           &
      bulk_visco_func

    REAL(dp), DIMENSION(nz) :: d_density

    INTEGER :: &
      n,       &
      iter

    REAL(dp), DIMENSION(nz) :: &
      strain_rate,             &
      d_density_prev,          &
      invariant,               &
      bulk_viscosity,          &
      shear_viscosity

    INTEGER, PARAMETER :: MAX_ITER = 1000

    !---------------------------------------------------------------------------

    ! first guess of the invariant
    invariant = invariant_func(nz, param_a, param_b)

    ! viscosities
    shear_viscosity = shear_visco_func(    &
    &  nz, param_a, rate_factor, invariant &
    )
    bulk_viscosity = bulk_visco_func(      &
    &  nz, param_b, rate_factor, invariant &
    )

    ! strain rate
    strain_rate = (                                                       &
    &  (1.0_dp/ (((4.0_dp / 3.0_dp) * shear_viscosity) + bulk_viscosity)) &
    &  * stress                                                           &
    )

    d_density      = -999999.9
    d_density_prev = +999999.9
    iter           = 0

    DO WHILE ( (maxval(abs(d_density - d_density_prev)) > 1.0e-2) .OR. (iter == MAX_ITER) )

      ! first guess of the invariant
      invariant = invariant_func(          &
      &  nz, param_a, param_b, strain_rate &
      )

      ! viscosities
      shear_viscosity = shear_visco_func(    &
      &  nz, param_a, rate_factor, invariant &
      )
      bulk_viscosity = bulk_visco_func(      &
      &  nz, param_b, rate_factor, invariant &
      )

      ! strain rate
      strain_rate = (                                                       &
      &  (1.0_dp/ (((4.0_dp / 3.0_dp) * shear_viscosity) + bulk_viscosity)) &
      &  * stress                                                           &
      )

      ! densification
      d_density_prev = (1.0_dp * d_density)
      d_density = (dt * strain_rate * density)

      iter = (iter + 1)
    END DO

    ! avoid the singularity at ice density
    DO n = 1, nz, 1
      IF ( (density(n) + d_density(n)) > (ICE_DENSITY - 10.0D-5) ) THEN
        d_density(n) = (ICE_DENSITY - density(n)) - 10.0D-5
      END IF
    END DO
  END FUNCTION tfm_density_gagliardiniSolve
END MODULE tfm_density_gagliardini



MODULE tfm_density
! ------------------------------------------------------------------------------
! Module: tfm_density
!
! This module is a collection of various different firn densification
! models. All functions compute the change in density. The arguments
! passed to the individual functions are always the same allowing to
! pass them via a suitable interface. This means that some arguments
! may not be used by the function.
!
! Dependencies: tfm_essentials, tfm_constants, tfm_density_tools,
!               tfm_density_herronLangway, tfm_density_stress,
!               tfm_density_fischmeisterArzt, tfm_density_processes,
!               tfm_density_gagliardini
!
! Functions:
!   tfm_density_depth: Change in depth due to densification.
!   tfm_density_gagliardini1998: Gagliardini & Meyssonier (1998).
!   tfm_density_timmsfit: Timm's version of the Gagliardini model.
!   tfm_density_greve2009: Greve & Blatter (2009).
!   tfm_density_zwinger2007: Zwinger et al. (2007).
!   tfm_density_breant2017: Breant et al. (2017).
!   tfm_density_medley2020: Medley et al. (2020) (preprint).
!   tfm_density_herron1980: Herron & Langway (1980) (transient).
!   tfm_density_arthern1998: Arthern & Wingham (1998).
!   tfm_density_li2003: Li & Zwally (2003).
!   tfm_density_helsen2008: Helsen et al. (2008).
!   tfm_density_arthern2010: Arthern et al. (2010).
!   tfm_density_ligtenberg2011: Ligtenberg et al. (2011).
!   tfm_density_simonsen2013: Simonsen et al. (2013).
!
! Subroutines:
!   tfm_density_herron1980_analytical: Herron & Langway (analytical).
! ------------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! modules used
  !-----------------------------------------------------------------------------
  USE, INTRINSIC :: iso_fortran_env, ONLY : dp => real64

  USE tfm_essentials, ONLY :   &
    tfm_essentials_do_nothing, &
    tfm_essentials_mean_acc

  USE tfm_constants, ONLY : &
    ICE_DENSITY,            &
    ICE_N,                  &
    PI,                     &
    GAS_CONST,              &
    CLOSEOFF_DENSITY,       &
    TEMP_OFFSET,            &
    CRITICAL_DENSITY,       &
    WATER_DENSITY,          &
    ACC_GRAVITY,            &
    SECONDS_YEAR

  USE tfm_density_tools, ONLY : &
    tfm_density_bagmean,        &
    tfm_density_lin_interp

  USE tfm_density_herronLangway, ONLY : &
    tfm_density_hltype

  USE tfm_density_stress, ONLY :   &
    tfm_density_rel_boylemariotte, &
    tfm_density_computeStress

  USE tfm_density_fischmeisterArzt, ONLY : &
    tfm_density_arztContactarea,           &
    tfm_density_arztCoordination

  USE tfm_density_processes, ONLY :     &
    tfm_density_CobleCreep,             &
    tfm_density_NabarroHerringCreepMod, &
    tfm_density_DislocationCreep,       &
    tfm_density_NabarroHerringCreep,    &
    tfm_density_DislocationCreepMod,    &
    tfm_density_GrainBoundarySliding

  USE tfm_density_gagliardini, ONLY : &
    tfm_density_gagliardiniSolve,     &
    tfm_density_gagliardiniParamA0,   &
    tfm_density_gagliardiniParamB0,   &
    tfm_density_gagliardiniRate

  !-----------------------------------------------------------------------------
  ! declarations
  !-----------------------------------------------------------------------------
  IMPLICIT NONE (TYPE, EXTERNAL)
  PRIVATE

  !-----------------------------------------------------------------------------
  ! public routines
  !-----------------------------------------------------------------------------
  PUBLIC ::                      &
    tfm_density_arthern1998,     &
    tfm_density_li2003,          &
    tfm_density_helsen2008,      &
    tfm_density_breant2017,      &
    tfm_density_zwinger2007,     &
    tfm_density_greve2009,       &
    tfm_density_gagliardini1998, &
    tfm_density_simonsen2013,    &
    tfm_density_medley2020,      &
    tfm_density_herron1980,      &
    tfm_density_arthern2010,     &
    tfm_density_ligtenberg2011,  &
    tfm_density_timmsfit,        &
    tfm_density_sintering,       &
    tfm_density_depth

  !-----------------------------------------------------------------------------
  ! routines of this module
  !-----------------------------------------------------------------------------
  CONTAINS


  FUNCTION tfm_density_depth( &
  &  nz,                      &
  &  depth,                   &
  &  density,                 &
  &  d_density                &
  ) RESULT(d_depth)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_depth
    !
    ! The function computes the change in depth along a firn profile from
    ! change in density.
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "depth", "density", "d_density",
    !       "d_depth".
    !   depth: Depth along the firn profile (m).
    !   density: Density along the firn profile (kg m**-3).
    !   d_density: Density change along the firn profile (kg m**-3).
    !
    ! Result:
    !   d_depth: Change in depth along the firn profile (m).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      density,                             &
      d_density

    INTEGER :: n
    REAL(dp), DIMENSION(nz) :: &
      d_depth,                 &
      dz,                      &
      ddz

    !---------------------------------------------------------------------------

    dz(1) = 0.0_dp
    dz(2:nz) = (depth(2:nz) - depth(1:nz-1))
    ddz = (dz * (density / (density + d_density))) - dz

    d_depth(1) = 0.0_dp
    DO n = 2, nz, 1
      d_depth(n) = d_depth(n-1) + ddz(n)
    END DO
  END FUNCTION tfm_density_depth


  FUNCTION tfm_density_gagliardini1998( &
  &  nz,                                &
  &  dt,                                &
  &  depth,                             &
  &  density,                           &
  &  temperature,                       &
  &  age,                               &
  &  grain_radius                       &
  ) RESULT(d_density)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_denisty_gagliardini1998
    !
    ! Gagliardini, O. and Meyssonnier, J. Flow simulation of a firn-covered
    ! cold glacier. Annals of Glaciology, 24, 242-248, (1997).
    ! https://doi.org/10.3189/S0260305500012246
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "depth", "density", "temperature", "age",
    !       "grain_radius", "d_density".
    !   dt: Time step (s).
    !   depth: Depth along the firn profile (m).
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temeperature along the firn profile (K).
    !   age: Age along the firn profile (s).
    !   grain_radius: Grain radius along the firn profile (m).
    !
    ! Result:
    !   d_density: Density change along the firn profile (kg m**-3).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt
    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      density,                             &
      temperature,                         &
      age,                                 &
      grain_radius

    INTEGER :: n
    REAL(dp) :: rel_density
    REAL(dp), DIMENSION(nz) :: &
      d_density,               &
      param_a,                 &
      param_b,                 &
      rate_factor,             &
      stress

    !---------------------------------------------------------------------------

    ! doing nothing
    CALL tfm_essentials_do_nothing(nz, age)
    CALL tfm_essentials_do_nothing(nz, grain_radius)

    ! stress
    CALL tfm_density_computeStress( &
    &  nz=nz,                       &
    &  depth=depth,                 &
    &  density=density,             &
    &  stress=stress                &
    )

    ! density dependent parameter of the model as defined by Greve & BLatter 2009
    DO n = 1, nz, 1

      rel_density = density(n) / ICE_DENSITY

      IF ( (rel_density > 0.0_dp) .AND. (rel_density < 0.785_dp) ) THEN
        param_a(n) = exp((-19.67_dp * rel_density) + 15.94_dp)
        param_b(n) = exp((-27.65_dp * rel_density) + 20.37_dp)

      ELSE IF ( (rel_density >= 0.785_dp) .and. (rel_density < 1.0_dp) ) THEN
        param_a(n) = tfm_density_gagliardiniParamA0(density(n))
        param_b(n) = tfm_density_gagliardiniParamB0(density(n))

      ELSE
        ! catch exception
        PRINT *, 'module: tfm_density'
        PRINT *, 'function: tfm_density_gagliardini1998'
        PRINT *, ''
        PRINT *, 'It seems the density exceeds the range of valid '
        PRINT *, 'values at some point!'
        PRINT *, ''
        PRINT *, 'Stopping right here!'
        STOP

      END IF
    END DO

    ! temperature dependent associated rate factor
    rate_factor = tfm_density_gagliardiniRate(nz, temperature)

    ! solving for the density change
    d_density = tfm_density_gagliardiniSolve( &
    &  nz=nz,                                 &
    &  density=density,                       &
    &  stress=stress,                         &
    &  dt=dt,                                 &
    &  param_a=param_a,                       &
    &  param_b=param_b,                       &
    &  rate_factor=rate_factor,               &
    &  invariant_func=invariant_func,         &
    &  shear_visco_func=visco_func,           &
    &  bulk_visco_func=visco_func             &
    )

    CONTAINS

    FUNCTION invariant_func( &
    &  nz,                   &
    &  param_a,              &
    &  param_b,              &
    &  strain_rate_inp       &
    ) RESULT(invariant)
      IMPLICIT NONE (TYPE, EXTERNAL)

      INTEGER, INTENT(IN)                 :: nz
      REAL(dp), DIMENSION(nz), INTENT(IN) :: &
        param_a,                             &
        param_b
      REAL(dp), DIMENSION(nz), INTENT(IN), OPTIONAL :: &
        strain_rate_inp

      REAL(dp), DIMENSION(nz) :: &
        invariant,               &
        strain_rate

      IF ( present(strain_rate_inp) ) THEN
        strain_rate = strain_rate_inp
      ELSE
        strain_rate = 1.0e-10_dp
      END IF

      invariant = (                                                       &
      &  strain_rate                                                      &
      &  * (((3.0_dp / (4.0_dp * param_a)) + (1.0_dp / param_b))**0.5_dp) &
      )
    END FUNCTION invariant_func


    FUNCTION visco_func( &
    &  nz,               &
    &  param,            &
    &  rate_factor,      &
    &  invariant         &
    ) RESULT(viscosity)
      IMPLICIT NONE (TYPE, EXTERNAL)

      INTEGER, INTENT(IN) :: nz
      REAL(dp), DIMENSION(nz), INTENT(IN) :: &
        param,                               &
        rate_factor,                         &
        invariant

      REAL(dp), DIMENSION(nz) :: viscosity

      viscosity = (                                        &
      &  (1.0_dp / param)                                &
      &  * rate_factor                                     &
      &  * (invariant**(-(1.0_dp - (1.0_dp / ICE_N)))) &
      )
    END FUNCTION visco_func
  END FUNCTION tfm_density_gagliardini1998


  FUNCTION tfm_density_timmsfit( &
  &  nz,                         &
  &  dt,                         &
  &  depth,                      &
  &  density,                    &
  &  temperature,                &
  &  age,                        &
  &  grain_radius                &
  ) RESULT(d_density)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_timmsfit
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "depth", "density", "temperature", "age",
    !       "grain_radius", "d_density".
    !   dt: Time step (s).
    !   depth: Depth along the firn profile (m).
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temeperature along the firn profile (K).
    !   age: Age along the firn profile (s).
    !   grain_radius: Grain radius along the firn profile (m).
    !
    ! Result:
    !   d_density: Density change along the firn profile (kg m**-3).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt
    REAL(dp), DIMENSION(nz), intent(in) :: &
      depth,                               &
      density,                             &
      temperature,                         &
      age,                                 &
      grain_radius

    INTEGER :: n
    REAL(dp) :: rel_density
    REAL(dp), DIMENSION(nz) :: &
      d_density,               &
      param_a,                 &
      param_b,                 &
      rate_factor,             &
      stress

    !---------------------------------------------------------------------------

    ! doing nothing
    CALL tfm_essentials_do_nothing(nz, age)
    CALL tfm_essentials_do_nothing(nz, grain_radius)

    ! stress
    CALL tfm_density_computeStress( &
    &  nz=nz,                       &
    &  depth=depth,                 &
    &  density=density,             &
    &  stress=stress                &
    )

    ! density dependent parameter of the model as defined by Greve & BLatter 2009
    DO n = 1, nz, 1
      rel_density = density(n) / ICE_DENSITY

      IF ( (rel_density > 0.0_dp) .AND. (rel_density <= 0.79_dp) ) THEN
        param_a(n) = exp(                       &
        &  24.60215_dp                          &
        &  - (58.573530_dp * rel_density)       &
        &  - (-35.5_dp * (rel_density**2.0_dp)) &
        )
        param_b(n) = (                                    &
        &  (                                              &
        &    tfm_density_gagliardiniParamB0(density(n))   &
        &    / tfm_density_gagliardiniParamA0(density(n)) &
        &  ) * param_a(n)                                 &
        )

      ELSE IF ( (rel_density > 0.79_dp) .AND. (rel_density < 1.0_dp) ) THEN
        param_a(n) = tfm_density_gagliardiniParamA0(density(n))
        param_b(n) = tfm_density_gagliardiniParamB0(density(n))

      ELSE
        ! catch exception
        PRINT *, 'module: tfm_density'
        PRINT *, 'function: tfm_density_timmsfit'
        PRINT *, ''
        PRINT *, 'It seems the density exceeds the range of valid '
        PRINT *, 'values at some point!'
        PRINT *, ''
        PRINT *, 'Stopping right here!'
        STOP
      END IF
    END DO

    ! temperature dependent associated rate factor
    rate_factor = tfm_density_gagliardiniRate( &
    &  nz=nz,                                  &
    &  temperature=temperature                 &
    )

    ! solving for the density change
    d_density = tfm_density_gagliardiniSolve( &
    &  nz=nz,                                 &
    &  density=density,                       &
    &  stress=stress,                         &
    &  dt=dt,                                 &
    &  param_a=param_a,                       &
    &  param_b=param_b,                       &
    &  rate_factor=rate_factor,               &
    &  invariant_func=invariant_func,         &
    &  shear_visco_func=visco_func,           &
    &  bulk_visco_func=visco_func             &
    )

    CONTAINS

    FUNCTION invariant_func( &
    &  nz,                   &
    &  param_a,              &
    &  param_b,              &
    &  strain_rate_inp       &
    ) RESULT(invariant)
      IMPLICIT NONE (TYPE, EXTERNAL)

      INTEGER, INTENT(IN)                 :: nz

      REAL(dp), DIMENSION(nz), INTENT(IN) :: &
        param_a,                             &
        param_b

      real(dp), dimension(nz), intent(in), optional :: &
        strain_rate_inp

      REAL(dp), DIMENSION(nz) :: &
        invariant,               &
        strain_rate

      IF ( present(strain_rate_inp) ) THEN
        strain_rate = strain_rate_inp
      ELSE
        strain_rate = 1.0e-10_dp
      END IF

      invariant = (                                                     &
      &  strain_rate                                                    &
      &  * ((                                                           &
      &    (1.0_dp / (3.0_dp * param_a)) + (1.0_dp / (4.0_dp * param_b) &
      &  ))**0.5_dp)                                                    &
      )
    END FUNCTION invariant_func


    FUNCTION visco_func( &
    &  nz,               &
    &  param,            &
    &  rate_factor,      &
    &  invariant         &
    ) RESULT(viscosity)

      INTEGER, INTENT(IN)                 :: nz
      REAL(dp), DIMENSION(nz), INTENT(IN) :: &
        param,                               &
        rate_factor,                         &
        invariant

      REAL(dp), DIMENSION(nz) :: viscosity

      viscosity = (                                    &
      &  (1.0_dp / (2.0 * param))                      &
      &  * rate_factor                                 &
      &  * (invariant**(-(1.0_dp - (1.0_dp / ICE_N)))) &
      )
    END FUNCTION visco_func
  END FUNCTION tfm_density_timmsfit


  FUNCTION tfm_density_greve2009( &
  &  nz,                          &
  &  dt,                          &
  &  depth,                       &
  &  density,                     &
  &  temperature,                 &
  &  age,                         &
  &  grain_radius                 &
  ) RESULT(d_density)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_greve2009
    !
    ! Greve, R. and Blatter, H. Dynamics of Ice Sheets and Glaciers.
    ! Springer, Berlin, (2009).
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "depth", "density", "temperature", "age",
    !       "grain_radius", "d_density".
    !   dt: Time step (s).
    !   depth: Depth along the firn profile (m).
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temeperature along the firn profile (K).
    !   age: Age along the firn profile (s).
    !   grain_radius: Grain radius along the firn profile (m).
    !
    ! Result:
    !   d_density: Density change along the firn profile (kg m**-3).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt
    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      density,                             &
      temperature,                         &
      age,                                 &
      grain_radius

    INTEGER :: n
    REAL(dp) :: rel_density
    REAL(dp), DIMENSION(nz) :: &
      d_density,               &
      param_a,                 &
      param_b,                 &
      rate_factor,             &
      stress

    !---------------------------------------------------------------------------

    ! doing nothing
    CALL tfm_essentials_do_nothing(nz, age)
    CALL tfm_essentials_do_nothing(nz, grain_radius)

    ! stress
    CALL tfm_density_computeStress( &
    &  nz=nz,                       &
    &  depth=depth,                 &
    &  density=density,             &
    &  stress=stress                &
    )

    ! density dependent parameter of the model as defined by Greve & BLatter 2009
    DO n = 1, nz, 1
      rel_density = density(n) / ICE_DENSITY

      IF ( (rel_density > 0.0_dp) .AND. (rel_density <= 0.81_dp) ) THEN
        param_a(n) = exp(13.22240_dp - (15.78652_dp * rel_density))
        param_b(n) = exp(15.09371_dp - (20.46489_dp * rel_density))

      ELSE IF ( (rel_density > 0.81_dp) .AND. (rel_density < 1.0_dp) ) THEN
        param_a(n) = tfm_density_gagliardiniParamA0(density(n))
        param_b(n) = tfm_density_gagliardiniParamB0(density(n))

      ELSE
        ! catch exception
        PRINT *, 'module: tfm_density'
        PRINT *, 'function: tfm_density_greve2009'
        PRINT *, ''
        PRINT *, 'It seems the density exceeds the range of valid '
        PRINT *, 'values at some point!'
        PRINT *, ''
        PRINT *, 'Stopping right here!'
        STOP
      END IF
    END DO

    ! temperature dependent associated rate factor
    rate_factor = tfm_density_gagliardiniRate(nz, temperature)

    ! solving for the density change
    d_density = tfm_density_gagliardiniSolve( &
    &  nz=nz,                                 &
    &  density=density,                       &
    &  stress=stress,                         &
    &  dt=dt,                                 &
    &  param_a=param_a,                       &
    &  param_b=param_b,                       &
    &  rate_factor=rate_factor,               &
    &  invariant_func=invariant_func,         &
    &  shear_visco_func=visco_func,           &
    &  bulk_visco_func=visco_func             &
    )

    CONTAINS

    FUNCTION invariant_func( &
    &  nz,                   &
    &  param_a,              &
    &  param_b,              &
    &  strain_rate_inp       &
    ) RESULT(invariant)
      IMPLICIT NONE (TYPE, EXTERNAL)

      INTEGER, INTENT(IN)                 :: nz
      REAL(dp), DIMENSION(nz), INTENT(IN) :: &
        param_a,                             &
        param_b
      REAL(dp), DIMENSION(nz), INTENT(IN), OPTIONAL :: &
        strain_rate_inp

      REAL(dp), DIMENSION(nz) :: &
        invariant,               &
        strain_rate

      IF ( present(strain_rate_inp) ) THEN
        strain_rate = strain_rate_inp
      ELSE
        strain_rate = 1.0D-10
      END IF

      invariant = (                                                      &
      &  strain_rate                                                     &
      &  * ((                                                            &
      &    (1.0_dp / (3.0_dp * param_a)) + (1.0_dp / (4.0_dp * param_b)) &
      &  )**0.5_dp)                                                      &
      )
    END FUNCTION invariant_func


    FUNCTION visco_func( &
    &  nz,               &
    &  param,            &
    &  rate_factor,      &
    &  invariant         &
    ) RESULT(viscosity)

      INTEGER, INTENT(IN)                 :: nz
      REAL(dp), DIMENSION(nz), INTENT(IN) :: &
        param,                               &
        rate_factor,                         &
        invariant

      REAL(dp), DIMENSION(nz) :: viscosity

      viscosity = (                                    &
      &  (1.0_dp / (2.0_dp * param))                   &
      &  * rate_factor                                 &
      &  * (invariant**(-(1.0_dp - (1.0_dp / ICE_N)))) &
      )
    END FUNCTION visco_func
  END FUNCTION tfm_density_greve2009


  FUNCTION tfm_density_zwinger2007( &
  &  nz,                            &
  &  dt,                            &
  &  depth,                         &
  &  density,                       &
  &  temperature,                   &
  &  age,                           &
  &  grain_radius                   &
  ) RESULT(d_density)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_zwinger2007
    !
    ! Zwinger, T., Greve, R., Gagliardini, O., Shiraiwa, T., and Lyly, M. A
    ! full Stokes-flow thermo-mechanical model for firn and ice applied to
    ! the Gorshkov crater glacier, Kamchatka. Annals of Glaciology, 45,
    ! 29-37, (2007). https://doi.org/10.3189/172756407782282543
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "depth", "density", "temperature", "age",
    !       "grain_radius", "d_density".
    !   dt: Time step (s).
    !   depth: Depth along the firn profile (m).
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temeperature along the firn profile (K).
    !   age: Age along the firn profile (s).
    !   grain_radius: Grain radius along the firn profile (m).
    !
    ! Result:
    !   d_density: Density change along the firn profile (kg m**-3).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt
    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      density,                             &
      temperature,                         &
      age,                                 &
      grain_radius

    INTEGER :: n
    REAL(dp) :: rel_density
    REAL(dp), DIMENSION(nz) :: &
      d_density,               &
      param_a,                 &
      param_b,                 &
      rate_factor,             &
      stress

    !---------------------------------------------------------------------------

    ! doing nothing
    CALL tfm_essentials_do_nothing(nz, age)
    CALL tfm_essentials_do_nothing(nz, grain_radius)

    ! stress
    CALL tfm_density_computeStress( &
    &  nz=nz,                       &
    &  depth=depth,                 &
    &  density=density,             &
    &  stress=stress                &
    )

    ! density dependent parameter of the model as defined by Greve & BLatter 2009
    DO n = 1, nz, 1
      rel_density = density(n) / ICE_DENSITY

      IF ( (rel_density > 0.0_dp) .AND. (rel_density <= 0.81_dp) ) THEN
        param_a(n) = exp(13.22240_dp - (15.78652_dp * rel_density))
        param_b(n) = exp(15.09371_dp - (20.46489_dp * rel_density))

      ELSE IF ( (rel_density > 0.81_dp) .AND. (rel_density < 1.0_dp) ) THEN
        param_a(n) = tfm_density_gagliardiniParamA0(density(n))
        param_b(n) = tfm_density_gagliardiniParamB0(density(n))

      ELSE
        ! catch exception
        PRINT *, 'module: tfm_density'
        PRINT *, 'function: tfm_density_zwinger2007'
        PRINT *, ''
        PRINT *, 'It seems the density exceeds the range of valid '
        PRINT *, 'values at some point!'
        PRINT *, ''
        PRINT *, 'Stopping right here!'
        STOP
      END IF
    END DO

    ! temperature dependent associated rate factor
    rate_factor = tfm_density_gagliardiniRate( &
    &  nz=nz,                                  &
    &  temperature=temperature                 &
    )

    ! solving for the density change
    d_density = tfm_density_gagliardiniSolve( &
    &  nz=nz,                                 &
    &  density=density,                       &
    &  stress=stress,                         &
    &  dt=dt,                                 &
    &  param_a=param_a,                       &
    &  param_b=param_b,                       &
    &  rate_factor=rate_factor,               &
    &  invariant_func=invariant_func,         &
    &  shear_visco_func=shear_visco_func,     &
    &  bulk_visco_func=bulk_visco_func        &
    )

    CONTAINS

    FUNCTION invariant_func( &
    &  nz,                   &
    &  param_a,              &
    &  param_b,              &
    &  strain_rate_inp       &
    ) RESULT(invariant)
      IMPLICIT NONE (TYPE, EXTERNAL)

      INTEGER, INTENT(IN)                 :: nz
      REAL(dp), DIMENSION(nz), INTENT(IN) :: &
        param_a,                             &
        param_b
      REAL(dp), DIMENSION(nz), INTENT(IN), OPTIONAL :: &
        strain_rate_inp

      REAL(dp), DIMENSION(nz) :: &
        invariant,               &
        strain_rate

      IF ( present(strain_rate_inp) ) THEN
        strain_rate = strain_rate_inp
      ELSE
        strain_rate = 1.0D-10
      END IF

      invariant = (                                                       &
      &  strain_rate                                                      &
      &  * (((3.0_dp / (4.0_dp * param_a)) + (1.0_dp / param_b))**0.5_dp) &
      )
    END FUNCTION invariant_func


    FUNCTION shear_visco_func( &
    &  nz,                     &
    &  param_a,                &
    &  rate_factor,            &
    &  invariant               &
    ) RESULT(shear_viscosity)

      INTEGER, INTENT(IN)                 :: nz
      REAL(dp), DIMENSION(nz), INTENT(IN) :: &
        param_a,                             &
        rate_factor,                         &
        invariant

      REAL(dp), DIMENSION(nz) :: shear_viscosity

      shear_viscosity = (                              &
      &  (2.0_dp / param_a)                            &
      &  * rate_factor                                 &
      &  * (invariant**(-(1.0_dp - (1.0_dp / ICE_N)))) &
      )
    END FUNCTION shear_visco_func


    FUNCTION bulk_visco_func( &
    &  nz,                    &
    &  param_b,               &
    &  rate_factor,           &
    &  invariant              &
    ) RESULT(bulk_viscosity)

      INTEGER, INTENT(IN)                 :: nz
      REAL(dp), DIMENSION(nz), INTENT(IN) :: &
        param_b,                             &
        rate_factor,                         &
        invariant

      REAL(dp), DIMENSION(nz) :: bulk_viscosity

      bulk_viscosity = (                               &
      &  (1.0_dp / param_b)                            &
      &  * rate_factor                                 &
      &  * (invariant**(-(1.0_dp - (1.0_dp / ICE_N)))) &
      )
    END FUNCTION bulk_visco_func
  END FUNCTION tfm_density_zwinger2007


  FUNCTION tfm_density_breant2017( &
  &  nz,                           &
  &  dt,                           &
  &  depth,                        &
  &  density,                      &
  &  temperature,                  &
  &  age,                          &
  &  grain_radius                  &
  ) RESULT(d_density)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_breant2017
    !
    ! Breant, C., Martinerie, P., Orsi, A., Arnaud, L., and Landais, A.
    ! Modelling firn thickness evolution during the last deglaciation:
    ! constraints on sensitivity to temperature and impurities. Clim. Past,
    ! 13, 833-853, (2017). https://doi.org/10.5194/cp-13-833-2017
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "depth", "density", "temperature", "age",
    !       "grain_radius", "d_density".
    !   dt: Time step (s).
    !   depth: Depth along the firn profile (m).
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temeperature along the firn profile (K).
    !   age: Age along the firn profile (s).
    !   grain_radius: Grain radius along the firn profile (m).
    !
    ! Result:
    !   d_density: Density change along the firn profile (kg m**-3).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt
    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      density,                             &
      temperature,                         &
      age,                                 &
      grain_radius

    integer :: &
      n,       &
      m,       &
      mz

    REAL(dp), DIMENSION(nz) :: &
      rel_density,             &
      d_density,               &
      arrhenius_dc,            &
      stress,                  &
      eff_stress,              &
      coord_number,            &
      contact_area,            &
      strain_rate

    REAL(dp) ::           &
      gamma_gbs,          &
      dz,                 &
      arrhenius_gamma,    &
      contact_area_gamma, &
      eff_stress_gamma,   &
      stress_gamma,       &
      temperature_gamma

    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: bagmean

    REAL(dp), PARAMETER :: &
      DZERO = 0.56_dp,     &
      MIN_STRESS = 0.1D5,  &
      A0 = 7.89D-15,       &
      A1 = 1.05D9,         &
      A2 = 1400.0_dp,      &
      A3 = 6.0D-15,        &
      Q1 = 110000.0_dp,    &
      Q2 =  75000.0_dp,    &
      Q3 =   1500.0_dp,    &
      QGBS = 49500.0_dp

    !---------------------------------------------------------------------------

    CALL tfm_essentials_do_nothing(nz, age)
    CALL tfm_essentials_do_nothing(nz, grain_radius)

    rel_density = (density / ICE_DENSITY)

    ! temperature dependence
    arrhenius_dc = A0 * (                            &
    &    (A1 * exp(-Q1 / (GAS_CONST * temperature))) &
    &  + (A2 * exp(-Q2 / (GAS_CONST * temperature))) &
    &  + (A3 * exp(-Q3 / (GAS_CONST * temperature))) &
    )

    ! model by Arzt
    coord_number = tfm_density_arztCoordination( &
    &  n=nz,                                     &
    &  rel_density=rel_density,                  &
    &  d_zero=DZERO                              &
    )
    contact_area = tfm_density_arztContactarea( &
    &  n=nz,                                    &
    &  rel_density=rel_density,                 &
    &  d_zero=DZERO                             &
    )

    ! stress
    CALL tfm_density_computeStress( &
    &  nz=nz,                       &
    &  depth=depth,                 &
    &  density=density,             &
    &  stress=stress                &
    )
    DO n = 1, nz, 1
      IF (stress(n) < MIN_STRESS) stress(n) = MIN_STRESS
    END DO
    eff_stress = stress * (                                        &
    &  (4.0_dp * PI) / (contact_area * coord_number * rel_density) &
    )


    dz = 1.0_dp
    mz = floor((depth(nz) - depth(1)) / dz)
    ALLOCATE(bagmean(2,mz))
    CALL tfm_density_bagmean( &
    &  nz=nz,                 &
    &  depth=depth,           &
    &  density=rel_density,   &
    &  dz=dz,                 &
    &  mz=mz,                 &
    &  bagmean=bagmean        &
    )

    DO m = mz, 1, -1
      IF ( bagmean(2,m) > 0.6_dp ) EXIT
    END DO

    DO n = nz, 1, -1
      IF ( depth(n) < bagmean(1,m+1) ) EXIT
    END DO

    CALL tfm_density_lin_interp(         &
    &  z0=depth(n+1),                    &
    &  z1=depth(n),                      &
    &  v0=stress(n+1),                   &
    &  v1=stress(n),                     &
    &  dz=(bagmean(1,m+1) - depth(n+1)), &
    &  v=stress_gamma                    &
    )
    CALL tfm_density_lin_interp(         &
    &  z0=depth(n+1),                    &
    &  z1=depth(n),                      &
    &  v0=temperature(n+1),              &
    &  v1=temperature(n),                &
    &  dz=(bagmean(1,m+1) - depth(n+1)), &
    &  v=temperature_gamma               &
    )

    do n = nz, 1, -1
      IF ( depth(n) < bagmean(1,m) ) EXIT
    end do

    CALL tfm_density_lin_interp(       &
    &  z0=depth(n+1),                  &
    &  z1=depth(n),                    &
    &  v0=arrhenius_dc(n+1),           &
    &  v1=arrhenius_dc(n),             &
    &  dz=(bagmean(1,m) - depth(n+1)), &
    &  v=arrhenius_gamma               &
    )
    CALL tfm_density_lin_interp(       &
    &  z0=depth(n+1),                  &
    &  z1=depth(n),                    &
    &  v0=contact_area(n+1),           &
    &  v1=contact_area(n),             &
    &  dz=(bagmean(1,m) - depth(n+1)), &
    &  v=contact_area_gamma            &
    )
    CALL tfm_density_lin_interp(       &
    &  z0=depth(n+1),                  &
    &  z1=depth(n),                    &
    &  v0=eff_stress(n+1),             &
    &  v1=eff_stress(n),               &
    &  dz=(bagmean(1,m) - depth(n+1)), &
    &  v=eff_stress_gamma              &
    )

    gamma_gbs = (                                              &
    &  (5.3_dp * arrhenius_gamma)                              &
    &  * ((DZERO * (bagmean(2,m)**2.0_dp))**(1.0_dp / 3.0_dp)) &
    &  * ((contact_area_gamma / PI)**0.5_dp)                   &
    &  * ((eff_stress_gamma / 3.0_dp)**ICE_N)                  &
    )
    gamma_gbs = gamma_gbs / (                                                &
    &  (stress_gamma / (bagmean(2,m+1)**2.0_dp))                             &
    &  * (1.0_dp + (0.5_dp / 6.0_dp) - ((5.0_dp / 3.0_dp) * bagmean(2,m+1))) &
    &  * (exp(-QGBS / (GAS_CONST * temperature_gamma)))                      &
    )


    DO n = 1, nz, 1

      IF ( rel_density(n) < 0.6_dp ) THEN

        strain_rate(n) = (                                                       &
        &  gamma_gbs                                                             &
        &  * (stress(n) / (rel_density(n)**2.0_dp))                              &
        &  * (1.0_dp + (0.5_dp / 6.0_dp) - ((5.0_dp / 3.0_dp) * rel_density(n))) &
        &  * (exp(-QGBS / (GAS_CONST * temperature(n))))                         &
        )

      ELSE IF ( ( rel_density(n) >= 0.6_dp ) .AND. ( rel_density(n) < 1.0_dp) ) THEN

        IF ( rel_density(n) > 0.9_dp ) THEN
          eff_stress(n) = tfm_density_rel_boylemariotte(rel_density(n))
        END IF

        strain_rate(n) = (                                           &
        &  (5.3_dp * arrhenius_dc(n))                                &
        &  * ((DZERO * (rel_density(n)**2.0_dp))**(1.0_dp / 3.0_dp)) &
        &  * ((contact_area(n) / PI)**0.5_dp)                        &
        &  * ((eff_stress(n) / 3.0_dp)**ICE_N)                       &
        )

      ELSE

        ! catch exception
        PRINT *, 'module: tfm_density, function tfm_density_brean2017'
        PRINT *, 'The density exceeds ice density!'
        PRINT *, 'Stopping right here!'
        STOP

      END IF
    END DO

    ! densification
    d_density = (dt * strain_rate * density)
  END FUNCTION tfm_density_breant2017


  FUNCTION tfm_density_medley2020( &
  &  nz,                           &
  &  dt,                           &
  &  depth,                        &
  &  density,                      &
  &  temperature,                  &
  &  age,                          &
  &  grain_radius                  &
  ) RESULT(d_density)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function tfm_density_medley2020
    !
    ! Medley, B., Neumann, T. A., Zwally, H. J., and Smith, B. E. Forty-year
    ! Simulations of Firn Processes over the Greenland and Antarctic Ice
    ! Sheets. The Cryosphere Discuss. [preprint], in review, (2020).
    ! https://doi.org/10.5194/tc-2020-266
    !
    ! NOTE! This functions is based on the preprint version of the paper.
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "depth", "density", "temperature", "age",
    !       "grain_radius", "d_density".
    !   dt: Time step (s).
    !   depth: Depth along the firn profile (m).
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temeperature along the firn profile (K).
    !   age: Age along the firn profile (s).
    !   grain_radius: Grain radius along the firn profile (m).
    !
    ! Result:
    !   d_density: Density change along the firn profile (kg m**-3).
    !---------------------------------------------------------------------------

    REAL(dp), PARAMETER :: &
      ! first stage parameters
      ALPHA0 = 0.9250_dp,  &
      EC0    = 60000.0_dp, &
      ! seconds stage parameters
      ALPHA1 = 0.6354_dp,  &
      EC1    = 56973.0_dp, &
      ! other parameters
      EG = 42400.0_dp

    REAL(dp) :: &
      a0,       &
      a1

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt
    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      density,                             &
      temperature,                         &
      age,                                 &
      grain_radius

    REAL(dp), DIMENSION(nz) :: d_density
    INTEGER :: n
    REAL(dp) :: mean_temperature

    !---------------------------------------------------------------------------

    CALL tfm_essentials_do_nothing(nz, grain_radius)

    ! 10 m temperature
    DO n = nz, 1, -1
      IF ( (depth(n) - depth(nz)) <= -10.0_dp ) EXIT
    END DO
    mean_temperature = temperature(n)

    ! factor depending on the mean annual temperature
    a0 = (0.07_dp * exp(EG / (GAS_CONST * mean_temperature)))
    a1 = (0.03_dp * exp(EG / (GAS_CONST * mean_temperature)))

    ! call Herron & Langway model
    CALL tfm_density_HLtype(              &
    &  nZ=nz,                             &
    &  dt=dt,                             &
    &  stage1_params=[ a0, ALPHA0, EC0 ], &
    &  stage2_params=[ a1, ALPHA1, EC1 ], &
    &  depth=depth,                       &
    &  temperature=temperature,           &
    &  density=density,                   &
    &  age=age,                           &
    &  d_density=d_density                &
    )
  END FUNCTION tfm_density_medley2020


  FUNCTION tfm_density_herron1980( &
  &  nz,                           &
  &  dt,                           &
  &  depth,                        &
  &  density,                      &
  &  temperature,                  &
  &  age,                          &
  &  grain_radius                  &
  ) RESULT(d_density)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_herron1980
    !
    ! Herron, M. M. and Langway, C. C. Firn Densification: An Empirical
    ! Model. Journal of Glaciology, 25 (93), 373-385, (1980).
    ! https://doi.org/10.3189/S0022143000015239
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "depth", "density", "temperature", "age",
    !       "grain_radius", "d_density".
    !   dt: Time step (s).
    !   depth: Depth along the firn profile (m).
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temeperature along the firn profile (K).
    !   age: Age along the firn profile (s).
    !   grain_radius: Grain radius along the firn profile (m).
    !
    ! Result:
    !   d_density: Density change along the firn profile (kg m**-3).
    !---------------------------------------------------------------------------

    REAL(dp), PARAMETER ::                               &
      ! first stage parameters
      ALPHA0 = 1.0_dp,                                   &
      EC0 = 10160.0_dp,                                  &
      A0 = (11.0_dp * (0.001_dp**ALPHA0) / ACC_GRAVITY), &
      ! second stage parameters
      ALPHA1 = 0.5_dp,                                   &
      EC1 = 21400.0_dp,                                  &
      A1 = (575.0_dp * (0.001_dp**ALPHA1) / ACC_GRAVITY)

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt
    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      density,                             &
      temperature,                         &
      age,                                 &
      grain_radius

    REAL(dp), DIMENSION(nz) :: d_density

    !---------------------------------------------------------------------------

    CALL tfm_essentials_do_nothing(nz, grain_radius)

    ! call Herron & Langway model
    CALL tfm_density_HLtype(              &
    &  nz=nz,                             &
    &  dt=dt,                             &
    &  stage1_params=[ A0, ALPHA0, EC0 ], &
    &  stage2_params=[ A1, ALPHA1, EC1 ], &
    &  depth=depth,                       &
    &  temperature=temperature,           &
    &  density=density,                   &
    &  age=age,                           &
    &  d_density=d_density                &
    )
  END FUNCTION tfm_density_herron1980


  FUNCTION tfm_density_sintering( &
  &  nz,                          &
  &  dt,                          &
  &  depth,                       &
  &  density,                     &
  &  temperature,                 &
  &  age,                         &
  &  grain_radius                 &
  ) RESULT(d_density)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt
    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      density,                             &
      temperature,                         &
      age,                                 &
      grain_radius

    REAL(dp), DIMENSION(nz) :: d_density
    INTEGER  :: n

    REAL(dp) ::         &
      neck_radius,      &
      pore_radius,      &
      driving_force,    &
      alt_grain_radius, &
      contact_area

    REAL(dp), DIMENSION(nz) :: &
      stress,                  &
      strain_rate

    REAL(dp), PARAMETER ::        &
      AMPLITUDE_GBO = 4.25E-6_dp, & ! m
      MU = 0.7_dp,                &
      AIR_PRESSURE = 101325.0_dp ! Pa

    !---------------------------------------------------------------------------

    CALL tfm_essentials_do_nothing(nz, age)

    ! computation of the current stress state
    CALL tfm_density_computeStress( &
    &  nz=nz,                       &
    &  depth=depth,                 &
    &  density=density,             &
    &  stress=stress                &
    )

    strain_rate(:) = 0.0_dp
    DO n = 1, nz, 1

      alt_grain_radius = (                            &
      &  ((density(n) / 360.0_dp)**(1.0_dp / 3.0_dp)) &
      &  * 0.0005_dp                                  &
      )
      contact_area = tfm_density_arztContactarea( &
      &  1,                                       &
      &  rel_density=(density(n) / ICE_DENSITY),  &
      &  d_zero=(359.0_dp / ICE_DENSITY)          &
      )
      contact_area = contact_area * (alt_grain_radius**2.0_dp)
      neck_radius = ((contact_area**0.5_dp) / PI)


      IF ( density(n) <= CRITICAL_DENSITY ) THEN

        !neck_radius = (MU * grain_radius(n))
        !neck_radius = (MU * alt_grain_radius)

        strain_rate(n) = strain_rate(n) + tfm_density_GrainBoundarySliding( &
        &  temperature=temperature(n),                                      &
        &  density=density(n),                                              &
        &  grain_radius=alt_grain_radius,                                   &
        &  neck_radius=neck_radius,                                         &
        &  amplitude_gbo=AMPLITUDE_GBO,                                     &
        &  driving_force=stress(n)                                          &
        )
      END IF


      IF ( density(n) <= CLOSEOFF_DENSITY ) THEN

        driving_force = ((ICE_DENSITY / density(n)) * stress(n))

        strain_rate(n) = strain_rate(n) + tfm_density_NabarroHerringCreep( &
        &  temperature=temperature(n),                                     &
        &  grain_radius=alt_grain_radius,                                  &
        &  driving_force=driving_force                                     &
        )
        strain_rate(n) = strain_rate(n) + tfm_density_CobleCreep( &
        &  temperature=temperature(n),                            &
        &  grain_radius=alt_grain_radius,                         &
        &  driving_force=driving_force                            &
        )

        driving_force = stress(n)

        strain_rate(n) = strain_rate(n) + tfm_density_DislocationCreep( &
        &  temperature=temperature(n),                                  &
        &  density=density(n),                                          &
        &  driving_force=driving_force                                  &
        )

      ELSE IF ( (density(n) > CLOSEOFF_DENSITY) .AND. (density(n) < ICE_DENSITY) ) THEN

        pore_radius = (                                               &
        &  ((1.0_dp - (density(n) / ICE_DENSITY))**(1.0_dp / 3.0_dp)) &
        &  * alt_grain_radius                                         &
        )

        !driving_force = ((ICE_DENSITY / density(n)) * stress(n))
        driving_force = (                                          &
        &    (ICE_DENSITY / density(n))                            &
        &    * (                                                   &
        &      stress(n)                                           &
        &      - (AIR_PRESSURE * (                                 &
        &        (density(n) * (ICE_DENSITY - CLOSEOFF_DENSITY))   &
        &        / (CLOSEOFF_DENSITY * (ICE_DENSITY - density(n))) &
        &      ))                                                  &
        &    )                                                     &
        )

        strain_rate(n) = strain_rate(n) + tfm_density_NabarroHerringCreepMod( &
        &  temperature=temperature(n),                                        &
        &  density=density(n),                                                &
        &  grain_radius=alt_grain_radius,                                     &
        &  pore_radius=pore_radius,                                           &
        &  driving_force=driving_force                                        &
        )
        strain_rate(n) = strain_rate(n) + tfm_density_CobleCreep( &
        &  temperature=temperature(n),                            &
        &  grain_radius=alt_grain_radius,                         &
        &  driving_force=driving_force                            &
        )

        !driving_force = stress(n)
        driving_force = (                                      &
        &  stress(n)                                           &
        &  - (AIR_PRESSURE * (                                 &
        &    (density(n) * (ICE_DENSITY - CLOSEOFF_DENSITY))   &
        &    / (CLOSEOFF_DENSITY * (ICE_DENSITY - density(n))) &
        &  ))                                                  &
        )

        strain_rate(n) = strain_rate(n) + tfm_density_DislocationCreepMod( &
        &  temperature=temperature(n),                                     &
        &  density=density(n),                                             &
        &  driving_force=driving_force                                     &
        )

      ELSE IF ( density(n) >= ICE_DENSITY ) THEN

        strain_rate(n) = strain_rate(n) + 0.0_dp

      ELSE
        PRINT *, 'Module: tfm_density'
        PRINT *, 'Function: tfm_density_sintering'
        PRINT *, ''
        PRINT *, 'The density seems to show and irregular value.'
        PRINT *, 'Stopping right here!'
        STOP
      END IF
    END DO

    ! density change
    d_density = -(dt * strain_rate * density)
  END FUNCTION tfm_density_sintering


  FUNCTION tfm_density_arthern1998( &
  &  nz,                            &
  &  dt,                            &
  &  depth,                         &
  &  density,                       &
  &  temperature,                   &
  &  age,                           &
  &  grain_radius                   &
  ) RESULT(d_density)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_arthern1998
    !
    ! Arthern, R. J. and Wingham, D. J. The Natural Fluctuations of Firn
    ! Densification and Their Effect on the Geodetic Determination of Ice
    ! Sheet Mass Balance. Climatic Change, 40, 605-624, (1998).
    ! https://doi.org/10.1023/A:1005320713306
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "depth", "density", "temperature", "age",
    !       "grain_radius", "d_density".
    !   dt: Time step (s).
    !   depth: Depth along the firn profile (m).
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temeperature along the firn profile (K).
    !   age: Age along the firn profile (s).
    !   grain_radius: Grain radius along the firn profile (m).
    !
    ! Parameters:
    !   AMPLITUDE_GBO: Amplitude of grain boundary obstructions (m), see
    !     Arthern & Wingham (1998).
    !   MU: Ratio of neck radius to grain radius.
    !   AIR_PRESSURE: Normal air pressure for the calculation of the
    !     Boyle-Mariotte law in the third stage of densification (Pa).
    !
    ! Result:
    !   d_density: Density change along the firn profile (kg m**-3).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt
    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      density,                             &
      temperature,                         &
      age,                                 &
      grain_radius

    INTEGER :: n
    REAL(dp), DIMENSION(nz) :: &
      d_density,               &
      stress,                  &
      strain_rate

    REAL(dp) ::      &
      driving_force, &
      neck_radius,   &
      pore_radius

    REAL(dp), PARAMETER ::       &
      AMPLITUDE_GBO = 4.0E-6_dp, & ! m
      MU = 0.7_dp,               &
      AIR_PRESSURE = 101325.0 ! (Pa)

    !---------------------------------------------------------------------------

    CALL tfm_essentials_do_nothing(nz, age)

    ! computation of the current stress state
    CALL tfm_density_computeStress( &
    &  nz=nz,                       &
    &  depth=depth,                 &
    &  density=density,             &
    &  stress=stress                &
    )

    ! computation of the strain rate
    strain_rate(:) = 0.0_dp

    DO n = 1, nz, 1

      IF ( (density(n) > 0.0_dp) .AND. (density(n) <= CRITICAL_DENSITY) ) THEN

        neck_radius = (MU * grain_radius(n))

        strain_rate(n) = strain_rate(n) + tfm_density_GrainBoundarySliding( &
        &  temperature=temperature(n),                                      &
        &  density=density(n),                                              &
        &  grain_radius=grain_radius(n),                                    &
        &  neck_radius=neck_radius,                                         &
        &  amplitude_gbo=AMPLITUDE_GBO,                                     &
        &  driving_force=stress(n)                                          &
        )

      ELSE IF ( (density(n) > CRITICAL_DENSITY) .AND. (density(n) <= CLOSEOFF_DENSITY) ) THEN

        driving_force = ((ICE_DENSITY / density(n)) * stress(n))

        strain_rate(n) = strain_rate(n) + tfm_density_NabarroHerringCreep( &
        &  temperature=temperature(n),                                     &
        &  grain_radius=grain_radius(n),                                   &
        &  driving_force=driving_force                                     &
        )

        strain_rate(n) = strain_rate(n) + tfm_density_CobleCreep( &
        &  temperature=temperature(n),                            &
        &  grain_radius=grain_radius(n),                          &
        &  driving_force=driving_force                            &
        )

        strain_rate(n) = strain_rate(n) + tfm_density_DislocationCreep( &
        &  temperature=temperature(n),                                  &
        &  density=density(n),                                          &
        &  driving_force=stress(n)                                      &
        )

      ELSE IF ( (density(n) > CLOSEOFF_DENSITY) .AND. (density(n) < ICE_DENSITY)) THEN

        pore_radius = (                                               &
        &  ((1.0_dp - (density(n) / ICE_DENSITY))**(1.0_dp / 3.0_dp)) &
        &  * grain_radius(n)                                          &
        )

        driving_force = (                                          &
        &    (ICE_DENSITY / density(n))                            &
        &    * (                                                   &
        &      stress(n)                                           &
        &      - (AIR_PRESSURE * (                                 &
        &        (density(n) * (ICE_DENSITY - CLOSEOFF_DENSITY))   &
        &        / (CLOSEOFF_DENSITY * (ICE_DENSITY - density(n))) &
        &      ))                                                  &
        &    )                                                     &
        )

        strain_rate(n) = strain_rate(n) + tfm_density_NabarroHerringCreepMod( &
        &  temperature=temperature(n),                                        &
        &  density=density(n),                                                &
        &  grain_radius=grain_radius(n),                                      &
        &  pore_radius=pore_radius,                                           &
        &  driving_force=driving_force                                        &
        )

        strain_rate(n) = strain_rate(n) + tfm_density_CobleCreep( &
        &  temperature=temperature(n),                            &
        &  grain_radius=grain_radius(n),                          &
        &  driving_force=driving_force                            &
        )

        driving_force = (                                      &
        &  stress(n)                                           &
        &  - (AIR_PRESSURE * (                                 &
        &    (density(n) * (ICE_DENSITY - CLOSEOFF_DENSITY))   &
        &    / (CLOSEOFF_DENSITY * (ICE_DENSITY - density(n))) &
        &  ))                                                  &
        )

        strain_rate(n) = strain_rate(n) + tfm_density_DislocationCreepMod( &
        &  temperature=temperature(n),                                     &
        &  density=density(n),                                             &
        &  driving_force=driving_force                                     &
        )

      ELSE IF ( density(n) >= ICE_DENSITY ) THEN

        strain_rate(n) = (strain_rate(n) + 0.0_dp)

      ELSE

        PRINT *, 'Module: tfm_density'
        PRINT *, 'Function: tfm_density_arthern1998'
        PRINT *, ''
        PRINT *, 'The density seems to show an irregular value.'
        PRINT *, 'Stopping right here!'
        STOP

      END IF
    END DO

    ! densification from strain rate
    d_density = (-dt * strain_rate * density)
  END FUNCTION tfm_density_arthern1998


  FUNCTION tfm_density_li2003( &
  &  nz,                       &
  &  dt,                       &
  &  depth,                    &
  &  density,                  &
  &  temperature,              &
  &  age,                      &
  &  grain_radius              &
  ) RESULT(d_density)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_li2003
    !
    ! Li, J., Zwally, H. J., Corneja, H., and Yi, D. Seasonal variation of
    ! snow-surface elevation in North Greenland as modeled and detected by
    ! satellite radar altimetry. Annals of Glaciology, 37, 223-238, (2003).
    ! https://doi.org/10.3189/172756403781815889
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "depth", "density", "temperature", "age",
    !       "grain_radius", "d_density".
    !   dt: Time step (s).
    !   depth: Depth along the firn profile (m).
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temeperature along the firn profile (K).
    !   age: Age along the firn profile (s).
    !   grain_radius: Grain radius along the firn profile (m).
    !
    ! Result:
    !   d_density: Density change along the firn profile (kg m**-3).
    !---------------------------------------------------------------------------

    REAL(dp), parameter :: &
      ! first stage parameters
      ALPHA0 = 1.0_dp,     &
      ! second stage parameters
      ALPHA1 = 1.0_dp

    REAL(dp) :: &
      ec0,      &
      a0,       &
      ec1,      &
      a1

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt
    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      density,                             &
      temperature,                         &
      age,                                 &
      grain_radius

    INTEGER :: n
    REAL(dp), DIMENSION(nz) :: d_density
    REAL(dp)                :: &
      mean_temperature,        &
      grain_growth_rate,       &
      beta

    !---------------------------------------------------------------------------

    CALL tfm_essentials_do_nothing(nz, grain_radius)

    ! 10 m temperature
    DO n = nz, 1, -1
      IF ( (depth(n) - depth(nz)) <= -10.0_dp ) EXIT
    END DO
    mean_temperature = temperature(n)

    ! grain growth rate and beta factor
    grain_growth_rate = (                                   &
    &  8.36_dp                                              &
    &  * (abs(mean_temperature - TEMP_OFFSET)**(-2.061_dp)) &
    )
    beta = (139.21_dp - (0.542_dp * mean_temperature))

    ! translation to a HL-model
    a0 = (grain_growth_rate * beta) / ICE_DENSITY
    a1 = a0

    ! activation energies
    ec0 = (                                                   &
    &  883.8_dp                               &
    &  * (abs(mean_temperature - TEMP_OFFSET)**(-0.885_dp)) &
    )
    ec1 = ec0

    ! call Herron & Langway model
    CALL tfm_density_HLtype(              &
    &  nz=nz,                             &
    &  dt=dt,                             &
    &  stage1_params=[ a0, ALPHA0, ec0 ], &
    &  stage2_params=[ a1, ALPHA1, ec1 ], &
    &  depth=depth,                       &
    &  temperature=temperature,           &
    &  density=density,                   &
    &  age=age,                           &
    &  d_density=d_density                &
    )
  END FUNCTION tfm_density_li2003


  FUNCTION tfm_density_helsen2008( &
  &  nz,                           &
  &  dt,                           &
  &  depth,                        &
  &  density,                      &
  &  temperature,                  &
  &  age,                          &
  &  grain_radius                  &
  ) RESULT(d_density)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_helsen2008
    !
    ! Helsen, M. M., van den Broeke, M. R., van den Wal, R. S. W.,
    ! van de Berg, W. J., van Meijgaard, E., Davis, C. H., Li, Y., and
    ! Goodwin, I. Elevation Changes in Antarctica Mainly Determined by
    ! Accumulation Variability. Science, 320 (5883), 1626-1629, (2008).
    ! https://doi.org/10.1126/science.1153894
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "depth", "density", "temperature", "age",
    !       "grain_radius", "d_density".
    !   dt: Time step (s).
    !   depth: Depth along the firn profile (m).
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temeperature along the firn profile (K).
    !   age: Age along the firn profile (s).
    !   grain_radius: Grain radius along the firn profile (m).
    !
    ! Result:
    !   d_density: Density change along the firn profile (kg m**-3).
    !---------------------------------------------------------------------------

    REAL(dp), PARAMETER :: &
      ! first stage parameters
      ALPHA0 = 1.0_dp,     &
      EC0    = 0.0_dp,     &
      ! second stage parameters
      ALPHA1 = 1.0_dp,     &
      EC1    = 0.0_dp

    REAL(dp) :: &
      a0,       &
      a1

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt
    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      density,                             &
      temperature,                         &
      age,                                 &
      grain_radius

    INTEGER :: n
    REAL(dp), DIMENSION(nz) :: d_density
    REAL(dp)                :: &
      mean_temperature,        &
      grain_growth_rate,       &
      beta

    !---------------------------------------------------------------------------

    CALL tfm_essentials_do_nothing(nz, grain_radius)

    ! 10 m temperature
    DO n = nz, 1, -1
      IF ( (depth(n) - depth(nz)) <= -10.0_dp ) EXIT
    END DO
    mean_temperature = temperature(n)

    ! grain growth rate and beta factor
    grain_growth_rate = (                                   &
    &  8.36_dp                                              &
    &  * (abs(mean_temperature - TEMP_OFFSET)**(-2.061_dp)) &
    )
    beta = (76.138_dp - (0.28965_dp * mean_temperature))

    ! translation to a HL-model
    a0 = (grain_growth_rate * beta) / ICE_DENSITY
    a1 = (grain_growth_rate * beta) / ICE_DENSITY

    ! call Herron & Langway model
    CALL tfm_density_HLtype(              &
    &  nz=nz,                             &
    &  dt=dt,                             &
    &  stage1_params=[ a0, ALPHA0, EC0 ], &
    &  stage2_params=[ a1, ALPHA1, EC1 ], &
    &  depth=depth,                       &
    &  temperature=temperature,           &
    &  density=density,                   &
    &  age=age,                           &
    &  d_density=d_density                &
    )
  END FUNCTION tfm_density_helsen2008


  FUNCTION tfm_density_arthern2010( &
  &  nz,                            &
  &  dt,                            &
  &  depth,                         &
  &  density,                       &
  &  temperature,                   &
  &  age,                           &
  &  grain_radius                   &
  ) RESULT(d_density)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_arthern2010
    !
    ! Arthern, R. J., Vaughan, D. G., Rankin, A. M., Mulvaney, R., and
    ! Thomas, E. R. In situ measurements of Antarctic snow compaction
    ! compared with predictions of models. Journal of Geophysical Research:
    ! Earth Surface, 115 (F3), (2010). https://doi.org/10.1029/2009JF001306
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "depth", "density", "temperature", "age",
    !       "grain_radius", "d_density".
    !   dt: Time step (s).
    !   depth: Depth along the firn profile (m).
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temeperature along the firn profile (K).
    !   age: Age along the firn profile (s).
    !   grain_radius: Grain radius along the firn profile (m).
    !
    ! Result:
    !   d_density: Density change along the firn profile (kg m**-3).
    !---------------------------------------------------------------------------

    REAL(dp), PARAMETER :: &
      ! first stage parameters
      ALPHA0 = 1.0_dp,     &
      EC0    = 60000.0_dp, &
      ! second stage parameters
      ALPHA1 = 1.0_dp,     &
      EC1    = 60000.0_dp, &
      ! further parameters
      EG = 42400.0_dp

    REAL(dp) :: &
      a0,       &
      a1

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt
    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      density,                             &
      temperature,                         &
      age,                                 &
      grain_radius

    REAL(dp), DIMENSION(nz) :: d_density
    REAL(dp)                :: mean_temperature
    INTEGER                 :: n

    !---------------------------------------------------------------------------

    CALL tfm_essentials_do_nothing(nz, grain_radius)

    ! 10 m temperature
    DO n = nz, 1, -1
      IF ( (depth(n) - depth(nz)) <= -10.0_dp ) EXIT
    END DO
    mean_temperature = temperature(n)

    ! factor depending on the mean annual temperature
    a0 = 0.07_dp * exp(EG / (GAS_CONST * mean_temperature))
    a1 = 0.03_dp * exp(EG / (GAS_CONST * mean_temperature))

    ! call Herron & Langway model
    CALL tfm_density_HLtype(              &
    &  nz=nz,                             &
    &  dt=dt,                             &
    &  stage1_params=[ a0, ALPHA0, EC0 ], &
    &  stage2_params=[ a1, ALPHA1, EC1 ], &
    &  depth=depth,                       &
    &  temperature=temperature,           &
    &  density=density,                   &
    &  age=age,                           &
    &  d_density=d_density                &
    )
  END FUNCTION tfm_density_arthern2010


  FUNCTION tfm_density_ligtenberg2011( &
  &  nz,                               &
  &  dt,                               &
  &  depth,                            &
  &  density,                          &
  &  temperature,                      &
  &  age,                              &
  &  grain_radius                      &
  ) RESULT(d_density)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_ligtenberg2011
    !
    ! Ligtenberg, S. R. M., Helsen, M. M. and van den Broeke, M. R. An
    ! improved semi-empirical model for the densification of Antarctic firn.
    ! The Cryosphere, 5, 809-819, (2011).
    ! https://doi.org/10.5194/tc-5-809-2011
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "depth", "density", "temperature", "age",
    !       "grain_radius", "d_density".
    !   dt: Time step (s).
    !   depth: Depth along the firn profile (m).
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temeperature along the firn profile (K).
    !   age: Age along the firn profile (s).
    !   grain_radius: Grain radius along the firn profile (m).
    !
    ! Result:
    !   d_density: Density change along the firn profile (kg m**-3).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt
    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      density,                             &
      temperature,                         &
      age,                                 &
      grain_radius

    INTEGER :: n
    REAL(dp), DIMENSION(nz) :: &
      d_density,               &
      mean_acc

    !---------------------------------------------------------------------------

    d_density = tfm_density_arthern2010( &
    &  nz=nz,                            &
    &  dt=dt,                            &
    &  depth=depth,                      &
    &  density=density,                  &
    &  temperature=temperature,          &
    &  age=age,                          &
    &  grain_radius=grain_radius         &
    )

    ! boundary between first and seconds stage
    DO n = nz, 1, -1
      IF ( density(n) >= 550.0_dp ) EXIT
    END DO

    ! change accumulation to (kg a-1 m-2)
    CALL tfm_essentials_mean_acc( &
    &  nz=nz,                     &
    &  depth=depth,               &
    &  density=density,           &
    &  age=age,                   &
    &  mean_acc=mean_acc          &
    )
    mean_acc = (mean_acc * WATER_DENSITY)

    d_density(n+1:nz-1) = (                                &
    &  d_density(n+1:nz-1)                                 &
    &  * (1.435_dp - (0.151_dp * log(mean_acc(n+1:nz-1)))) &
    )
    d_density(1:n) = (                                &
    &  d_density(1:n)                                 &
    &  * (2.366_dp - (0.293_dp * log(mean_acc(1:n)))) &
    )
  END FUNCTION tfm_density_ligtenberg2011


  FUNCTION tfm_density_simonsen2013( &
  &  nz,                             &
  &  dt,                             &
  &  depth,                          &
  &  density,                        &
  &  temperature,                    &
  &  age,                            &
  &  grain_radius                    &
  ) RESULT(d_density)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_density_simonsen2013
    !
    ! Simonsen, S. B., Stenseng, L., Adalgeisdottir, G., Fausto, R. S.,
    ! Hvidberg, C. S., and Lucas-Picher, P. Assessing a multilayered dynamic
    ! firn-compaction model for Greenland with ASIRAS radar measurements.
    ! Journal of Glaciology, 59 (215), 545-558, (2013).
    ! https://doi.org/10.3189/2013JoG12J158
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "depth", "density", "temperature", "age",
    !       "grain_radius", "d_density".
    !   dt: Time step (s).
    !   depth: Depth along the firn profile (m).
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temeperature along the firn profile (K).
    !   age: Age along the firn profile (s).
    !   grain_radius: Grain radius along the firn profile (m).
    !
    ! Result:
    !   d_density: Density change along the firn profile (kg m**-3).
    !---------------------------------------------------------------------------

    REAL(dp), PARAMETER :: &
      ! fist stage parameters (as implemented in CFM)
      F0 = 0.8_dp,         &
      ! second stage parameters (as implemented in CFm)
      F1 = 1.25_dp

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt
    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      density,                             &
      temperature,                         &
      age,                                 &
      grain_radius

    INTEGER :: n
    REAL(dp) :: mean_temperature
    REAL(dp), DIMENSION(nz) :: &
      d_density,               &
      mean_acc

    !---------------------------------------------------------------------------

    d_density = tfm_density_arthern2010( &
    &  nz=nz,                            &
    &  dt=dt,                            &
    &  depth=depth,                      &
    &  density=density,                  &
    &  temperature=temperature,          &
    &  age=age,                          &
    &  grain_radius=grain_radius         &
    )

    ! boundary between first and seconds stage
    DO n = nz, 1, -1
      IF ( density(n) >= 550.0_dp ) EXIT
    END DO

    ! change accumulation to (kg a-1 m-2)
    CALL tfm_essentials_mean_acc( &
    &  nz=nz,                     &
    &  depth=depth,               &
    &  density=density,           &
    &  age=age,                   &
    &  mean_acc=mean_acc          &
    )
    mean_acc = (mean_acc * WATER_DENSITY)

    ! 10 m temperature
    DO n = nz, 1, -1
      IF ( (depth(n) - depth(nz)) <= -10.0_dp ) EXIT
    END DO
    mean_temperature = temperature(n)

    d_density(n+1:nz-1) = F0 * d_density(n+1:nz-1)
    d_density(1:n) = (                                    &
    &  d_density(1:n)                                     &
    &  * (61.7_dp / (mean_acc(1:n)**0.5_dp))              &
    &  * exp(-3800.0_dp / (GAS_CONST * mean_temperature)) &
    )
  END FUNCTION tfm_density_simonsen2013
END MODULE tfm_density

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
MODULE tfm_constants
  !-----------------------------------------------------------------------------
  ! modules used
  !-----------------------------------------------------------------------------
  USE, INTRINSIC :: iso_fortran_env, ONLY: dp => real64

  !-----------------------------------------------------------------------------
  ! declarations
  !-----------------------------------------------------------------------------
  IMPLICIT NONE (TYPE, EXTERNAL)
  PRIVATE

  !-----------------------------------------------------------------------------
  ! public routines
  !-----------------------------------------------------------------------------
  PUBLIC ::                                  &
    ACC_GRAVITY,                             &
    GAS_CONST,                               &
    ICE_DENSITY,                             &
    WATER_DENSITY,                           &
    SECONDS_YEAR,                            &
    SECONDS_DAY,                             &
    SPECIFIC_HEAT_ICE,                       &
    LATENT_HEAT,                             &
    MELT_TEMP,                               &
    BOLTZMANN,                               &
    ICE_N,                                   &
    CRITICAL_DENSITY,                        &
    CLOSEOFF_DENSITY,                        &
    TEMP_OFFSET,                             &
    PI,                                      &
    VOLUME_H2O,                              &
    BURGERS_VECTOR,                          &
    PRE_FACTOR_LATTICE_DIFFUSION,            &
    ACTIVATION_ENERGY_LATTICE_DIFFUSION,     &
    PRE_FACTOR_BOUNDARY_DIFFUSION,           &
    ACTIVATION_ENERGY_BOUNDARY_DIFFUSION,    &
    PRE_FACTOR_DISLOCATION_CREEP_LOW,        &
    ACTIVATION_ENERGY_DISLOCATION_CREEP_LOW, &
    PRE_FACTOR_DISLOCATION_CREEP_HIGH,       &
    ACTIVATION_ENERGY_DISLOCATION_CREEP_HIGH


  REAL(dp), PARAMETER :: &
    ! acceleration due to gravity (kg / m s**2)
    ACC_GRAVITY = 9.81_dp, &
    ! gas constant (J / K mol)
    GAS_CONST = 8.31446261815324_dp, &
    ! ice density (kg / m**3)
    ICE_DENSITY = 917.0_dp, &
    ! water density (kg / m**3)
    WATER_DENSITY = 1000.0_dp, &
    ! seconds per year (s / yr)
    SECONDS_YEAR = (3600.0_dp * 24.0_dp * 365.0_dp), &
    ! seconds per day (s / d)
    SECONDS_DAY = (3600.0_dp * 24.0_dp), &
    ! specific heat capacity of ice (J / kg K) Reijmer et al. 2012
    SPECIFIC_HEAT_ICE = 2050.0_dp, &
    ! latent heat of ice (J / kg) Reijmer et al. 2012
    LATENT_HEAT = 334000.0_dp, &
    ! melt tempreature (K)
    MELT_TEMP = 273.15_dp, &
    ! Boltzmann constant
    BOLTZMANN = 1.380649E-23_dp, &
    ! ice flow exponent
    ICE_N = 3.0_dp, &
    ! critical density
    CRITICAL_DENSITY = 550.0_dp, &
    ! bubble close off density
    CLOSEOFF_DENSITY = 834.0_dp, &
    ! tempreature offset to convert K and C
    TEMP_OFFSET = 273.15_dp, &
    ! pi
    PI = 4.0d0 * datan(1.d0), &
    ! Volume of a H2O molecule (Maeno & Ebinuma, 1983)
    VOLUME_H2O = 3.27E-29_dp, &
    ! m ! Burgers vector of ice (Maeno & Ebinuma, 1983)
    BURGERS_VECTOR = 0.45E-9_dp

  ! pre factor and activation energy for lattice diffusion in ice (Maeno & Ebinnuma, 1983)
  REAL(dp), PARAMETER ::                             &
    PRE_FACTOR_LATTICE_DIFFUSION = 0.03_dp,          &
    ACTIVATION_ENERGY_LATTICE_DIFFUSION = 66200.0_dp ! J mol**-1

  ! pre factor and activatio energy for boundary diffusion in ice (Maeno & Ebinuma, 1983)
  REAL(dp), PARAMETER ::                              &
    PRE_FACTOR_BOUNDARY_DIFFUSION = 0.03_dp,          &
    ACTIVATION_ENERGY_BOUNDARY_DIFFUSION = 44100.0_dp ! J mol**-1

  ! pre factor and activation energy for dislocation creep in ice (Greve & Blatter, 2009)
  REAL(dp), PARAMETER ::                                   &
    PRE_FACTOR_DISLOCATION_CREEP_LOW = 3.985E-13_dp,       &
    ACTIVATION_ENERGY_DISLOCATION_CREEP_LOW = 60000.0_dp,  &
    PRE_FACTOR_DISLOCATION_CREEP_HIGH = 1.916E3_dp,        &
    ACTIVATION_ENERGY_DISLOCATION_CREEP_HIGH = 139000.0_dp
END MODULE tfm_constants

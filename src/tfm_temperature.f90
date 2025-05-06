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
MODULE tfm_temperature
  !-----------------------------------------------------------------------------
  ! modules used
  !-----------------------------------------------------------------------------
  USE, INTRINSIC :: iso_fortran_env, ONLY : dp => real64

  USE tfm_constants, ONLY : &
    MELT_TEMP,              &
    ICE_DENSITY

  USE tfm_essentials, ONLY : &
    tfm_essentials_do_nothing

  !-----------------------------------------------------------------------------
  ! declarations
  !-----------------------------------------------------------------------------
  IMPLICIT NONE (TYPE, EXTERNAL)
  PRIVATE

  !-----------------------------------------------------------------------------
  ! public routines
  !-----------------------------------------------------------------------------
  PUBLIC ::                                        &
    tfm_temperature_liquid_cond_geomMean,          &
    tfm_temperature_liquid_cond_voigt,             &
    tfm_temperature_sat_cond_geomMean,             &
    tfm_temperature_sat_cond_Voigt,                &
    tfm_temperature_sat_cond_Reuss,                &
    tfm_temperature_sat_cond_miller1969UpperBound, &
    tfm_temperature_conduct_Miller1969LowerBound,  &
    tfm_temperature_conduct_geomMean,              &
    tfm_temperature_capacity_Cuffey2010,           &
    tfm_temperature_conduct_Sturm1997,             &
    tfm_temperature_diffusion,                     &
    tfm_temperature_conduct_Marchenko2019,         &
    tfm_temperature_capacity_Paterson1994,         &
    tfm_temperature_conduct_Calonne2019,           &
    tfm_temperature_sat_cond_Miller1969LowerBound, &
    tfm_temperature_conduct_Miller1969UpperBound

  !-----------------------------------------------------------------------------
  ! interfaces
  !-----------------------------------------------------------------------------
  INTERFACE
    FUNCTION dry_thermcond_inter( &
    &  nz,                        &
    &  density,                   &
    &  temperature                &
    )
      USE, INTRINSIC :: iso_fortran_env, ONLY : dp => real64
      IMPLICIT NONE (TYPE, EXTERNAL)

      INTEGER, INTENT(IN)                 :: nz
      REAL(dp), DIMENSION(nz), INTENT(IN) :: &
        density,                             &
        temperature
      REAL(dp), DIMENSION(nz)             :: dry_thermcond_inter
    END FUNCTION dry_thermcond_inter
  END INTERFACE


  INTERFACE
    FUNCTION sat_thermcond_inter( &
    &  nz,                        &
    &  density                    &
    )
      USE, INTRINSIC :: iso_fortran_env, ONLY : dp => real64
      IMPLICIT NONE (TYPE, EXTERNAL)

      INTEGER, INTENT(IN)                 :: nz
      REAL(dp), DIMENSION(nz), INTENT(IN) :: density
      REAL(dp), DIMENSION(nz)             :: sat_thermcond_inter
    END FUNCTION sat_thermcond_inter
  END INTERFACE

  !-----------------------------------------------------------------------------
  ! routines of this module
  !-----------------------------------------------------------------------------
  CONTAINS


  FUNCTION tfm_temperature_diffusion( &
  &  nz,                              &
  &  dt,                              &
  &  depth,                           &
  &  density,                         &
  &  temperature,                     &
  &  heat_capacity,                   &
  &  thermal_conductivity             &
  ) result(d_temperature)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      density,                             &
      temperature,                         &
      heat_capacity,                       &
      thermal_conductivity

    REAL(dp), DIMENSION(nz) :: d_temperature
    INTEGER :: n
    REAL(dp) :: w

    REAL(dp), DIMENSION(nz)   :: &
      q,                         &
      ap,                        &
      n_temperature

    REAL(dp), DIMENSION(nz-1) :: &
      dz,                        &
      at,                        &
      ab

    REAL(dp), DIMENSION(nz-2) :: &
      dzt,                       &
      dzb,                       &
      dzc,                       &
      kp,                        &
      kt,                        &
      kb,                        &
      tkt,                       &
      tkb

    !---------------------------------------------------------------------------

    ! height changes
    dz = depth(2:nz) - depth(1:nz-1)
    dzt = dz(2:nz-1)
    dzb = dz(1:nz-2)
    dzc = (0.5_dp * (dzt + dzb))

    ! interface conductivity
    kp = thermal_conductivity(2:nz-1)
    tkt = thermal_conductivity(3:nz)
    tkb = thermal_conductivity(1:nz-2)
    kt = (2.0_dp * kp * tkt) / (kp + tkt)
    kb = (2.0_dp * kp * tkb) / (kp + tkb)

    ! coefficients
    at(2:nz-1) = -(kt / (dzt * dzc)) * (dt / (density(2:nz-1) * heat_capacity(2:nz-1)))
    ab(1:nz-2) = -(kb / (dzb * dzc)) * (dt / (density(2:nz-1) * heat_capacity(2:nz-1)))
    ap(2:nz-1) = 1.0 - (at(2:nz-1) + ab(1:nz-2))

    ! right hand side
    q = (1.0_dp * temperature)

    ! upper boundary condition
    ap(1) = 1.0_dp
    at(1) = 0.0_dp

    ! lower boundary condition
    ap(nz) = 1.0_dp
    ab(nz-1) = 0.0_dp

    ! Gauß elemination
    DO n = 2, nz, 1
      w = ab(n-1) / ap(n-1)
      ap(n) = ap(n) - (w * at(n-1))
      q(n) = q(n) - (w * q(n-1))
    END DO

    ! substitution
    n_temperature(nz) = q(nz) / ap(nz)
    DO n = nz - 1, 1, -1
      n_temperature(n) = (q(n) - (at(n) * n_temperature(n + 1))) / ap(n)
    END DO
    d_temperature = (n_temperature - temperature)
  END FUNCTION tfm_temperature_diffusion


  FUNCTION tfm_temperature_airCapacity_CRC( &
  &  nz,                                    &
  &  temperature                            &
  ) RESULT(air_heat_capacity)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_temperature_airCapacity_CRC
    !
    ! Heat capacity of air at given temperature and pressure of p = 0.1 MPa,
    ! computed by linear interplation from the vlaues listed in the CRC
    ! Handbook of Chemistriy and Physics.
    !
    ! Haynes, W. M., Lide, D. R., and Bruno, T. J. (Editors) (2016). Fluid
    ! Properties. In: CRC Handbook of Chemistry and Physics. CRC Press.
    !
    ! Author: Timm Schultz
    !
    ! Argument:
    !   Temperature (K).
    !
    ! Result:
    !   air_heat_capacity: Heat capacity of air at given temperature and
    !   p = 0.1 MPa (J kg**-1 K**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN)                 :: nz
    REAL(dp), DIMENSION(nz), INTENT(IN) :: temperature
    REAL(dp), DIMENSION(nz)             :: air_heat_capacity

    INTEGER :: &
      n,       &
      m

    REAL(dp), DIMENSION(25), PARAMETER :: REF_TEMP = [ &
    &  60.0_dp,                                        &
    &  78.79_dp,                                       &
    &  81.61_dp,                                       &
    &  100.0_dp,                                       &
    &  120.0_dp,                                       &
    &  140.0_dp,                                       &
    &  160.0_dp,                                       &
    &  180.0_dp,                                       &
    &  200.0_dp,                                       &
    &  220.0_dp,                                       &
    &  240.0_dp,                                       &
    &  260.0_dp,                                       &
    &  280.0_dp,                                       &
    &  300.0_dp,                                       &
    &  320.0_dp,                                       &
    &  340.0_dp,                                       &
    &  360.0_dp,                                       &
    &  380.0_dp,                                       &
    &  400.0_dp,                                       &
    &  500.0_dp,                                       &
    &  600.0_dp,                                       &
    &  700.0_dp,                                       &
    &  800.0_dp,                                       &
    &  900.0_dp,                                       &
    &  1000.0_dp                                       &
    ]
    REAL(dp), DIMENSION(25), PARAMETER :: REF_CAP = [ &
    &  1901.0_dp,                                     &
    &  1933.0_dp,                                     &
    &  1089.0_dp,                                     &
    &  1040.0_dp,                                     &
    &  1022.0_dp,                                     &
    &  1014.0_dp,                                     &
    &  1011.0_dp,                                     &
    &  1008.0_dp,                                     &
    &  1007.0_dp,                                     &
    &  1006.0_dp,                                     &
    &  1006.0_dp,                                     &
    &  1006.0_dp,                                     &
    &  1006.0_dp,                                     &
    &  1007.0_dp,                                     &
    &  1007.0_dp,                                     &
    &  1009.0_dp,                                     &
    &  1010.0_dp,                                     &
    &  1012.0_dp,                                     &
    &  1014.0_dp,                                     &
    &  1030.0_dp,                                     &
    &  1051.0_dp,                                     &
    &  1075.0_dp,                                     &
    &  1099.0_dp,                                     &
    &  1121.0_dp,                                     &
    &  1141.0_dp                                      &
    ]

    !---------------------------------------------------------------------------

    ! loop over all temperature values
    DO n = 1, nz, 1

      ! find relevant values for the interpolation
      DO m = 1, size(REF_TEMP), 1
        IF (REF_TEMP(m) > temperature(n)) THEN
          EXIT
        END IF
      END DO

      ! linear interpolation
      air_heat_capacity = (                                              &
      &  REF_CAP(m-1)                                                    &
      &  + (                                                             &
      &    ((REF_CAP(m) - REF_CAP(m-1)) / (REF_TEMP(m) - REF_TEMP(m-1))) &
      &    * (temperature(n) - REF_TEMP(m-1))                            &
      &  )                                                               &
      )
    END DO
  END FUNCTION tfm_temperature_airCapacity_CRC


  FUNCTION tfm_temperature_capacity_paterson1994( &
  &  nz,                                          &
  &  density,                                     &
  &  temperature,                                 &
  &  liquid_water                                 &
  ) RESULT(n_heat_capacity)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function tfm_temperature_capacity_Paterson1994
    !
    ! Heat capacity of firn following the mixture approach described in
    ! Kaviany (1991) (Voigt model weighting) and using the constant heat
    ! capacity of ice given by Paterson (1994).
    ! The heat capacity of water is that at 273.16 K and is taken from the
    ! CRC Handbook of Chemistry and Physics (Haynes, 2016).
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "density", "temperature", and
    !     "liquid_water".
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temperature along the firn profile (K).
    !   liquid_water: Volumetric liquid water content (1).
    !
    ! Result:
    !   n_heat_capacity: Effective heat capacity of firn (J kg**-1 K**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), intent(in) :: &
      density,                             &
      temperature,                         &
      liquid_water

    REAL(dp), DIMENSION(nz) :: n_heat_capacity

    REAL(dp), PARAMETER :: WATER_HEAT_CAPACITY_273 = 4219.4_dp

    REAL(dp), dimension(nz) :: &
      ice_cap,                 &
      air_cap,                 &
      wat_cap

    !---------------------------------------------------------------------------

    ice_cap = 2009.0_dp

    wat_cap = WATER_HEAT_CAPACITY_273

    air_cap = tfm_temperature_airCapacity_CRC( &
    &  nz=nz,                                  &
    &  temperature=temperature                 &
    )

    n_heat_capacity = (                                                &
    &  (ice_cap * (density / ICE_DENSITY))                             &
    &  + (wat_cap * liquid_water)                                      &
    &  + (air_cap * (1.0_dp - (density / ICE_DENSITY) - liquid_water)) &
    )
  END FUNCTION tfm_temperature_capacity_paterson1994


  FUNCTION tfm_temperature_capacity_Cuffey2010( &
  &  nz,                                        &
  &  density,                                   &
  &  temperature,                               &
  &  liquid_water                               &
  ) RESULT(n_heat_capacity)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function tfm_temperature_capacity_Cuffey2010
    !
    ! Heat capacity of firn following the mixture approach described in
    ! Kaviany (1991) (Voigt model weighting) and using the temperature
    ! dependentheat capacity of ice given by Cuffey & Paterson (2010).
    ! The heat capacity of water is that at 273.16 K and is taken from the
    ! CRC Handbook of Chemistry and Physics (Haynes, 2016).
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "density", "temperature", and
    !     "liquid_water".
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temperature along the firn profile (K).
    !   liquid_water: Volumetric liquid water content (1).
    !
    ! Result:
    !   n_heat_capacity: Effective heat capacity of firn (J kg**-1 K**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      density,                             &
      temperature,                         &
      liquid_water

    REAL(dp), DIMENSION(nz) :: n_heat_capacity

    REAL(dp), PARAMETER :: WATER_HEAT_CAPACITY_273 = 4219.4_dp

    REAL(dp), DIMENSION(nz) :: &
      ice_cap,                 &
      wat_cap,                 &
      air_cap

    !---------------------------------------------------------------------------

    ice_cap = (152.5_dp + (7.122_dp * temperature))

    wat_cap = WATER_HEAT_CAPACITY_273

    air_cap = tfm_temperature_airCapacity_CRC( &
    &  nz=nz,                                  &
    &  temperature=temperature                 &
    )

    n_heat_capacity = (                                                &
    &  (ice_cap * (density / ICE_DENSITY))                             &
    &  + (wat_cap * liquid_water)                                      &
    &  + (air_cap * (1.0_dp - (density / ICE_DENSITY) - liquid_water)) &
    )
  END FUNCTION tfm_temperature_capacity_Cuffey2010


  FUNCTION tfm_temperature_airConductivity_CRC( &
  &  nz,                                        &
  &  temperature                                &
  ) RESULT(air_thermal_conductivity)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_temperature_airConductivity_CRC
    !
    ! Thermal conductivity at given temperature and pressure of p = 0.1 MPa,
    ! computed by linear interpolation from the values listed in the CRC
    ! Handbook of Chemistry and Physics.
    !
    ! Haynes, W. M., Lide, D. R., and Bruno, T. J. (Editors) (2016). Fluid
    ! Properties. In: CRC Handbook of Chemistry and Physics. CRC Press.
    !
    ! Author: Timm Schultz
    !
    ! Argument:
    !   Temperature (K).
    !
    ! Result:
    !   air_thermal_conductivity: Thermal conductivity of air at given
    !     temperature and p = 0.1 MPa (W m**-1 K**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN)                 :: nz
    REAL(dp), DIMENSION(nz), INTENT(IN) :: temperature
    REAL(dp), DIMENSION(nz)             :: air_thermal_conductivity

    integer :: &
      n,       &
      m

    REAL(dp), DIMENSION(25), PARAMETER :: REF_TEMP = [ &
    &  60.0_dp,                                        &
    &  78.79_dp,                                       &
    &  81.61_dp,                                       &
    &  100.0_dp,                                       &
    &  120.0_dp,                                       &
    &  140.0_dp,                                       &
    &  160.0_dp,                                       &
    &  180.0_dp,                                       &
    &  200.0_dp,                                       &
    &  220.0_dp,                                       &
    &  240.0_dp,                                       &
    &  260.0_dp,                                       &
    &  280.0_dp,                                       &
    &  300.0_dp,                                       &
    &  320.0_dp,                                       &
    &  340.0_dp,                                       &
    &  360.0_dp,                                       &
    &  380.0_dp,                                       &
    &  400.0_dp,                                       &
    &  500.0_dp,                                       &
    &  600.0_dp,                                       &
    &  700.0_dp,                                       &
    &  800.0_dp,                                       &
    &  900.0_dp,                                       &
    &  1000.0_dp                                       &
    ]
    REAL(dp), DIMENSION(25), PARAMETER :: REF_COND = [ &
    &  0.1711_dp,                                      &
    &  0.1402_dp,                                      &
    &  0.7673_dp,                                      &
    &  0.9469_dp,                                      &
    &  1.138_dp,                                       &
    &  1.324_dp,                                       &
    &  1.505_dp,                                       &
    &  1.68_dp,                                        &
    &  1.85_dp,                                        &
    &  2.016_dp,                                       &
    &  2.177_dp,                                       &
    &  2.335_dp,                                       &
    &  2.488_dp,                                       &
    &  2.638_dp,                                       &
    &  2.785_dp,                                       &
    &  2.929_dp,                                       &
    &  3.071_dp,                                       &
    &  3.209_dp,                                       &
    &  3.345_dp,                                       &
    &  3.994_dp,                                       &
    &  4.601_dp,                                       &
    &  5.176_dp,                                       &
    &  5.725_dp,                                       &
    &  6.254_dp,                                       &
    &  6.768_dp                                        &
    ]

    !---------------------------------------------------------------------------

    ! loop over all temperature values
    DO n = 1, nz, 1

      ! find relevant values for the interpolation
      DO m = 1, size(REF_TEMP), 1
        IF (REF_TEMP(m) > temperature(n)) THEN
          EXIT
        END IF
      END DO

      ! linear interpolation
      air_thermal_conductivity(n) = (                                      &
      &  REF_COND(m-1)                                                     &
      &  + (                                                               &
      &    ((REF_COND(m) - REF_COND(m-1)) / (REF_TEMP(m) - REF_TEMP(m-1))) &
      &    * (temperature(n) - REF_TEMP(m-1))                              &
      &  )                                                                 &
      )
    END DO
  END FUNCTION tfm_temperature_airConductivity_CRC


  FUNCTION tfm_temperature_iceConductivity_Cuffey2010( &
  &  nz,                                               &
  &  temperature                                       &
  ) RESULT(thermal_conductivity)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_temperature_iceConductivity_Cuffey2010
    !
    ! Cuffey, K. M. and Paterson, W. S. B. (2010). The Physics of Glaciers.
    ! p. 400.
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimenstion of variable "temperature".
    !   temperature: Temperature (K).
    !
    ! Result:
    !   thermal_conductivity: Thermal conductivity of ice at given
    !     temperature (W m**-1 K**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN)                 :: nz
    REAL(dp), DIMENSION(nz), INTENT(IN) :: temperature
    REAL(dp), DIMENSION(nz)             :: thermal_conductivity

    !---------------------------------------------------------------------------

    thermal_conductivity = (9.828_dp * exp(-4.7D-3 * temperature))
  END FUNCTION tfm_temperature_iceConductivity_Cuffey2010


  FUNCTION tfm_temperature_mixture_miller1969UpperBound( &
  &  nz,                                                 &
  &  porosity,                                           &
  &  solid_cond,                                         &
  &  liquid_cond                                         &
  ) RESULT(eff_cond)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_temperature_mixutre_miller1969UpperBound
    !
    ! Miller, M. N. (1969). Bounds for Effective Electrical, Thermal, and
    ! Magenetic Properties of Heterogeneous Materials. Journal of
    ! Mathematical Physics, Volume 10, Number 11, pp. 1988-2004,
    ! https://doi.org/10.1063/1.1664794
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variabels "density" and "temperature".
    !   porosity: Porosity (1).
    !   solid_cond: Thermal conductivity of the solid phase (W m**-1 K**-1).
    !   liquid_cond: Thermal conductvitiy of the liquid phase
    !     (W m**-1 K**-1).
    !
    ! Result:
    !   eff_cond: Effective thermal conductivity (W m**-1 K**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      porosity,                            &
      solid_cond,                          &
      liquid_cond

    REAL(dp), DIMENSION(nz) :: eff_cond

    REAL(dp), PARAMETER     :: PARAM_G = (1.0_dp / 9.0_dp)
    REAL(dp), DIMENSION(nz) :: fu

    !---------------------------------------------------------------------------

    fu = (1.0_dp - (                                                        &
    &  (                                                                    &
    &    (1.0_dp - porosity)                                                &
    &    * (((solid_cond / liquid_cond) - 1.0_dp)**2.0_dp)                  &
    &    * porosity                                                         &
    &  )                                                                    &
    &  / (                                                                  &
    &    3.0_dp                                                             &
    &    * (                                                                &
    &      1.0_dp                                                           &
    &      + ((1.0_dp - porosity) * ((solid_cond / liquid_cond) - 1.0_dp))  &
    &    )                                                                  &
    &    * (                                                                &
    &      1.0_dp                                                           &
    &      + (                                                              &
    &        (1.0_dp - porosity)                                            &
    &        + ((3.0_dp * ((2.0_dp * porosity) - 1.0_dp)) * PARAM_G)        &
    &      )                                                                &
    &      * ((solid_cond / liquid_cond) - 1.0_dp)                          &
    &    )                                                                  &
    &  )                                                                    &
    ))

    eff_cond = (fu * (                                                 &
    &  1.0_dp                                                          &
    &  + ((1.0_dp - porosity) * ((solid_cond / liquid_cond) - 1.0_dp)) &
    ))
  END FUNCTION tfm_temperature_mixture_miller1969UpperBound


  FUNCTION tfm_temperature_mixture_miller1969LowerBound( &
  &  nz,                                                 &
  &  porosity,                                           &
  &  solid_cond,                                         &
  &  liquid_cond                                         &
  ) RESULT(eff_cond)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_temperature_mixutre_miller1969LowerBound
    !
    ! Miller, M. N. (1969). Bounds for Effective Electrical, Thermal, and
    ! Magenetic Properties of Heterogeneous Materials. Journal of
    ! Mathematical Physics, Volume 10, Number 11, pp. 1988-2004,
    ! https://doi.org/10.1063/1.1664794
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variabels "density" and "temperature".
    !   porosity: Porosity (1).
    !   solid_cond: Thermal conductivity of the solid phase (W m**-1 K**-1).
    !   liquid_cond: Thermal conductvitiy of the liquid phase
    !     (W m**-1 K**-1).
    !
    ! Result:
    !   eff_cond: Effective thermal conductivity (W m**-1 K**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      porosity,                            &
      solid_cond,                          &
      liquid_cond

    REAL(dp), DIMENSION(nz) :: eff_cond

    REAL(dp), PARAMETER :: PARAM_G = (1.0_dp / 9.0_dp)

    !---------------------------------------------------------------------------

    eff_cond = (                                                         &
    &  solid_cond * ((                                                   &
    &    (solid_cond / liquid_cond)                                      &
    &    - ((1.0_dp - porosity) * ((solid_cond / liquid_cond) - 1.0_dp)) &
    &    - (                                                             &
    &      (                                                             &
    &        (4.0_dp / 3.0_dp)                                           &
    &        * (((solid_cond / liquid_cond) - 1.0_dp)**2.0_dp)           &
    &        * (1.0_dp - porosity)                                       &
    &        * porosity                                                  &
    &      )                                                             &
    &      / (                                                           &
    &        1.0_dp                                                      &
    &        + (solid_cond / liquid_cond)                                &
    &        + (                                                         &
    &          (3.0_dp * (1.0_dp - (2.0_dp * porosity)))                 &
    &          * ((solid_cond / liquid_cond) - 1.0_dp)                   &
    &          * PARAM_G                                                 &
    &        )                                                           &
    &      )                                                             &
    &    )                                                               &
    &  )**(-1.0_dp))                                                     &
    )
  END FUNCTION tfm_temperature_mixture_miller1969LowerBound


  FUNCTION tfm_temperature_mixture_geomMean( &
  &  nz,                                     &
  &  porosity,                               &
  &  solid_cond,                             &
  &  liquid_cond                             &
  ) RESULT(eff_cond)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_temperature_mixture_geomMean
    !
    ! See:
    ! Kaviany, M. (1991). Principles of Heat Transfer in Porous Media. in:
    ! Ling, F. F. (edt.) Mechanical Engineering Series, Springer-Verlag, New
    ! York, Berlin, Heidelberg, https://doi.org/10.1063/1.1664794
    ! p. 126
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variabels "density" and "temperature".
    !   porosity: Porosity (1).
    !   solid_cond: Thermal conductivity of the solid phase (W m**-1 K**-1).
    !   liquid_cond: Thermal conductvitiy of the liquid phase
    !     (W m**-1 K**-1).
    !
    ! Result:
    !   eff_cond: Effective thermal conductivity (W m**-1 K**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      porosity,                            &
      solid_cond,                          &
      liquid_cond

    REAL(dp), DIMENSION(nz) :: eff_cond

    !---------------------------------------------------------------------------

    eff_cond = (                            &
    &  (liquid_cond**porosity)              &
    &  * (solid_cond**(1.0_dp - porosity))  &
    )
  END FUNCTION tfm_temperature_mixture_geomMean


  FUNCTION tfm_temperature_mixture_voigt( &
  &  nz,                                  &
  &  porosity,                            &
  &  solid_cond,                          &
  &  liquid_cond                          &
  ) RESULT(eff_cond)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function tfm_temperature_mixture_voigt
    !
    ! Voigt, W. (1889). Ueber die Beziehung zwischen den beiden
    ! Elasticitätsconstanten Isotroper Körper. Annalen der Physik, Volume
    ! 274, Issue 12, pp. 573-587,
    ! https://doi.org/10.1002/andp.18892741206
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variabels "density" and "temperature".
    !   porosity: Porosity (1).
    !   solid_cond: Thermal conductivity of the solid phase (W m**-1 K**-1).
    !   liquid_cond: Thermal conductvitiy of the liquid phase
    !     (W m**-1 K**-1).
    !
    ! Result:
    !   eff_cond: Effective thermal conductivity (W m**-1 K**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      porosity,                            &
      solid_cond,                          &
      liquid_cond

    REAL(dp), DIMENSION(nz) :: eff_cond

    !---------------------------------------------------------------------------

    eff_cond = (                            &
    &  (porosity * liquid_cond)             &
    &  + ((1.0_dp - porosity) * solid_cond) &
    )
  END FUNCTION tfm_temperature_mixture_voigt


  FUNCTION tfm_temperature_mixture_reuss( &
  &  nz,                                  &
  &  porosity,                            &
  &  solid_cond,                          &
  &  liquid_cond                          &
  ) RESULT(eff_cond)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function tfm_temperature_mixture_reuss
    !
    ! Reuss, A. (1929). Berechnung der Fließgrenze von Mischkristallen auf
    ! Grund der Platizitätsbedingung für Einkristalle. ZAMM - Journal of
    ! Applied Mathematics and Mechanics, Volume 9, Issue 1, pp. 49-58,
    ! https://doi.org/10.1002/zamm.19290090104
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variabels "density" and "temperature".
    !   porosity: Porosity (1).
    !   solid_cond: Thermal conductivity of the solid phase (W m**-1 K**-1).
    !   liquid_cond: Thermal conductvitiy of the liquid phase
    !     (W m**-1 K**-1).
    !
    ! Result:
    !   eff_cond: Effective thermal conductivity (W m**-1 K**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      porosity,                            &
      solid_cond,                          &
      liquid_cond

    REAL(dp), DIMENSION(nz) :: eff_cond

    !---------------------------------------------------------------------------

    eff_cond = (                                                       &
    &  solid_cond                                                      &
    &  / ((porosity * (solid_cond / liquid_cond)) + 1.0_dp - porosity) &
    )
  END FUNCTION tfm_temperature_mixture_reuss


  FUNCTION tfm_temperature_conduct_calonne2019( &
  &  nz,                                        &
  &  density,                                   &
  &  temperature                                &
  ) RESULT(n_thermal_conductivity)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_temperature_conduct_calonne2019
    !
    ! Calonne, N., Milliancourt, L., Burr, A., Philip, A., Martin, C. L.,
    ! Flin, F., and Geindreau, C. (2019). Thermal Conductivity of Snow,
    ! Firn, and Porous Ice From 3-D Image-Based Computations, Volume 46,
    ! Issue 22, pp. 13079-13089, https://doi.org/10.1029/2019GL085228
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variabels "density" and "temperature".
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temperature along the firn profile (K).
    !
    ! Result:
    !   n_thermal_conductivity: Thermal conductvitiy along the firn profile
    !     (W m**-1 K**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      density,                             &
      temperature

    REAL(dp), DIMENSION(nz) :: n_thermal_conductivity

    REAL(dp), DIMENSION(nz) :: &
      theta,                   &
      firn_cond,               &
      snow_cond,               &
      air_cond,                &
      ice_cond

    REAL(dp), PARAMETER ::           &
      TRANSITION_DENSITY = 450.0_dp, & ! kg m**-3
      ICE_COND_REF = 2.107_dp,       & ! W m**-1 K**-1
      AIR_COND_REF = 0.024_dp,       & ! W m**-1 K**-1
      PARAM_A = 0.02_dp                ! m**3 kf**-1

    !---------------------------------------------------------------------------

    theta = (                                                      &
    &  1.0_dp                                                      &
    &  / (                                                         &
    &    1.0_dp                                                    &
    &    + exp(-2.0_dp * PARAM_A * (density - TRANSITION_DENSITY)) &
    &  )                                                           &
    )

    snow_cond = (                     &
    &  0.024_dp                       &
    &  - (1.23D-4 * density)          &
    &  + (2.4D-6 * (density**2.0_dp)) &
    )

    firn_cond = (                                &
    &  2.107_dp                                  &
    &  + (0.003618_dp * (density - ICE_DENSITY)) &
    )

    air_cond = tfm_temperature_airConductivity_CRC( &
    &  nz=nz,                                       &
    &  temperature=temperature                      &
    )

    ice_cond = tfm_temperature_iceConductivity_Cuffey2010( &
    &  nz=nz,                                              &
    &  temperature=temperature                             &
    )

    n_thermal_conductivity = (                                   &
    &  (                                                         &
    &    (1.0_dp - theta) * snow_cond * (                        &
    &      (ice_cond * air_cond) / (ICE_COND_REF * AIR_COND_REF) &
    &    )                                                       &
    &  )                                                         &
    &  + (theta * firn_cond * (ice_cond / ICE_COND_REF))         &
    )
  END FUNCTION tfm_temperature_conduct_calonne2019


  FUNCTION tfm_temperature_conduct_marchenko2019( &
  &  nz,                                          &
  &  density,                                     &
  &  temperature                                  &
  ) RESULT(n_thermal_conductivity)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_temperature_conduct_marchenko2019
    !
    ! Marchenko, S., Cheng, G., Lötstedt, P., Pohjola, V., Pettersson, R.,
    ! van Pelt, W., and Reijmer, C. (2019). Thermal conductivity of firn at
    ! Lomonosovfonna, Svalbard, derived from subsurface temperature
    ! measurements. The Cryosphere, Volume 13, Issue 7, pp. 1843-1859,
    ! https://doi.org/10.5194/tc-13-1843-2019
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variabels "density" and "temperature".
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temperature along the firn profile (K).
    !
    ! Result:
    !   n_thermal_conductivity: Thermal conductvitiy along the firn profile
    !     (W m**-1 K**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      density,                             &
      temperature

    REAL(dp), DIMENSION(nz) :: n_thermal_conductivity

    !---------------------------------------------------------------------------

    CALL tfm_essentials_do_nothing(nz, temperature)

    n_thermal_conductivity = ((0301D-2 * density) - 0.724_dp)
  END FUNCTION tfm_temperature_conduct_marchenko2019


  FUNCTION tfm_temperature_conduct_sturm1997( &
  &  nz,                                      &
  &  density,                                 &
  &  temperature                              &
  ) RESULT(n_thermal_conductivity)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_temperature_conduct_sturm1997
    !
    ! Sturm, M., Holmgren, J., König, M., and Morris, K. (1997). The thermal
    ! conductivity of seasonal snow. Journal of Glaciology, Volume 43, Issue
    ! 143, pp. 26-41, https://doi.org/10.3189/S0022143000002781
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variabels "density" and "temperature".
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temperature along the firn profile (K).
    !
    ! Result:
    !   n_thermal_conductivity: Thermal conductvitiy along the firn profile
    !     (W m**-1 K**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      density,                             &
      temperature

    REAL(dp), DIMENSION(nz) :: n_thermal_conductivity

    !---------------------------------------------------------------------------

    CALL tfm_essentials_do_nothing(nz, temperature)

    n_thermal_conductivity = (            &
    &  (0.138_dp)                         &
    &  - ((1.010D-3) * density)           &
    &  + ((3.233D-6) * (density**2.0_dp)) &
    )
  END FUNCTION tfm_temperature_conduct_sturm1997


  FUNCTION tfm_temperature_conduct_miller1969UpperBound( &
  &  nz,                                                 &
  &  density,                                            &
  &  temperature                                         &
  ) RESULT(n_thermal_conductivity)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_temperature_conduct_miller1969UpperBound
    !
    ! Miller, M. N. (1969). Bounds for Effective Electrical, Thermal, and
    ! Magenetic Properties of Heterogeneous Materials. Journal of
    ! Mathematical Physics, Volume 10, Number 11, pp. 1988-2004,
    ! https://doi.org/10.1063/1.1664794
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variabels "density" and "temperature".
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temperature along the firn profile (K).
    !
    ! Result:
    !   n_thermal_conductivity: Thermal conductvitiy along the firn profile
    !     (W m**-1 K**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      density,                             &
      temperature

    REAL(dp), DIMENSION(nz) :: n_thermal_conductivity

    REAL(dp), DIMENSION(nz) :: &
      air_cond,                &
      ice_cond,                &
      porosity

    !---------------------------------------------------------------------------

    air_cond = tfm_temperature_airConductivity_CRC( &
    &  nz=nz,                                       &
    &  temperature=temperature                      &
    )

    ice_cond = tfm_temperature_iceConductivity_Cuffey2010( &
    &  nz=nz,                                              &
    &  temperature=temperature                             &
    )

    porosity = (1.0_dp - (density / ICE_DENSITY))

    n_thermal_conductivity = tfm_temperature_mixture_miller1969UpperBound( &
    &  nz=nz,                                                              &
    &  porosity=porosity,                                                  &
    &  solid_cond=ice_cond,                                                &
    &  liquid_cond=air_cond                                                &
    )
  END FUNCTION tfm_temperature_conduct_miller1969UpperBound


  FUNCTION tfm_temperature_conduct_miller1969LowerBound( &
  &  nz,                                                 &
  &  density,                                            &
  &  temperature                                         &
  ) RESULT(n_thermal_conductivity)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_temperature_conduct_miller1969LowerBound
    !
    ! Miller, M. N. (1969). Bounds for Effective Electrical, Thermal, and
    ! Magenetic Properties of Heterogeneous Materials. Journal of
    ! Mathematical Physics, Volume 10, Number 11, pp. 1988-2004,
    ! https://doi.org/10.1063/1.1664794
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variabels "density" and "temperature".
    !   porosity: Porosity (1).
    !   solid_cond: Thermal conductivity of the solid phase (W m**-1 K**-1).
    !   liquid_cond: Thermal conductvitiy of the liquid phase
    !     (W m**-1 K**-1).
    !
    ! Result:
    !   eff_cond: Effective thermal conductivity (W m**-1 K**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      density,                             &
      temperature

    REAL(dp), DIMENSION(nz) :: n_thermal_conductivity

    REAL(dp), PARAMETER :: PARAM_G = (1.0_dp / 9.0_dp)

    REAL(dp), DIMENSION(nz) :: &
      air_cond,                &
      ice_cond,                &
      porosity

    !---------------------------------------------------------------------------

    air_cond = tfm_temperature_airConductivity_CRC( &
    &  nz=nz,                                       &
    &  temperature=temperature                      &
    )

    ice_cond = tfm_temperature_iceConductivity_Cuffey2010( &
    &  nz=nz,                                              &
    &  temperature=temperature                             &
    )

    porosity = (1.0_dp - (density / ICE_DENSITY))

    n_thermal_conductivity = tfm_temperature_mixture_miller1969LowerBound( &
    &  nz=nz,                                                              &
    &  porosity=porosity,                                                  &
    &  solid_cond=ice_cond,                                                &
    &  liquid_cond=air_cond                                                &
    )
  END FUNCTION tfm_temperature_conduct_miller1969LowerBound


  FUNCTION tfm_temperature_conduct_geomMean( &
  &  nz,                                     &
  &  density,                                &
  &  temperature                             &
  ) RESULT(n_thermal_conductivity)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_temperature_conduct_geomMean
    !
    ! See:
    ! Kaviany, M. (1991). Principles of Heat Transfer in Porous Media. in:
    ! Ling, F. F. (edt.) Mechanical Engineering Series, Springer-Verlag, New
    ! York, Berlin, Heidelberg, https://doi.org/10.1063/1.1664794
    ! p. 126
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variabels "density" and "temperature".
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temperature along the firn profile (K).
    !
    ! Result:
    !   n_thermal_conductivity: Thermal conductvitiy along the firn profile
    !     (W m**-1 K**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      density,                             &
      temperature

    REAL(dp), DIMENSION(nz) :: n_thermal_conductivity

    REAL(dp), DIMENSION(nz) :: &
      air_cond,                &
      ice_cond,                &
      porosity

    !---------------------------------------------------------------------------

    air_cond = tfm_temperature_airConductivity_CRC( &
    &  nz=nz,                                       &
    &  temperature=temperature                      &
    )

    ice_cond = tfm_temperature_iceConductivity_Cuffey2010( &
    &  nz=nz,                                              &
    &  temperature=temperature                             &
    )

    porosity = (1.0_dp - (density / ICE_DENSITY))

    n_thermal_conductivity = tfm_temperature_mixture_geomMean( &
    &  nz=nz,                                                  &
    &  porosity=porosity,                                      &
    &  solid_cond=ice_cond,                                    &
    &  liquid_cond=air_cond                                    &
    )
  END FUNCTION tfm_temperature_conduct_geomMean


  FUNCTION tfm_temperature_sat_cond_Voigt( &
  &  nz,                                   &
  &  density                               &
  ) RESULT (sat_cond)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_temperature_sat_Voigt
    !
    ! Function calculating the effective thermal conductivity of firn at
    ! water saturation (needed for the computation of the effective thermal
    ! conductivity at unsaturated water levels). The effective thermal
    ! conductivity is computed using a the model of Voigt (1889). The firn
    ! temperature is assumed to be the melt temperature.
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "density" and "temperature".
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temperature along the firn profile (K).
    !
    ! Result:
    !   sat_cond: Effective thermal conductivity of firn at water
    !     saturation (W m**-1 K**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN)                 :: nz
    REAL(dp), DIMENSION(nz), INTENT(IN) :: density
    REAL(dp), DIMENSION(nz)             :: sat_cond

    ! Thermal conductivity of water at 273.15 K (W m**-1 K**-1)
    REAL(dp), PARAMETER :: WATER_COND_273 = 0.561_dp

    REAL(dp), DIMENSION(nz) :: &
      porosity,                &
      temperature,             &
      ice_cond,                &
      water_cond

    !---------------------------------------------------------------------------

    porosity = (1.0_dp - (density / ICE_DENSITY))
    temperature(:) = MELT_TEMP

    ice_cond = tfm_temperature_iceConductivity_Cuffey2010( &
    &  nz=nz,                                              &
    &  temperature=temperature                             &
    )

    water_cond(:) = WATER_COND_273

    sat_cond = tfm_temperature_mixture_Voigt( &
    &  nz=nz,                                 &
    &  porosity=porosity,                     &
    &  solid_cond=ice_cond,                   &
    &  liquid_cond=WATER_COND                 &
    )
  END FUNCTION tfm_temperature_sat_cond_Voigt


  FUNCTION tfm_temperature_sat_cond_Reuss( &
  &  nz,                                   &
  &  density                               &
  ) RESULT (sat_cond)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_temperature_sat_Reuss
    !
    ! Function calculating the effective thermal conductivity of firn at
    ! water saturation (needed for the computation of the effective thermal
    ! conductivity at unsaturated water levels). The effective thermal
    ! conductivity is computed using a the model of Reuss (1929). The firn
    ! temperature is assumed to be the melt temperature.
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "density" and "temperature".
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temperature along the firn profile (K).
    !
    ! Result:
    !   sat_cond: Effective thermal conductivity of firn at water
    !     saturation (W m**-1 K**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN)                 :: nz
    REAL(dp), DIMENSION(nz), INTENT(IN) :: density
    REAL(dp), DIMENSION(nz)             :: sat_cond

    ! Thermal conductivity of water at 273.15 K (W m**-1 K**-1)
    REAL(dp), PARAMETER :: WATER_COND_273 = 0.561_dp

    REAL(dp), DIMENSION(nz) :: &
      porosity,                &
      temperature,             &
      ice_cond,                &
      water_cond

    !---------------------------------------------------------------------------

    porosity = (1.0_dp - (density / ICE_DENSITY))
    temperature(:) = MELT_TEMP

    ice_cond = tfm_temperature_iceConductivity_Cuffey2010( &
    &  nz=nz,                                              &
    &  temperature=temperature                             &
    )

    water_cond(:) = WATER_COND_273

    sat_cond = tfm_temperature_mixture_Reuss( &
    &  nz=nz,                                 &
    &  porosity=porosity,                     &
    &  solid_cond=ice_cond,                   &
    &  liquid_cond=WATER_COND                 &
    )
  END FUNCTION tfm_temperature_sat_cond_Reuss


  FUNCTION tfm_temperature_sat_cond_Miller1969UpperBound( &
  &  nz,                                                  &
  &  density                                              &
  ) RESULT (sat_cond)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_temperature_sat_Miller1969UpperBound
    !
    ! Function calculating the effective thermal conductivity of firn at
    ! water saturation (needed for the computation of the effective thermal
    ! conductivity at unsaturated water levels). The effective thermal
    ! conductivity is computed using a the upper bound model of Miller
    ! (1969). The firn temperature is assumed to be the melt temperature.
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "density" and "temperature".
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temperature along the firn profile (K).
    !
    ! Result:
    !   sat_cond: Effective thermal conductivity of firn at water
    !     saturation (W m**-1 K**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN)                 :: nz
    REAL(dp), DIMENSION(nz), INTENT(IN) :: density
    REAL(dp), DIMENSION(nz)             :: sat_cond

    ! Thermal conductivity of water at 273.15 K (W m**-1 K**-1)
    REAL(dp), PARAMETER :: WATER_COND_273 = 0.561_dp

    REAL(dp), DIMENSION(nz) :: &
      porosity,                &
      temperature,             &
      ice_cond,                &
      water_cond

    !---------------------------------------------------------------------------

    porosity = (1.0_dp - (density / ICE_DENSITY))
    temperature(:) = MELT_TEMP
    ice_cond = tfm_temperature_iceConductivity_Cuffey2010( &
    &  nz=nz,                                              &
    &  temperature=temperature                             &
    )
    water_cond(:) = WATER_COND_273

    sat_cond = tfm_temperature_mixture_Miller1969UpperBound( &
    &  nz=nz,                                                &
    &  porosity=porosity,                                    &
    &  solid_cond=ice_cond,                                  &
    &  liquid_cond=WATER_COND                                &
    )
  END FUNCTION tfm_temperature_sat_cond_Miller1969UpperBound


  FUNCTION tfm_temperature_sat_cond_Miller1969LowerBound( &
  &  nz,                                                  &
  &  density                                              &
  ) RESULT (sat_cond)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_temperature_sat_Miller1969LowerBound
    !
    ! Function calculating the effective thermal conductivity of firn at
    ! water saturation (needed for the computation of the effective thermal
    ! conductivity at unsaturated water levels). The effective thermal
    ! conductivity is computed using a the lower bound model of Miller
    ! (1969). The firn temperature is assumed to be the melt temperature.
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "density" and "temperature".
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temperature along the firn profile (K).
    !
    ! Result:
    !   sat_cond: Effective thermal conductivity of firn at water
    !     saturation (W m**-1 K**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN)                 :: nz
    REAL(dp), DIMENSION(nz), INTENT(IN) :: density
    REAL(dp), DIMENSION(nz)             :: sat_cond

    ! Thermal conductivity of water at 273.15 K (W m**-1 K**-1)
    REAL(dp), PARAMETER :: WATER_COND_273 = 0.561_dp

    REAL(dp), DIMENSION(nz) :: &
      porosity,                &
      temperature,             &
      ice_cond,                &
      water_cond

    !---------------------------------------------------------------------------

    porosity = (1.0_dp - (density / ICE_DENSITY))
    temperature(:) = MELT_TEMP
    ice_cond = tfm_temperature_iceConductivity_Cuffey2010( &
    &  nz=nz,                                              &
    &  temperature=temperature                             &
    )
    water_cond(:) = WATER_COND_273

    sat_cond = tfm_temperature_mixture_Miller1969LowerBound( &
    &  nz=nz,                                                &
    &  porosity=porosity,                                    &
    &  solid_cond=ice_cond,                                  &
    &  liquid_cond=WATER_COND                                &
    )
  END FUNCTION tfm_temperature_sat_cond_Miller1969LowerBound


  FUNCTION tfm_temperature_sat_cond_geomMean( &
  &  nz,                                      &
  &  density                                  &
  ) RESULT (sat_cond)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_temperature_sat_geomMean
    !
    ! Function calculating the effective thermal conductivity of firn at
    ! water saturation (needed for the computation of the effective thermal
    ! conductivity at unsaturated water levels). The effective thermal
    ! conductivity is computed using a geometric mean weighting approach
    ! (following Kaviany, 1991). The firn temperature is assumed to be the
    ! melt temperature.
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "density" and "temperature".
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temperature along the firn profile (K).
    !
    ! Result:
    !   sat_cond: Effective thermal conductivity of firn at water
    !     saturation (W m**-1 K**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN)                 :: nz
    REAL(dp), DIMENSION(nz), INTENT(IN) :: density
    REAL(dp), DIMENSION(nz)             :: sat_cond

    ! Thermal conductivity of water at 273.15 K (W m**-1 K**-1)
    REAL(dp), PARAMETER :: WATER_COND_273 = 0.561_dp

    REAL(dp), DIMENSION(nz) :: &
      porosity,                &
      temperature,             &
      ice_cond,                &
      water_cond

    !---------------------------------------------------------------------------

    porosity = (1.0_dp - (density / ICE_DENSITY))
    temperature(:) = MELT_TEMP
    ice_cond = tfm_temperature_iceConductivity_Cuffey2010( &
    &  nz=nz,                                              &
    &  temperature=temperature                             &
    )
    water_cond(:) = WATER_COND_273

    sat_cond = tfm_temperature_mixture_geomMean( &
    &  nz=nz,                                    &
    &  porosity=porosity,                        &
    &  solid_cond=ice_cond,                      &
    &  liquid_cond=WATER_COND                    &
    )
  END FUNCTION tfm_temperature_sat_cond_geomMean


  FUNCTION tfm_temperature_liquid_cond_geomMean( &
  &  nz,                                         &
  &  density,                                    &
  &  temperature,                                &
  &  liquid_water,                               &
  &  dry_thermcond_model,                        &
  &  sat_thermcond_model                         &
  ) RESULT(n_thermal_conductivity)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_temperature_liquid_cond_geomMean
    !
    ! The function computes the effective thermal conductivity of
    ! unsaturated firn following Kaviany (1991), using geometric mean
    ! weighting between dry and saturated conditions.
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "density", "temperature", and
    !     "liuqid_water".
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temperature along the firn profile (K).
    !   liquid_water: Liquid water content along the profile (1).
    !   dry_thermcond_model: Model to use for the computation of effective
    !     thermal conductvity at dry conditions.
    !   sat_thermcond_model: Model to use for the computation of effective
    !     thermal conductviity at saturated conditions.
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      density,                             &
      temperature,                         &
      liquid_water

    PROCEDURE(dry_thermcond_inter), POINTER :: dry_thermcond_model
    PROCEDURE(sat_thermcond_inter), POINTER :: sat_thermcond_model

    REAL(dp), DIMENSION(nz) :: n_thermal_conductivity

    REAL(dp), DIMENSION(nz) :: &
      dry_cond,                &
      sat_cond

    !---------------------------------------------------------------------------

    dry_cond = dry_thermcond_model( &
    &  nz=nz,                       &
    &  density=density,             &
    &  temperature=temperature      &
    )
    sat_cond = sat_thermcond_model( &
    &  nz=nz,                       &
    &  density=density              &
    )
    n_thermal_conductivity = tfm_temperature_mixture_geomMean( &
    &  nz=nz,                                                  &
    &  porosity=liquid_water,                                  &
    &  solid_cond=dry_cond,                                    &
    &  liquid_cond=sat_cond                                    &
    )
  END FUNCTION tfm_temperature_liquid_cond_geomMean


  FUNCTION tfm_temperature_liquid_cond_Voigt( &
  &  nz,                                      &
  &  density,                                 &
  &  temperature,                             &
  &  liquid_water,                            &
  &  dry_thermcond_model,                     &
  &  sat_thermcond_model                      &
  ) RESULT(n_thermal_conductivity)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_temperature_liquid_cond_Voigt
    !
    ! The function computes the effective thermal conductivity of
    ! unsaturated firn following Kaviany (1991), using the model of
    ! Voigt (1889).
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "density", "temperature", and
    !     "liuqid_water".
    !   density: Density along the firn profile (kg m**-3).
    !   temperature: Temperature along the firn profile (K).
    !   liquid_water: Liquid water content along the profile (1).
    !   dry_thermcond_model: Model to use for the computation of effective
    !     thermal conductvity at dry conditions.
    !   sat_thermcond_model: Model to use for the computation of effective
    !     thermal conductviity at saturated conditions.
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      density,                             &
      temperature,                         &
      liquid_water

    PROCEDURE(dry_thermcond_inter), POINTER :: dry_thermcond_model
    PROCEDURE(sat_thermcond_inter), POINTER :: sat_thermcond_model

    REAL(dp), DIMENSION(nz) :: n_thermal_conductivity

    REAL(dp), DIMENSION(nz) :: &
      dry_cond,                &
      sat_cond

    !---------------------------------------------------------------------------

    dry_cond = dry_thermcond_model( &
    &  nz=nz,                       &
    &  density=density,             &
    &  temperature=temperature      &
    )
    sat_cond = sat_thermcond_model( &
    &  nz=nz,                       &
    &  density=density              &
    )
    n_thermal_conductivity = tfm_temperature_mixture_Voigt( &
    &  nz=nz,                                               &
    &  porosity=liquid_water,                               &
    &  solid_cond=dry_cond,                                 &
    &  liquid_cond=sat_cond                                 &
    )
  END FUNCTION tfm_temperature_liquid_cond_Voigt
END MODULE Tfm_temperature

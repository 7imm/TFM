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
MODULE tfm_num
  !-----------------------------------------------------------------------------
  ! modules used
  !-----------------------------------------------------------------------------
  USE, INTRINSIC :: iso_fortran_env, ONLY: dp => real64

  USE tfm_liquid, ONLY :         &
    van_genuchten_inter,         &
    tfm_liquid_bucket,           &
    tfm_liquid_RichardsEquation, &
    vgParametersDaanen2009,      &
    vgParametersYamaguchi2012,   &
    vgParametersYamaguchi2010

  USE tfm_temperature, ONLY :                      &
    tfm_temperature_diffusion,                     &
    tfm_temperature_capacity_paterson1994,         &
    tfm_temperature_capacity_Cuffey2010,           &
    tfm_temperature_conduct_sturm1997,             &
    tfm_temperature_conduct_calonne2019,           &
    tfm_temperature_conduct_marchenko2019,         &
    tfm_temperature_conduct_miller1969UpperBound,  &
    tfm_temperature_conduct_miller1969LowerBound,  &
    tfm_temperature_conduct_geomMean,              &
    tfm_temperature_liquid_cond_geomMean,          &
    tfm_temperature_liquid_cond_voigt,             &
    tfm_temperature_liquid_cond_geomMean,          &
    tfm_temperature_sat_cond_geomMean,             &
    tfm_temperature_sat_cond_Voigt,                &
    tfm_temperature_sat_cond_reuss,                &
    tfm_temperature_sat_cond_Miller1969UpperBound, &
    tfm_temperature_sat_cond_Miller1969LowerBound

  USE tfm_grain, ONLY :      &
    tfm_grain_arthern2010,   &
    tfm_grain_zwally2002,    &
    tfm_grain_brun1989,      &
    tfm_grain_tusima1978,    &
    tfm_grain_katsushima2009

  USE tfm_density, ONLY :        &
    tfm_density_medley2020,      &
    tfm_density_herron1980,      &
    tfm_density_arthern2010,     &
    tfm_density_ligtenberg2011,  &
    tfm_density_simonsen2013,    &
    tfm_density_arthern1998,     &
    tfm_density_li2003,          &
    tfm_density_helsen2008,      &
    tfm_density_breant2017,      &
    tfm_density_zwinger2007,     &
    tfm_density_greve2009,       &
    tfm_density_gagliardini1998, &
    tfm_density_timmsfit,        &
    tfm_density_sintering,       &
    tfm_density_depth

  USE tfm_constants, ONLY : &
    SECONDS_YEAR,           &
    CLOSEOFF_DENSITY,       &
    WATER_DENSITY

  USE tfm_llStructure, ONLY : &
    llProps,                  &
    llGetLast,                &
    llGetData,                &
    llUpdateList,             &
    llDropData,               &
    llAppendData,             &
    llPropsDropData

  !-----------------------------------------------------------------------------
  ! declarations
  !-----------------------------------------------------------------------------
  IMPLICIT NONE (TYPE, EXTERNAL)
  PRIVATE

  !-----------------------------------------------------------------------------
  ! public routines
  !-----------------------------------------------------------------------------
  PUBLIC ::                    &
    sim_models,                &
    tfm_num_modelinit,         &
    tfm_num_trimprofilelength, &
    tfm_num_trimprofileage,    &
    tfm_num_surface,           &
    tfm_num_step

  !-----------------------------------------------------------------------------
  ! interfaces
  !-----------------------------------------------------------------------------
  ! interface for liquid water
  INTERFACE
    SUBROUTINE liquid_inter( &
    &  nz,                   &
    &  dt,                   &
    &  depth,                &
    &  density,              &
    &  temperature,          &
    &  grain_radius,         &
    &  water_content,        &
    &  liquid_accumulation,  &
    &  runoff,               &
    &  van_genuchten_model   &
    )
      USE, INTRINSIC :: iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE (TYPE, EXTERNAL)

      INTEGER, INTENT(IN) :: &
        nz

      REAL(dp), INTENT(IN) :: &
        dt

      REAL(dp), DIMENSION(nz), INTENT(IN) :: &
        depth,                               &
        grain_radius

      REAL(dp), INTENT(IN) :: &
        liquid_accumulation

      PROCEDURE(van_genuchten_inter), POINTER :: &
        van_genuchten_model

      REAL(dp), DIMENSION(nz), INTENT(INOUT) :: &
        density,                                &
        temperature,                            &
        water_content

      REAL(dp), INTENT(INOUT) :: &
        runoff
    END SUBROUTINE liquid_inter
  END INTERFACE


  ! interface for densification function
  INTERFACE
    FUNCTION density_inter( &
    &  nz,                  &
    &  dz,                  &
    &  depth,               &
    &  density,             &
    &  temperature,         &
    &  age,                 &
    &  grain_radius         &
    )
      USE, INTRINSIC :: iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE (TYPE, EXTERNAL)

      INTEGER, INTENT(IN) :: &
        nz

      REAL(dp), INTENT(IN) :: &
        dz

      REAL(dp), DIMENSION(nz), INTENT(IN) :: &
        depth,                               &
        density,                             &
        temperature,                         &
        age,                                 &
        grain_radius

      REAL(dp), DIMENSION(nz) :: &
        density_inter
    END FUNCTION density_inter
  END INTERFACE


  ! interface for temperature
  INTERFACE
    FUNCTION temperature_inter( &
    &  nz,                      &
    &  dt,                      &
    &  depth,                   &
    &  density,                 &
    &  temperature,             &
    &  heat_capacity,           &
    &  thermal_conductivity     &
    )
      USE, INTRINSIC :: iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE (TYPE, EXTERNAL)

      INTEGER, INTENT(IN) :: &
        nz

      REAL(dp), INTENT(IN) :: &
        dt

      REAL(dp), DIMENSION(nz), INTENT(IN) :: &
        depth,                               &
        density,                             &
        temperature,                         &
        heat_capacity,                       &
        thermal_conductivity

      REAL(dp), DIMENSION(nz) :: &
        temperature_inter
    END FUNCTION temperature_inter
  END INTERFACE


  ! interface for heat capacity
  INTERFACE
    FUNCTION heatcap_inter( &
    &  nz,                  &
    &  density,             &
    &  temperature,         &
    &  liquid_water         &
    )
      USE, INTRINSIC :: iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE (TYPE, EXTERNAL)

      INTEGER, INTENT(IN) :: &
        nz

      REAL(dp), DIMENSION(nz), intent(in) :: &
        density,                             &
        temperature,                         &
        liquid_water

      REAL(dp), DIMENSION(nz) :: &
        heatcap_inter
    END FUNCTION heatcap_inter
  END INTERFACE


  ! interface for dry firn thermal conductivity
  INTERFACE
    FUNCTION thermcond_inter( &
    &  nz,                    &
    &  density,               &
    &  temperature            &
    )
      USE, INTRINSIC :: iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE (TYPE, EXTERNAL)

      INTEGER, INTENT(IN) :: &
        nz

      REAL(dp), DIMENSION(nz), INTENT(IN) :: &
        density,                             &
        temperature

      REAL(dp), DIMENSION(nz) :: &
        thermcond_inter
    END FUNCTION thermcond_inter
  END INTERFACE


  ! interface for saturated thermal conductivity
  INTERFACE
    FUNCTION saturation_thermcond_inter( &
    &  nz,                               &
    &  density                           &
    )
      USE, INTRINSIC :: iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE (TYPE, EXTERNAL)

      INTEGER, INTENT(IN) :: &
        nz

      REAL(dp), DIMENSION(nz), INTENT(IN) :: &
        density

      REAL(dp), DIMENSION(nz) :: &
        saturation_thermcond_inter
    END FUNCTION saturation_thermcond_inter
  END INTERFACE


  ! interface for liquid thermal conductivity
  INTERFACE
    FUNCTION liquid_thermcond_inter( &
    &  nz,                           &
    &  density,                      &
    &  temperature,                  &
    &  liquid_water,                 &
    &  thermcond_model,              &
    &  sat_thermcond_model           &
    )
      USE, INTRINSIC :: iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE (TYPE, EXTERNAL)

      INTEGER, INTENT(IN) :: &
        nz

      REAL(dp), DIMENSION(nz), INTENT(IN) :: &
        density,                             &
        temperature,                         &
        liquid_water

      PROCEDURE(thermcond_inter), POINTER :: &
        thermcond_model

      PROCEDURE(saturation_thermcond_inter), POINTER :: &
        sat_thermcond_model

      REAL(dp), DIMENSION(nz) :: &
        liquid_thermcond_inter
    END FUNCTION liquid_thermcond_inter
  END INTERFACE


  ! interface for grain growth
  INTERFACE
    FUNCTION grain_growth_inter( &
    &  nz,                       &
    &  dt,                       &
    &  temperature,              &
    &  density,                  &
    &  liquid_water,             &
    &  grain_radius              &
    )
      USE, INTRINSIC :: iso_fortran_env, ONLY: dp => real64
      IMPLICIT NONE (TYPE, EXTERNAL)

      INTEGER, INTENT(IN) :: &
        nz

      REAL(dp), INTENT(IN) :: &
        dt

      REAL(dp), DIMENSION(nz), INTENT(IN) :: &
        temperature,                         &
        density,                             &
        liquid_water,                        &
        grain_radius

      REAL(dp), DIMENSION(nz) :: &
        grain_growth_inter
    END FUNCTION grain_growth_inter
  END INTERFACE


  !-----------------------------------------------------------------------------
  ! types
  !-----------------------------------------------------------------------------
  TYPE sim_models
    PROCEDURE(density_inter),              POINTER, NOPASS :: dens_model             => null()
    PROCEDURE(temperature_inter),          POINTER, NOPASS :: temp_model             => null()
    PROCEDURE(heatcap_inter),              POINTER, NOPASS :: heatcap_model          => null()
    PROCEDURE(thermcond_inter),            POINTER, NOPASS :: thermcond_model        => null()
    PROCEDURE(liquid_thermcond_inter),     POINTER, NOPASS :: liquid_thermcond_model => null()
    PROCEDURE(saturation_thermcond_inter), POINTER, NOPASS :: sat_thermcond_model    => null()
    PROCEDURE(liquid_inter),               POINTER, NOPASS :: liquid_model           => null()
    PROCEDURE(grain_growth_inter),         POINTER, NOPASS :: grain_model            => null()
    PROCEDURE(van_genuchten_inter),        POINTER, NOPASS :: van_genuchten_model    => null()
  END TYPE sim_models


  !-----------------------------------------------------------------------------
  ! routines of this module
  !-----------------------------------------------------------------------------
  contains


  SUBROUTINE tfm_num_modelinit(             &
  &  solve_density,                         &
  &  solve_temperature,                     &
  &  solve_heat_capacity,                   &
  &  solve_thermal_conductivity,            &
  &  solve_liquid_thermal_conductivity,     &
  &  solve_saturation_thermal_conductivity, &
  &  solve_liquid,                          &
  &  solve_van_genuchten,                   &
  &  solve_grain_growth,                    &
  &  models                                 &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)

    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: &
      solve_density,                          &
      solve_temperature,                      &
      solve_heat_capacity,                    &
      solve_thermal_conductivity,             &
      solve_liquid_thermal_conductivity,      &
      solve_saturation_thermal_conductivity,  &
      solve_liquid,                           &
      solve_van_genuchten,                    &
      solve_grain_growth

    ! pointer definition
    TYPE(sim_models), INTENT(INOUT) :: models

    ! default / fallback
    models%dens_model             => null()
    models%temp_model             => null()
    models%grain_model            => null()
    models%heatcap_model          => null()
    models%thermcond_model        => null()
    models%liquid_thermcond_model => null()
    models%sat_thermcond_model    => null()
    models%liquid_model           => null()
    models%van_genuchten_model    => null()

    IF ( present(solve_density) ) THEN
      IF ( solve_density == 'false' ) THEN
        models%dens_model => null()
      ELSE IF ( solve_density == 'medley2020' ) THEN
        models%dens_model => tfm_density_medley2020
      ELSE IF ( solve_density == 'herron1980' ) THEN
        models%dens_model => tfm_density_herron1980
      ELSE IF ( solve_density == 'arthern2010' ) THEN
        models%dens_model => tfm_density_arthern2010
      ELSE IF ( solve_density == 'ligtenberg2011' ) THEN
        models%dens_model => tfm_density_ligtenberg2011
      ELSE IF ( solve_density == 'simonsen2013' ) THEN
        models%dens_model => tfm_density_simonsen2013
      ELSE IF ( solve_density == 'arthern1998' ) THEN
        models%dens_model => tfm_density_arthern1998
      ELSE IF ( solve_density == 'li2003' ) THEN
        models%dens_model => tfm_density_li2003
      ELSE IF ( solve_density == 'helsen2008' ) THEN
        models%dens_model => tfm_density_helsen2008
      ELSE IF ( solve_density == 'breant2017' ) THEN
        models%dens_model => tfm_density_breant2017
      ELSE IF ( solve_density == 'zwinger2007' ) THEN
        models%dens_model => tfm_density_zwinger2007
      ELSE IF ( solve_density == 'greve2009' ) THEN
        models%dens_model => tfm_density_greve2009
      ELSE IF ( solve_density == 'gagliardini1998' ) THEN
        models%dens_model => tfm_density_gagliardini1998
      ELSE IF ( solve_density == 'timmsfit' ) THEN
        models%dens_model => tfm_density_timmsfit
      ELSE IF ( solve_density == 'sintering' ) THEN
        models%dens_model => tfm_density_sintering
      ELSE
        PRINT *, 'module: tfm_num'
        PRINT *, 'subroutine: tfm_num_modelinit'
        PRINT *, 'can not find densification model: ', solve_density
        PRINT *, 'stopping right here!'
        STOP
      END IF
    END IF

    IF ( present(solve_temperature) ) THEN
      IF ( solve_temperature == 'false' ) THEN
        models%temp_model => null()
      ELSE IF ( solve_temperature == 'true' ) THEN
        models%temp_model => tfm_temperature_diffusion
      ELSE
        PRINT *, 'module: tfm_num'
        PRINT *, 'subroutine: tfm_num_modelinit'
        PRINT *, 'can not find temperature model: ', solve_temperature
        PRINT *, 'stopping right here!'
        STOP
      END IF
    END IF

    IF ( present(solve_liquid) ) THEN
      IF ( solve_liquid == 'false' ) THEN
        models%liquid_model => null()
      ELSE IF ( solve_liquid == 'bucket' ) THEN
        models%liquid_model => tfm_liquid_bucket
      ELSE IF ( solve_liquid == 'richards_equation' ) THEN
        models%liquid_model => tfm_liquid_RichardsEquation

        IF ( present(solve_van_genuchten) ) THEN
          IF ( solve_van_genuchten == 'false' ) THEN
            PRINT *, 'module: tfm_num'
            PRINT *, 'subroutine: tfm_num_modelinit'
            PRINT *, 'method "richards_equation" is defined for'
            PRINT *, 'solving liquid water, but no model for the'
            PRINT *, 'van Genuchten parameters is given!'
            PRINT *, 'stopping right here!'
            STOP
          ELSE IF ( solve_van_genuchten == 'daanen2009' ) THEN
            models%van_genuchten_model => vgParametersDaanen2009
          ELSE IF ( solve_van_genuchten == 'yamaguchi2010' ) THEN
            models%van_genuchten_model => vgParametersYamaguchi2010
          ELSE IF ( solve_van_genuchten == 'yamaguchi2012' ) THEN
            models%van_genuchten_model => vgParametersYamaguchi2012
          ELSE
            PRINT *, 'module: tfm_num'
            PRINT *, 'subroutine: tfm_num_modelinit'
            PRINT *, 'can not find van Genuchten model:', solve_van_genuchten
            PRINT *, 'stopping right here!'
            STOP
          END IF

        ELSE
          PRINT *, 'module: tfm_num'
          PRINT *, 'subroutine: tfm_num_modelinit'
          PRINT *, 'methods "richards_equation" is defined for'
          PRINT *, 'solving liquid water, but no model for the'
          PRINT *, 'van Genuchten parameters is given!'
          PRINT *, 'stopping right here!'
          STOP
        END IF

      ELSE
        PRINT *, 'module: tfm_num'
        PRINT *, 'subroutine: tfm_num_modelinit'
        PRINT *, 'can not find liquid model: ', solve_liquid
        PRINT *, 'stopping right here!'
        STOP
      END IF
    END IF

    IF ( present(solve_heat_capacity) ) THEN
      IF ( solve_heat_capacity == 'false' ) THEN
        models%heatcap_model => null()
      ELSE IF ( solve_heat_capacity == 'paterson1994' ) THEN
        models%heatcap_model => tfm_temperature_capacity_paterson1994
      ELSE IF ( solve_heat_capacity == 'cuffey2010' ) THEN
        models%heatcap_model => tfm_temperature_capacity_Cuffey2010
      ELSE
        PRINT *, 'module: tfm_num'
        PRINT *, 'subroutine: tfm_num_modelinit'
        PRINT *, 'can not find heat capacity model: ', solve_heat_capacity
        PRINT *, 'stopping right here!'
        STOP
      END IF
    END IF

    IF ( present(solve_thermal_conductivity) ) THEN
      IF ( solve_thermal_conductivity == 'false' ) THEN
        models%thermcond_model => null()
      ELSE IF ( solve_thermal_conductivity == 'sturm1997' ) THEN
        models%thermcond_model => tfm_temperature_conduct_sturm1997
      ELSE IF ( solve_thermal_conductivity == 'calonne2019' ) THEN
        models%thermcond_model => tfm_temperature_conduct_calonne2019
      ELSE IF ( solve_thermal_conductivity == 'marchenko2019' ) THEN
        models%thermcond_model => tfm_temperature_conduct_marchenko2019
      ELSE IF ( solve_thermal_conductivity == 'miller1969upperbound' ) THEN
        models%thermcond_model => tfm_temperature_conduct_miller1969UpperBound
      ELSE IF ( solve_thermal_conductivity == 'miller1969lowerbound' ) THEN
        models%thermcond_model => tfm_temperature_conduct_miller1969LowerBound
      ELSE IF ( solve_thermal_conductivity == 'geometricmean' ) THEN
        models%thermcond_model => tfm_temperature_conduct_geomMean
      ELSE
        PRINT *, 'module: tfm_num'
        PRINT *, 'subroutine: tfm_num_modelinit'
        PRINT *, 'can not find thermal conductivity model: ', solve_heat_capacity
        PRINT *, 'stopping right here!'
        STOP
      END IF
    END IF

    IF ( present(solve_liquid_thermal_conductivity) ) THEN
      IF ( solve_liquid_thermal_conductivity == 'false' ) THEN
        models%liquid_thermcond_model => null()
      ELSE IF ( solve_liquid_thermal_conductivity == 'geometricmean' ) THEN
        models%liquid_thermcond_model => tfm_temperature_liquid_cond_geomMean
      ELSE IF ( solve_liquid_thermal_conductivity == 'voigt' ) THEN
        models%liquid_thermcond_model => tfm_temperature_liquid_cond_voigt
      ELSE
        PRINT *, 'module: tfm_num'
        PRINT *, 'subroutine: tfm_num_modelinit'
        PRINT *, 'can not find liquid thermal conductivity mode: ', &
          solve_liquid_thermal_conductivity
        PRINT *, 'stopping right here!'
        STOP
      END IF

      IF (                                                    &
      &  present(solve_saturation_thermal_conductivity)       &
      &  .AND. (solve_liquid_thermal_conductivity /= 'false') &
      ) THEN

        IF ( solve_saturation_thermal_conductivity == 'geometricmean' ) THEN
          models%sat_thermcond_model => tfm_temperature_sat_cond_geomMean
        ELSE IF ( solve_saturation_thermal_conductivity == 'voigt' ) THEN
          models%sat_thermcond_model => tfm_temperature_sat_cond_Voigt
        ELSE IF ( solve_saturation_thermal_conductivity == 'reuss' ) THEN
          models%sat_thermcond_model => tfm_temperature_sat_cond_Reuss
        ELSE IF ( solve_saturation_thermal_conductivity == 'miller1969upperbound' ) THEN
          models%sat_thermcond_model => tfm_temperature_sat_cond_Miller1969UpperBound
        ELSE IF ( solve_saturation_thermal_conductivity == 'miller1969lowerbound' ) THEN
          models%sat_thermcond_model => tfm_temperature_sat_cond_Miller1969LowerBound
        ELSE
          PRINT *, 'module: tfm_num'
          PRINT *, 'subroutine: tfm_num_modelinit'
          PRINT *, 'There is a model defined for solving the thermal'
          PRINT *, 'conductivity in presence of liquid water. But '
          PRINT *, 'there is no model defined how to solve the thermal'
          PRINT *, 'conductivity at water saturation.'
          PRINT *, 'Stopping right here!'
          STOP
        END IF
      END IF
    END IF

    IF ( present(solve_grain_growth) ) THEN
      IF ( solve_grain_growth == 'false' ) THEN
        models%grain_model => null()
      ELSE IF ( solve_grain_growth == 'arthern2010' ) THEN
        models%grain_model => tfm_grain_arthern2010
      ELSE IF ( solve_grain_growth == 'zwally2002' ) THEN
        models%grain_model => tfm_grain_zwally2002
      ELSE IF ( solve_grain_growth == 'brun1989' ) THEN
        models%grain_model => tfm_grain_brun1989
      ELSE IF ( solve_grain_growth == 'tusima1978' ) THEN
        models%grain_model => tfm_grain_tusima1978
      ELSE IF ( solve_grain_growth == 'katsushima2009' ) THEN
        models%grain_model => tfm_grain_katsushima2009
      ELSE
        PRINT *, 'module: tfm_num'
        PRINT *, 'subroutine: tfm_num_modelinit'
        PRINT *, 'can not find grain growth model: ', solve_grain_growth
        PRINT *, 'stopping right here!'
        STOP
      END IF
    END IF
  END SUBROUTINE tfm_num_modelinit


  SUBROUTINE tfm_num_step( &
  &  dt,                   &
  &  models,               &
  &  props,                &
  &  runoff,               &
  &  liquid_acc            &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)

    REAL(dp), INTENT(IN)              :: dt
    TYPE(sim_models), INTENT(IN)      :: models
    TYPE(llProps), INTENT(INOUT)      :: props
    REAL(dp), INTENT(INOUT), OPTIONAL :: runoff
    REAL(dp), INTENT(IN), OPTIONAL    :: liquid_acc

    INTEGER :: &
      nz,      &
      n,       &
      m
    REAL(dp), DIMENSION(6)  :: residuum

    REAL(dp), DIMENSION(props%depth%length) :: &
      depth,                                   &
      density,                                 &
      temperature,                             &
      grain_radius,                            &
      heat_capacity,                           &
      thermal_conductivity,                    &
      liquidwater,                             &
      age,                                     &
      !
      n_depth,                                 &
      n_density,                               &
      n_temperature,                           &
      n_grain_radius,                          &
      n_heat_capacity,                         &
      n_thermal_conductivity,                  &
      n_liquidwater,                           &
      n_age,                                   &
      !
      d_depth,                                 &
      d_density,                               &
      d_temperature,                           &
      d_grain_radius,                          &
      dens_residuum

    ! initialization
    nz                   = props%depth%length
    depth                = llGetData(props%depth)
    density              = llGetData(props%density)
    temperature          = llGetData(props%temperature)
    grain_radius         = llGetData(props%grain_radius)
    heat_capacity        = llGetData(props%heatcap)
    thermal_conductivity = llGetData(props%thermcond)
    liquidwater          = llGetDAta(props%liquidwater)
    age                  = llGetData(props%age)

    ! liquid model
    !if ( (associated(models%liquid_model)) .and. (liquid_acc > 0.0) ) then
    IF ( associated(models%liquid_model) ) THEN
      CALL models%liquid_model(                         &
      &  nz=nz,                                         &
      &  dt=dt,                                         &
      &  depth=depth,                                   &
      &  density=density,                               &
      &  temperature=temperature,                       &
      &  grain_radius=grain_radius,                     &
      &  water_content=liquidwater,                     &
      &  liquid_accumulation=liquid_acc,                &
      &  runoff=runoff,                                 &
      &  van_genuchten_model=models%van_genuchten_model &
      )
    END IF

    n_depth                = depth
    n_density              = density
    n_temperature          = temperature
    n_grain_radius         = grain_radius
    n_heat_capacity        = heat_capacity
    n_thermal_conductivity = thermal_conductivity
    n_liquidwater          = liquidwater
    n_age                  = age

    d_density      = 0.0_dp
    d_temperature  = 0.0_dp
    d_grain_radius = 0.0_dp
    residuum       = 0.0_dp
    residuum(1)    = -9999.9_dp

    ! Picard loop
    n = 0
    DO WHILE ( (maxval(abs(residuum)) > 1.0e-2_dp) .AND. (n < 100) )

      ! density model
      IF ( associated(models%dens_model) ) THEN
        d_density = models%dens_model( &
        &  nz, dt,                     &
        &  depth=n_depth,              &
        &  density=n_density,          &
        &  temperature=n_temperature,  &
        &  age=n_age,                  &
        &  grain_radius=grain_radius   &
        )

        ! There is the possibility that the residuum of the density is
        ! always high because despite the density is converging. This
        ! happens due to the discontinuous function describing
        ! densification. The density "flickers" around the value of
        ! 550 kg m-3. Therefore the residuum at this density is forced
        ! to zero, (which is not ideal).
        dens_residuum = abs(n_density - (density + d_density))
        DO m = 1, nz, 1
          IF ( floor(n_density(m)) == 550.0_dp       ) dens_residuum(m) = 0.0_dp
          IF ( floor(n_density(m)) == CLOSEOFF_DENSITY ) dens_residuum(m) = 0.0_dp
        END DO
        residuum(1) = maxval(dens_residuum)
      END IF

      ! heat capacity model
      IF ( associated(models%heatcap_model) ) THEN
        n_heat_capacity = models%heatcap_model( &
        &  nz,                                  &
        &  n_density,                           &
        &  n_temperature,                       &
        &  n_liquidwater                        &
        )
        !residuum(2) = maxval(abs(n_heat_capacity - heatcap))
      END IF

      ! thermal conductivity model
      IF ( associated(models%thermcond_model) ) THEN
        IF (                                                 &
        &  (any(n_liquidwater > 0.0_dp))                     &
        &  .AND. (associated(models%liquid_thermcond_model)) &
        ) THEN
          n_thermal_conductivity = models%liquid_thermcond_model( &
          &  nz,                                                  &
          &  n_density,                                           &
          &  n_temperature,                                       &
          &  n_liquidwater,                                       &
          &  models%thermcond_model,                              &
          &  models%sat_thermcond_model                           &
          )
        ELSE
          n_thermal_conductivity = models%thermcond_model( &
          &  nz,                                           &
          &  n_density,                                    &
          &  n_temperature                                 &
          )
        END IF
        !residuum(3) = maxval(abs(n_thermal_conductivity - thermcond))
      END IF

      ! temperature model
      IF ( associated(models%temp_model) ) THEN
        d_temperature = models%temp_model(             &
        &  nz, dt,                                     &
        &  depth=n_depth,                              &
        &  density=n_density,                          &
        &  temperature=temperature,                    &
        &  heat_capacity=n_heat_capacity,              &
        &  thermal_conductivity=n_thermal_conductivity &
        )
        residuum(4) = maxval(abs(                        &
        &  n_temperature - (temperature + d_temperature) &
        ))
      END IF

      ! grain growth model
      IF ( associated(models%grain_model) ) THEN
        d_grain_radius = models%grain_model( &
        &  nz, dt,                           &
        &  temperature=n_temperature,        &
        &  density=n_density,                &
        &  liquid_water=n_liquidwater,       &
        &  grain_radius=grain_radius         &
        )
        residuum(5) = maxval(abs(                           &
        &  n_grain_radius - (grain_radius + d_grain_radius) &
        ))
      END IF

      ! depth evolution
      d_depth = tfm_density_depth( &
      &  nz,                       &
      &  depth=depth,              &
      &  density=density,          &
      &  d_density=d_density       &
      )
      residuum(6) = maxval(abs(n_depth - (depth + d_depth)))

      ! reassignment
      n_depth                = (depth + d_depth)
      n_density              = (density + d_density)
      n_temperature          = (temperature + d_temperature)
      n_grain_radius         = (grain_radius + d_grain_radius)
      n_heat_capacity        = n_heat_capacity
      n_thermal_conductivity = n_thermal_conductivity

      n = n + 1
    END DO

    ! raise the age
    n_age = tfm_num_age(nz, dt, age)

    ! new value
    CALL llUpdateList(   &
    &  self=props%depth, &
    &  data=n_depth      &
    )
    CALL llUpdateList(     &
    &  self=props%density, &
    &  data=n_density      &
    )
    CALL llUpdateList(         &
    &  self=props%temperature, &
    &  data=n_temperature      &
    )
    CALL llUpdateList(          &
    &  self=props%grain_radius, &
    &  data=n_grain_radius      &
    )
    CALL llUpdateList(      &
    &  self=props%heatcap,  &
    &  data=n_heat_capacity &
    )
    CALL llUpdateList(             &
    &  self=props%thermcond,       &
    &  data=n_thermal_conductivity &
    )
    CALL llUpdateList( &
    &  self=props%age, &
    &  data=n_age      &
    )
    CALL llUpdateList(         &
    &  self=props%liquidwater, &
    &  data=n_liquidwater      &
    )
  END SUBROUTINE tfm_num_step


  SUBROUTINE tfm_num_surface( &
  &  dt,                      &
  &  forcing,                 &
  &  models,                  &
  &  props                    &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)

    REAL(dp), INTENT(IN)               :: dt
    REAL(dp), DIMENSION(7), INTENT(IN) :: forcing
    TYPE(sim_models), INTENT(IN)       :: models

    TYPE(llProps), INTENT(INOUT) :: props

    REAL(dp), DIMENSION(props%depth%length) :: &
      depth,                                   &
      density,                                 &
      temperature,                             &
      heatcap,                                 &
      thermcond,                               &
      grain_radius,                            &
      liquidwater,                             &
      age

    INTEGER :: &
      nz,      &
      n

    REAL(dp) ::   &
      dz,         &
      dm,         &
      am,         &
      surf_dens,  &
      surf_temp,  &
      surf_grain, &
      solid_acc

    nz         = props%depth%length
    surf_temp  = forcing(3)
    surf_dens  = forcing(4)
    solid_acc  = forcing(5)
    surf_grain = forcing(7)

    ! mass to be removed or added
    dm = solid_acc * dt * WATER_DENSITY

    ! the accumulation is zero
    IF ( dm == 0.0_dp ) THEN
      RETURN

    ! theres accumulation and there are still elements available
    ELSE IF ( dm > 0.0_dp ) THEN

      ! height change computed from surface density
      dz = dm / surf_dens

      ! with very small accumulation there can occure dpison problem
      ! causing two layers with the same depth
      IF ( ((llGetLast(props%depth) + dz) - llGetLast(props%depth)) <= 1.0E-10_dp ) THEN
      !if ( ((llGetLast(props%depth) + dz) - llGetLast(props%depth)) == 0.0_dp ) then
        RETURN
      END IF

      CALL llAppendData(                      &
      &  self=props%depth,                    &
      &  n=1,                                 &
      &  data=[ llGetLast(props%depth) + dz ] &
      )
      CALL llAppendData(     &
      &  self=props%density, &
      &  n=1,                &
      &  data=[ surf_dens ]  &
      )
      CALL llAppendData(         &
      &  self=props%temperature, &
      &  n=1,                    &
      &  data=[ surf_temp ]      &
      )
      CALL llAppendData(          &
      &  self=props%grain_radius, &
      &  n=1,                     &
      &  data=[ surf_grain ]      &
      )
      CALL llAppendData(                                 &
      &  self=props%heatcap,                             &
      &  n=1,                                            &
      &  data=[ models%heatcap_model(                    &
      &    nz=1,                                         &
      &    density=[ llGetLast(props%density) ],         &
      &    temperature=[ llGetLast(props%temperature) ], &
      &    liquid_water=[ llGetLast(props%liquidwater) ] &
      &  ) ]                                             &
      )
      CALL llAppendData(                                 &
      &  self=props%thermcond,                           &
      &  n=1,                                            &
      &  data=[ models%thermcond_model(                  &
      &     nz=1,                                        &
      &     density=[ llGetLast(props%density) ],        &
      &     temperature=[ llGetLast(props%temperature) ] &
      &  ) ]                                             &
      )
      CALL llAppendData(         &
      &  self=props%liquidwater, &
      &  n=1,                    &
      &  data=[ 0.0_dp ]         &
      )
      CALL llAppendData( &
      &  self=props%age, &
      &  n=1,            &
      &  data=[ 0.0_dp ] &
      )

    ! theres ablation
    ELSE IF ( dm < 0.0_dp ) THEN

      depth        = llGetData(props%depth)
      density      = llGetData(props%density)
      temperature  = llGetData(props%temperature)
      heatcap      = llGetData(props%heatcap)
      thermcond    = llGetData(props%thermcond)
      grain_radius = llGetData(props%grain_radius)
      liquidwater  = llGetData(props%liquidwater)
      age          = llGetData(props%age)

      ! removal of layers
      am = -dm
      DO n = nz - 1, 1, -1
        am = am - ((depth(n+1) - depth(n)) * density(n))
        if ( am <= 0.0_dp ) EXIT
        dm = dm + ((depth(n+1) - depth(n)) * density(n))
      END DO
      dz = dm / density(n)
      n = n + 1

      ! interpolation
      CALL tfm_num_lin_interp( &
      &  z0=depth(n),          &
      &  z1=depth(n-1),        &
      &  v0=temperature(n),    &
      &  v1=temperature(n-1),  &
      &  dz=dz,                &
      &  v=temperature(n)      &
      )
      CALL tfm_num_lin_interp( &
      &  z0=depth(n),          &
      &  z1=depth(n-1),        &
      &  v0=grain_radius(n),   &
      &  v1=grain_radius(n-1), &
      &  dz=dz,                &
      &  v=grain_radius(n)     &
      )
      CALL tfm_num_lin_interp( &
      &  z0=depth(n),          &
      &  z1=depth(n-1),        &
      &  v0=liquidwater(n),    &
      &  v1=liquidwater(n-1),  &
      &  dz=dz,                &
      &  v=liquidwater(n)      &
      )
      CALL tfm_num_lin_interp( &
      &  z0=depth(n),          &
      &  z1=depth(n-1),        &
      &  v0=heatcap(n),        &
      &  v1=heatcap(n-1),      &
      &  dz=dz,                &
      &  v=heatcap(n)          &
      )
      CALL tfm_num_lin_interp( &
      &  z0=depth(n),          &
      &  z1=depth(n-1),        &
      &  v0=thermcond(n),      &
      &  v1=thermcond(n-1),    &
      &  dz=dz,                &
      &  v=thermcond(n)        &
      )
      CALL tfm_num_lin_interp( &
      &  z0=depth(n),          &
      &  z1=depth(n-1),        &
      &  v0=age(n),            &
      &  v1=age(n-1),          &
      &  dz=dz,                &
      &  v=age(n)              &
      )

      ! new value
      CALL llUpdateList(   &
      &  self=props%depth, &
      &  data=depth        &
      )
      CALL llUpdateList(     &
      &  self=props%density, &
      &  data=density        &
      )
      CALL llUpdateList(        &
      &  self=props%temperature, &
      &  data=temperature        &
      )
      CALL llUpdateList(          &
      &  self=props%grain_radius, &
      &  data=grain_radius        &
      )
      CALL llUpdateList(     &
      &  self=props%heatcap, &
      &  data=heatcap        &
      )
      CALL llUpdateList(       &
      &  self=props%thermcond, &
      &  data=thermcond        &
      )
      CALL llUpdateList( &
      &  self=props%age, &
      &  data=age        &
      )
      CALL llUpdateList(         &
      &  self=props%liquidwater, &
      &  data=liquidwater        &
      )

      ! fill removed layers with NaN value
      CALL llDropData(     &
      &  self=props%depth, &
      &  n=-(nz - n)       &
      )
      CALL llDropData(       &
      &  self=props%density, &
      &  n=-(nz - n)         &
      )
      CALL llDropData(           &
      &  self=props%temperature, &
      &  n=-(nz - n)             &
      )
      CALL llDropData(       &
      &  self=props%heatcap, &
      &  n=-(nz - n)         &
      )
      CALL llDropData(         &
      &  self=props%thermcond, &
      &  n=-(nz - n)           &
      )
      CALL llDropData(            &
      &  self=props%grain_radius, &
      &  n=-(nz - n)              &
      )
      CALL llDropData(           &
      &  self=props%liquidwater, &
      &  n=-(nz - n)             &
      )
      CALL llDropData(   &
      &  self=props%age, &
      &  n=-(nz - n)     &
      )

      ! new depth / height of the uppermost layer
      !props%depth%tail%data(props%depth%tind - 1) = (     &
      !&  props%depth%tail%data(props%depth%tind - 1) + dz &
      !)
      depth(n) = (depth(n) + dz)
      CALL llUpdateList(   &
      &  self=props%depth, &
      &  data=depth        &
      )
    END IF
  END SUBROUTINE tfm_num_surface


  SUBROUTINE tfm_num_lin_interp( &
  &  z0,                         &
  &  z1,                         &
  &  v0,                         &
  &  v1,                         &
  &  dz,                         &
  &  v                           &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)

    REAL(dp), INTENT(IN) :: &
      z0,                   &
      z1,                   &
      v0,                   &
      v1,                   &
      dz

    REAL(dp), INTENT(INOUT) :: v

    v = v0 + ((v1 - v0) / (z1 - z0)) * dz
  END SUBROUTINE tfm_num_lin_interp


  FUNCTION tfm_num_age( &
  &  nz,                &
  &  dt,                &
  &  age                &
  ) RESULT (n_age)
    IMPLICIT NONE (TYPE, EXTERNAL)

    INTEGER, INTENT(IN)                 :: nz
    REAL(dp), INTENT(IN)                :: dt
    REAL(dp), DIMENSION(nz), INTENT(IN) :: age
    REAL(dp), DIMENSION(nz)             :: n_age

    n_age = (age + dt)
  END FUNCTION tfm_num_age


  SUBROUTINE tfm_num_trimProfileLength( &
  &  props,                             &
  &  length                             &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)

    TYPE(llProps), INTENT(INOUT) :: props
    REAL(dp), INTENT(IN)         :: length

    INTEGER                                 :: n
    REAL(dp), DIMENSION(props%depth%length) :: depth

    depth = llGetData(props%depth)
    DO n = 1, size(depth), 2
      IF ( (depth(size(depth)) - depth(n)) <= length ) EXIT
    END DO
    n = (n - 1)

    CALL llPropsDropData( &
    &  self=props,        &
    &  n=n                &
    )
  END SUBROUTINE tfm_num_trimProfileLength


  SUBROUTINE tfm_num_trimProfileAge( &
  &  props,                          &
  &  max_age                         &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)

    TYPE(llProps), INTENT(INOUT) :: props
    REAL(dp), INTENT(IN)         :: max_age

    INTEGER                               :: n
    REAL(dp), DIMENSION(props%age%length) :: age

    age = (llGetData(props%age) / SECONDS_YEAR)
    DO n = 1, size(age), 2
      IF ( age(n) <= max_age ) EXIT
    END DO
    n = (n - 1)

    CALL llPropsDropData( &
    &  self=props,        &
    &  n=n                &
    )
  END SUBROUTINE tfm_num_trimProfileAge
END MODULE tfm_num

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
MODULE tfm_liquid
  !-----------------------------------------------------------------------------
  ! modules used
  !-----------------------------------------------------------------------------
  USE, INTRINSIC :: iso_fortran_env, ONLY: dp => real64

  USE tfm_essentials, ONLY : &
    tfm_essentials_do_nothing

  USE tfm_constants, ONLY : &
    WATER_DENSITY,          &
    ACC_GRAVITY,            &
    ICE_DENSITY,            &
    LATENT_HEAT,            &
    SPECIFIC_HEAT_ICE,      &
    MELT_TEMP

  !-----------------------------------------------------------------------------
  ! declarations
  !-----------------------------------------------------------------------------
  IMPLICIT NONE (TYPE, EXTERNAL)
  PRIVATE

  !-----------------------------------------------------------------------------
  ! public routines
  !-----------------------------------------------------------------------------
  PUBLIC ::                      &
    van_genuchten_inter,         &
    tfm_liquid_richardsequation, &
    tfm_liquid_bucket,           &
    vgParametersDaanen2009,      &
    vgParametersYamaguchi2010,   &
    vgParametersYamaguchi2012

  !-----------------------------------------------------------------------------
  ! types
  !-----------------------------------------------------------------------------
  TYPE vanGenuchtenParameters
    REAL(dp), DIMENSION(:), ALLOCATABLE :: &
      alpha,                               &
      n,                                   &
      m
  END TYPE vanGenuchtenParameters

  !-----------------------------------------------------------------------------
  ! interfaces
  !-----------------------------------------------------------------------------
  ! interface for van Genuchten models
  INTERFACE
    SUBROUTINE van_genuchten_inter( &
    &  nz,                          &
    &  density,                     &
    &  grain_radius,                &
    &  vg_params                    &
    )
      USE, INTRINSIC :: iso_fortran_env, ONLY: dp => real64
      IMPORT vanGenuchtenParameters
      IMPLICIT NONE (TYPE, EXTERNAL)

      INTEGER, INTENT(IN)                         :: nz
      REAL(dp), DIMENSION(nz), INTENT(IN)         :: &
        density,                                     &
        grain_radius
      TYPE(vanGenuchtenParameters), INTENT(INOUT) :: vg_params
    END SUBROUTINE van_genuchten_inter
  END INTERFACE

  !-----------------------------------------------------------------------------
  ! parameters
  !-----------------------------------------------------------------------------
  REAL(dp), PARAMETER :: IMP_DENSITY = 830.0

  !-----------------------------------------------------------------------------
  ! routines of this module
  !-----------------------------------------------------------------------------
  CONTAINS


  SUBROUTINE vgAllocateParams( &
  &  nz,                       &
  &  vg_params                 &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN)                         :: nz
    TYPE(vanGenuchtenParameters), INTENT(INOUT) :: vg_params

    !---------------------------------------------------------------------------

    ALLOCATE(vg_params%alpha(nz))
    ALLOCATE(vg_params%n(nz))
    ALLOCATE(vg_params%m(nz))
  END SUBROUTINE vgAllocateParams


  SUBROUTINE vgDeallocateParams(vg_params)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    TYPE(vanGenuchtenParameters), INTENT(INOUT) :: vg_params

    !---------------------------------------------------------------------------

    DEALLOCATE(vg_params%alpha)
    DEALLOCATE(vg_params%n)
    DEALLOCATE(vg_params%m)
  END SUBROUTINE vgDeallocateParams


  SUBROUTINE vgParametersYamaguchi2012( &
  &  nz,                                &
  &  density,                           &
  &  grain_radius,                      &
  &  vg_params                          &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Subroutine: vgParametersYamaguchi2012
    !
    ! van Geuchten parameters for snow according to Yamaguchi et al. 2012.
    !
    ! Yamaguchi, S., Watanabe, K., Katsushima, T., Sato, A., and Kumakura
    ! (2012). Dependence of the water retention curve of snow on snow
    ! characteristics. Annals of Glaciology, 53 (61), pp. 6-12,
    ! https://doi.org/10.3189/2012AoG61A001
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "density", "grain radius", and
    !     "vg_params".
    !   density: Density along the profile (m/s).
    !   grain_radius: Grain radius along the profile (kg/m**3).
    !   vg_params - on input: van Genuchten parameters
    !     (of type vanGenuchtenParameters).
    !
    ! Result:
    !   vg_params - on output: van Genuchten parameters
    !     (of type vanGenuchtenParameters).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      density,                             &
      grain_radius

    TYPE(vanGenuchtenParameters), INTENT(INOUT) :: vg_params

    !---------------------------------------------------------------------------

    ! alpha parameter
    vg_params%alpha = (                    &
    &  4.4E+6_dp                           &
    &  * ((                                &
    &    density / (2.0_dp * grain_radius) &
    &  )**(-0.98_dp))                      &
    )

    ! n parameter
    vg_params%n = (                          &
    &  1.0_dp + (                            &
    &    2.7E-3_dp                           &
    &    * ((                                &
    &      density / (2.0_dp * grain_radius) &
    &    )**(0.61_dp))                       &
    &  )                                     &
    )

    ! m parameter
    vg_params%m = (1.0_dp - (1.0_dp / vg_params%n))
  END SUBROUTINE vgParametersYamaguchi2012


  SUBROUTINE vgParametersDaanen2009( &
  &  nz,                             &
  &  density,                        &
  &  grain_radius,                   &
  &  vg_params                       &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Subroutine: vgParametersDaanen2009
    !
    ! van Geuchten parameters for snow according to Daanen & Nieber 2009.
    !
    ! Daanen, R. P. and Nieber, J. L. (2009). Model for Coupled Liquid Water
    ! Flow and Heat Transport with Phase Change in a Snowpack. Journal of
    ! Cold Regions Engineering, 23 (2), pp. 43-68,
    ! https://doi.org/10.1061/(ASCE)0887-381X(2009)23:2(43)
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "density", "grain radius", and
    !     "vg_params".
    !   density: Density along the profile (m/s).
    !   grain_radius: Grain radius along the profile (kg/m**3).
    !   vg_params - on input: van Genuchten parameters
    !     (of type vanGenuchtenParameters).
    !
    ! Result:
    !   vg_params - on output: van Genuchten parameters
    !     (of type vanGenuchtenParameters).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      density,                             &
      grain_radius

    TYPE(vanGenuchtenParameters), INTENT(INOUT) :: vg_params

    !---------------------------------------------------------------------------

    CALL tfm_essentials_do_nothing(nz, density)

    vg_params%alpha = (                                  &
    &  (30.0_dp * (2.0_dp * (grain_radius * 1000.0_dp))) &
    &  + 12.0_dp                                         &
    )

    vg_params%n = (                                     &
    &  (0.8_dp * (2.0_dp * (grain_radius * 1000.0_dp))) &
    &  + 3.0_dp                                         &
    )

    vg_params%m = (1.0_dp - (1.0_dp / vg_params%n))
  END SUBROUTINE vgParametersDaanen2009


  SUBROUTINE vgParametersYamaguchi2010( &
  &  nz,                                &
  &  density,                           &
  &  grain_radius,                      &
  &  vg_params                          &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Subroutine: vgParametersYamaguchi2010
    !
    ! van Geuchten parameters for snow according to Yamaguchi et al. 2010.
    !
    ! Yamaguchi, S., Katsushima, T., Sato, A. and Kumakura, T. (2010). Water
    ! retention curve of snow with different grain sizes. Cold Regions
    ! Science and Technology, 64 (2), pp. 87-93,
    ! https://doi.org/10.1016/j.coldregions.2010.05.008
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "density", "grain radius", and
    !     "vg_params".
    !   density: Density along the profile (m/s).
    !   grain_radius: Grain radius along the profile (kg/m**3).
    !   vg_params - on input: van Genuchten parameters
    !     (of type vanGenuchtenParameters).
    !
    ! Result:
    !   vg_params - on output: van Genuchten parameters
    !     (of type vanGenuchtenParameters).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      density,                             &
      grain_radius

    TYPE(vanGenuchtenParameters), INTENT(INOUT) :: vg_params

    !---------------------------------------------------------------------------

    CALL tfm_essentials_do_nothing(nz, density)

    vg_params%alpha = (                              &
    &  (7.3_dp * (2.0 * (grain_radius * 1000.0_dp))) &
    &  + 1.9_dp                                      &
    )

    vg_params%n = (                                   &
    &  (-3.3_dp * (2.0 * (grain_radius * 1000.0_dp))) &
    &  + 14.4_dp                                      &
    )

    vg_params%m = (1.0_dp - (1.0_dp / vg_params%n))
  END SUBROUTINE vgParametersYamaguchi2010


  FUNCTION vgDryLayers( &
  &  nz,                &
  &  water_content      &
  ) RESULT(n_water_content)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN)                 :: nz
    REAL(dp), DIMENSION(nz), INTENT(IN) :: water_content
    REAL(dp), DIMENSION(nz)             :: n_water_content

    !---------------------------------------------------------------------------

    n_water_content = max(water_content, 1.0D-6)
  END FUNCTION vgDryLayers


  FUNCTION vgSaturationWC( &
  &  nz,                   &
  &  density               &
  ) RESULT(saturation_wc)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN)                 :: nz
    REAL(dp), DIMENSION(nz), INTENT(IN) :: density
    REAL(dp), DIMENSION(nz)             :: saturation_wc

    !---------------------------------------------------------------------------

    saturation_wc = ((ICE_DENSITY - density) / WATER_DENSITY)
  END FUNCTION vgSaturationWC


  FUNCTION vgSaturationCondCalonne2012( &
  &  nz,                                &
  &  density,                           &
  &  grain_radius                       &
  ) RESULT(saturation_cond)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: vgSaturationCondCalonne2012
    !
    ! Hydraulic conductivity at saturation according to Calonne et al.
    ! (2012).
    !
    ! Calonne, N., Geindreau, C., Flin, F., Morin, S., Lesaffre, B.,
    ! Rolland du Roscoat, S., and Charrier, P. (2012). 3-D image based
    ! numerical computations of snow permeability: links to specific surface
    ! area, density, and microstructual anisotropy. The Cryosphere, 6,
    ! pp. 939-951, https://doi.org/10.5194/tc-6-939-2012
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "density" and "grain_radius".
    !   density: Density along the firn profile (kg/m**3).
    !   grain_raius: Equivalent grain radius along the profile (m).
    !
    ! Result:
    !  saturation_cond: Hydraulic conductivity at saturation (m/s).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      density,                             &
      grain_radius

    REAL(dp), DIMENSION(nz) :: saturation_cond
    REAL(dp), PARAMETER :: DYNAMIC_VISCOSITY = 0.001792_dp

    !---------------------------------------------------------------------------

    saturation_cond = (                                    &
    &  ((WATER_DENSITY * ACC_GRAVITY) / DYNAMIC_VISCOSITY) &
    &  * (                                                 &
    &    (3.0_dp * (grain_radius**2.0_dp))                 &
    &    * exp(-0.013_dp * density)                        &
    &  )                                                   &
    )
  END FUNCTION vgSaturationCondCalonne2012


  FUNCTION vgSaturationCondShimizu1970( &
  &  nz,                                &
  &  density,                           &
  &  grain_radius                       &
  ) RESULT(saturation_cond)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: vgSaturationCondShimizu1970
    !
    ! Hydraulic conductivity at saturation according to Shimizu (1970).
    !
    ! Shimizu, H. (1970). Air Permeability of Deposited Snow. Contributions
    ! from the Institute of Low Temperature Science, A22, pp. 1-32,
    ! http://hdl.handle.net/2115/20234
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "density" and "grain_radius".
    !   density: Density along the firn profile (kg/m**3).
    !   grain_raius: Equivalent grain radius along the profile (m).
    !
    ! Result:
    !  saturation_cond: Hydraulic conductivity at saturation (m/s).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), intent(in) :: &
      density,                             &
      grain_radius

    REAL(dp), DIMENSION(nz) :: saturation_cond
    REAL(dp), PARAMETER :: DYNAMIC_VISCOSITY = 0.001792_dp

    !---------------------------------------------------------------------------

    saturation_cond = (                                      &
    &  0.077_dp * (grain_radius**2.0_dp)                     &
    &  * exp(-0.0078_dp * density)                           &
    &  * ((WATER_DENSITY * ACC_GRAVITY) / DYNAMIC_VISCOSITY) &
    )
  END FUNCTION vgSaturationCondShimizu1970


  FUNCTION vgSaturationSelfConsistent( &
  &  nz,                               &
  &  density,                          &
  &  grain_radius                      &
  ) RESULT(saturation_cond)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), intent(in) :: &
      density,                             &
      grain_radius

    REAL(dp), DIMENSION(nz) :: &
      saturation_cond,         &
      beta

    REAL(dp), PARAMETER :: DYNAMIC_VISCOSITY = 0.001792_dp

    !---------------------------------------------------------------------------

    beta = (1.0_dp - (1.0_dp - (density / ICE_DENSITY)))**(1.0_dp / 3.0_dp)

    saturation_cond = (                                     &
    &  (WATER_DENSITY * ACC_GRAVITY) / DYNAMIC_VISCOSITY    &
    &  * (                                                  &
    &    (grain_radius**2.0_dp / (3.0_dp * (beta**2.0_dp))) &
    &    * (-1.0_dp + (                                     &
    &      (2.0_dp + (3.0_dp * (beta**5.0_dp)))             &
    &      / (beta * (3.0_dp + (2.0_dp * (beta**5.0_dp))))  &
    &    ))                                                 &
    &  )                                                    &
    )
  END FUNCTION vgSaturationSelfConsistent


  FUNCTION vgSaturationCarmanKozeny( &
  &  nz,                             &
  &  density,                        &
  &  grain_radius                    &
  ) RESULT(saturation_cond)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN)                 :: nz
    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      density,                             &
      grain_radius

    REAL(dp), DIMENSION(nz) :: &
      saturation_cond,         &
      porosity

    REAL(dp), PARAMETER :: DYNAMIC_VISCOSITY = 0.001792_dp

    !---------------------------------------------------------------------------

    porosity = (1.0_dp - (density / ICE_DENSITY))

    saturation_cond = (                                      &
    &  ((WATER_DENSITY * ACC_GRAVITY) / DYNAMIC_VISCOSITY)   &
    &  * (                                                   &
    &    ((4.0_dp * grain_radius**2.0) * (porosity**3.0_dp)) &
    &    / (180.0_dp * ((1.0_dp - porosity)**2.0_dp))        &
    &  )                                                     &
    )
  END FUNCTION vgSaturationCarmanKozeny


  FUNCTION vgResidualWCWever2014( &
  &  nz,                          &
  &  water_content,               &
  &  residual_wc                  &
  ) RESULT(n_residual_wc)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function vgResidualWCWever2014
    !
    ! Residual water content as defined by Wever et al. (2014).
    !
    ! Wever, N., Fierz, N., Hirashima, H., and Lehnin, M. (2014). Solving
    ! Richards Equation for snow improves snowpack meltwater runoff
    ! estimations in detailed multi-layer snowpack model. The Cryosphere,
    ! 8, pp. 257-274, https://doi.org/10.5194/tc-8-257-2014
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "water_content" and "residual_wc".
    !   water_content: Volumetric water content along the profile (1).
    !   residual_wc: Volumetric residual water content along the
    !     profile (1).
    !
    ! Result:
    !   n_residual_wc: Volumetric residual water content along the
    !     profile (1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      water_content,                       &
      residual_wc

    REAL(dp), DIMENSION(nz) :: n_residual_wc

    !---------------------------------------------------------------------------

    n_residual_wc = min(                             &
    &  0.02_dp,                                    &
    &  max((0.75_dp * water_content), residual_wc) &
    )
  END FUNCTION vgResidualWCWever2014


  FUNCTION vgRelativeHydraulicCond( &
  &  nz,                            &
  &  head,                          &
  &  vg_params                      &
  ) RESULT(rel_hydraulic_cond)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN)                      :: nz
    REAL(dp), DIMENSION(nz), INTENT(IN)      :: head
    TYPE(vanGenuchtenParameters), INTENT(IN) :: vg_params
    REAL(dp), DIMENSION(nz)                  :: rel_hydraulic_cond

    !---------------------------------------------------------------------------

    rel_hydraulic_cond = (                                           &
    &  ((                                                            &
    &    1.0_dp                                                      &
    &    - (                                                         &
    &      ((vg_params%alpha * abs(head))**(vg_params%n - 1.0_dp))   &
    &      * ((                                                      &
    &        1.0_dp + ((vg_params%alpha * abs(head))**(vg_params%n)) &
    &      )**(-vg_params%m))                                        &
    &    )                                                           &
    &  )**2.0_dp)                                                    &
    &  / ((                                                          &
    &    1.0_dp + ((vg_params%alpha * abs(head))**(vg_params%n))     &
    &  )**(vg_params%m / 2.0_dp))                                    &
    )

    WHERE ( head >= 0.0_dp )
      rel_hydraulic_cond = 1.0_dp
    END WHERE
  END FUNCTION vgRelativeHydraulicCond


  FUNCTION vgSpecificMoistureCapNum( &
  &  nz,                             &
  &  head,                           &
  &  residual_wc,                    &
  &  saturation_wc,                  &
  &  vg_params                       &
  ) RESULT(specific_cap)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), intent(in)    :: &
      head,                                   &
      residual_wc,                            &
      saturation_wc

    TYPE(vanGenuchtenParameters), INTENT(IN) :: vg_params

    REAL(dp), DIMENSION(nz) :: &
      specific_cap,            &
      cplus,                   &
      cminus

    REAL(dp), PARAMETER :: DHEAD = 1.0E-4_dp

    !---------------------------------------------------------------------------

    cplus = (                                                          &
    &  (saturation_wc - residual_wc)                                   &
    &  / ((                                                            &
    &    1.0_dp + ((vg_params%alpha * abs(head + DHEAD))**vg_params%n) &
    &  )**vg_params%m)                                                 &
    )
    cminus = (                                                         &
    &  (saturation_wc - residual_wc)                                   &
    &  / ((                                                            &
    &    1.0_dp + ((vg_params%alpha * abs(head - DHEAD))**vg_params%n) &
    &  )**vg_params%m)                                                 &
    )
    specific_cap = ((cplus - cminus) / (2.0_dp * DHEAD))

    WHERE ( head >= 0.0_dp )
      specific_cap = 1.0_dp
    END WHERE
  END FUNCTION vgSpecificMoistureCapNum


  FUNCTION vgWaterContent( &
  &  nz,                   &
  &  head,                 &
  &  saturation_wc,        &
  &  residual_wc,          &
  &  vg_params             &
  ) RESULT(n_water_content)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      head,                                &
      saturation_wc,                       &
      residual_wc

    TYPE(vanGenuchtenParameters), INTENT(IN) :: vg_params
    REAL(dp), DIMENSION(nz) :: n_water_content

    !---------------------------------------------------------------------------

    n_water_content = (                                          &
    &  residual_wc                                               &
    &  + (                                                       &
    &    (saturation_wc - residual_wc)                           &
    &    / ((                                                    &
    &      1.0_dp + ((vg_params%alpha * abs(head))**vg_params%n) &
    &    )**vg_params%m)                                         &
    &  )                                                         &
    )

    WHERE ( head >= 0.0_dp )
      n_water_content = saturation_wc
    END WHERE
  END FUNCTION vgWaterContent


  FUNCTION vgHead(  &
  &  nz,            &
  &  water_content, &
  &  saturation_wc, &
  &  residual_wc,   &
  &  vg_params      &
  ) RESULT(head)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      water_content,                       &
      saturation_wc,                       &
      residual_wc

    TYPE(vanGenuchtenParameters), INTENT(IN) :: vg_params
    REAL(dp), DIMENSION(nz) :: head

    !---------------------------------------------------------------------------

    head = (                                 &
    &  -(1.0_dp / vg_params%alpha)           &
    &  * ((                                  &
    &    ((                                  &
    &      (saturation_wc - residual_wc)     &
    &      / (water_content - residual_wc)   &
    &    )**(1.0_dp / vg_params%m)) - 1.0_dp &
    &  )**(1.0_dp / vg_params%n))            &
    )
  END FUNCTION vgHead


  FUNCTION vgLiquidMass( &
  &  nz,                 &
  &  depth,              &
  &  water_content       &
  ) RESULT(liquid_mass)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      water_content

    REAL(dp) :: liquid_mass

    !---------------------------------------------------------------------------

    liquid_mass = sum(               &
    &  (depth(2:nz) - depth(1:nz-1)) &
    &  * water_content(2:nz)         &
    &  * WATER_DENSITY               &
    )
  END FUNCTION vgLiquidMass


  SUBROUTINE initKavetski2001( &
  &  nz,                       &
  &  water_content,            &
  &  saturation_wc,            &
  &  residual_wc,              &
  &  saturation_cond,          &
  &  vg_params,                &
  &  dt,                       &
  &  truncerr_tolerance,       &
  &  head,                     &
  &  dhdt,                     &
  &  dwcdt,                    &
  &  local_dt,                 &
  &  safety_inp,               &
  &  eps_inp                   &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), DIMENSION(nz), INTENT(IN)    :: &
      water_content,                          &
      saturation_wc,                          &
      residual_wc,                            &
      saturation_cond

    TYPE(vanGenuchtenParameters), INTENT(IN) :: vg_params
    REAL(dp), INTENT(IN)                     :: dt
    REAL(dp), DIMENSION(2), INTENT(IN)       :: truncerr_tolerance

    REAL(dp), INTENT(IN), OPTIONAL :: &
      safety_inp,                     &
      eps_inp

    REAL(dp), DIMENSION(nz), INTENT(INOUT) :: &
      head,                                   &
      dhdt,                                   &
      dwcdt

    REAL(dp), INTENT(INOUT) :: local_dt

    REAL(dp), DIMENSION(nz) :: &
      specific_cap,            &
      hydraulic_cond

    REAL(dp) :: &
      safety,   &
      eps

    !---------------------------------------------------------------------------

    IF ( present(safety_inp) ) THEN
      safety = safety_inp
    ELSE
      safety = 0.8_dp
    END IF

    IF ( present(eps_inp) ) THEN
      eps = eps_inp
    ELSE
      eps = 1.0E-10_dp
    END IF

    head = vgHead(                  &
    &  nz=nz,                       &
    &  water_content=water_content, &
    &  saturation_wc=saturation_wc, &
    &  residual_wc=residual_wc,     &
    &  vg_params=vg_params          &
    )

    hydraulic_cond = (                                &
    &  saturation_cond                                &
    &  * vgRelativeHydraulicCond(nz, head, vg_params) &
    )

    specific_cap = vgSpecificMoistureCapNum( &
    &  nz=nz,                                &
    &  head=head,                            &
    &  residual_wc=residual_wc,              &
    &  saturation_wc=saturation_wc,          &
    &  vg_params=vg_params                   &
    )

    dhdt = ((-hydraulic_cond * head) / specific_cap)
    dwcdt = (specific_cap * dhdt)

    local_dt = minval(safety * (                       &
    &  (                                               &
    &    ((truncerr_tolerance(2)**0.5_dp) * abs(head)) &
    &    + (truncerr_tolerance(1)**0.5_dp)             &
    &  )                                               &
    &  /(max(abs(dhdt), eps))                          &
    ))
    local_dt = min(dt, local_dt)
  END SUBROUTINE initKavetski2001


  SUBROUTINE timeStepKavetski2001( &
  &  nz,                           &
  &  local_dt,                     &
  &  elapsed_time,                 &
  &  head,                         &
  &  last_head,                    &
  &  last_dhdt,                    &
  &  n_water_content,              &
  &  c_water_content,              &
  &  last_dwcdt,                   &
  &  truncerr_tolerance,           &
  &  safety_inp,                   &
  &  rmin_inp,                     &
  &  rmax_inp,                     &
  &  eps_inp                       &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz

    REAL(dp), INTENT(INOUT) :: &
      local_dt,                &
      elapsed_time

    REAL(dp), DIMENSION(nz), INTENT(INOUT) :: &
      head,                                   &
      last_head,                              &
      last_dhdt,                              &
      n_water_content,                        &
      c_water_content,                        &
      last_dwcdt

    REAL(dp), DIMENSION(2), INTENT(IN) :: truncerr_tolerance

    REAL(dp), INTENT(IN), OPTIONAL :: &
      safety_inp,                     &
      rmin_inp,                       &
      rmax_inp,                       &
      eps_inp

    REAL(dp) :: &
      safety,   &
      rmin,     &
      rmax,     &
      eps,      &
      dt_modifier

    INTEGER :: icrit
    REAL(dp), DIMENSION(nz) :: &
      error,                   &
      condition

    !---------------------------------------------------------------------------

    ! optional arguemnts
    IF ( present(safety_inp) ) THEN
      safety = safety_inp
    ELSE
      safety = 0.8_dp
    END IF

    IF ( present(rmin_inp) ) THEN
      rmin = rmin_inp
    ELSE
      rmin = 0.1_dp
    END IF

    IF ( present(rmax_inp) ) THEN
      rmax = rmax_inp
    ELSE
      rmax = 2.0_dp
    END IF

    IF ( present(eps_inp) ) THEN
      eps = eps_inp
    ELSE
      eps = 1.0E-10_dp
    END IF

    ! time step adjustment
    error = (                                             &
    &  0.5_dp * local_dt                                  &
    &  * abs(last_dhdt - ((head - last_head) / local_dt)) &
    )
    error(1) = 0.0_dp

    condition = (                            &
    &  error                                 &
    &  - (truncerr_tolerance(2) * abs(head)) &
    &  - (truncerr_tolerance(1))             &
    )

    icrit = maxloc(condition, 1)

    dt_modifier = safety * ((                                               &
    &  ((truncerr_tolerance(2) * abs(head(icrit))) + truncerr_tolerance(1)) &
    &  / (max(error(icrit), eps))                                           &
    )**0.5_dp)

    IF ( condition(icrit) < 0.0_dp ) THEN

      last_dwcdt = (                                      &
      &  (                                                &
      &    (2.0_dp * (c_water_content - n_water_content)) &
      &    - (local_dt * last_dwcdt)                      &
      &  )                                                &
      &  / local_dt                                       &
      )
      last_dhdt = ((head - last_head) / local_dt)
      n_water_content = c_water_content
      c_water_content = n_water_content
      last_head = head
      elapsed_time = (elapsed_time + local_dt)
      local_dt = local_dt * min(dt_modifier, rmax)

    ELSE IF ( condition(icrit) >= 0.0_dp ) THEN

      head = last_head
      c_water_content = n_water_content
      local_dt = local_dt * max(dt_modifier, rmin)

    ELSE
      PRINT *, ''
      PRINT *, '********************************************************'
      PRINT *, '* Module: vgRichards                                   *'
      PRINT *, '* Function: timeStepKavetski2011                       *'
      PRINT *, '*                                                      *'
      PRINT *, '* The mixed absolute-relative error, used to determine *'
      PRINT *, '* whether the last time step is accepted or not, seems *'
      PRINT *, '* to show an irregular value!                          *'
      PRINT *, '* Stopping right here!                                 *'
      PRINT *, '********************************************************'
      STOP
    END IF
  END SUBROUTINE timeStepKavetski2001


  FUNCTION solveRichardsEquation( &
  &  nz,                          &
  &  dt,                          &
  &  depth,                       &
  &  head,                        &
  &  n_water_content,             &
  &  c_water_content,             &
  &  last_dwcdt,                  &
  &  hydraulic_cond,              &
  &  specific_cap,                &
  &  influx,                      &
  &  outflux                      &
  ) RESULT(d_head)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      head,                                &
      n_water_content,                     &
      c_water_content,                     &
      last_dwcdt,                          &
      hydraulic_cond,                      &
      specific_cap

    REAL(dp), INTENT(IN) :: &
      influx,               &
      outflux

    REAL(dp), DIMENSION(nz) :: &
      d_head,                  &
      alpha_u,                 &
      alpha_l,                 &
      am,                      &
      au,                      &
      al,                      &
      b

    INTEGER  :: m
    REAL(dp) :: &
      w,        &
      inter_cond

    !---------------------------------------------------------------------------

    alpha_u(2:nz-1) = (                                &
    &  (hydraulic_cond(2:nz-1) + hydraulic_cond(3:nz)) &
    &  / (                                             &
    &    (depth(3:nz) - depth(1:nz-2))                 &
    &    * (depth(3:nz) - depth(2:nz-1))               &
    &  )                                               &
    )
    alpha_u(1)  = -999999.9_dp
    alpha_u(nz) = -999999.9_dp

    alpha_l(2:nz-1) = (                                  &
    &  (hydraulic_cond(2:nz-1) + hydraulic_cond(1:nz-2)) &
    &  / (                                               &
    &    (depth(3:nz) - depth(1:nz-2))                   &
    &    * (depth(2:nz-1) - depth(1:nz-2))               &
    &  )                                                 &
    )
    alpha_l(1)  = -999999.9_dp
    alpha_l(nz) = -999999.9_dp

    ! main diagonal
    am(2:nz-1) = (                            &
    &  (2.0_dp * (specific_cap(2:nz-1)) / dt) &
    &  + alpha_u(2:nz-1)                      &
    &  + alpha_l(2:nz-1)                      &
    )

    ! upper diagonal
    au(2:nz-1) = (-alpha_u(2:nz-1))
    au(1)  = 0.0_dp
    au(nz) = -999999.9_dp

    ! lower diagonal
    al(2:nz-1) = (-alpha_l(2:nz-1))
    al(1)  = -999999.9_dp
    al(nz) = 0.0_dp

    ! right hand side
    b(2:nz-1) = (                                                            &
    &  ((2.0_dp * (n_water_content(2:nz-1) - c_water_content(2:nz-1))) / dt) &
    &  + (last_dwcdt(2:nz-1))                                                &
    &  + (                                                                   &
    &    (hydraulic_cond(3:nz) - hydraulic_cond(1:nz-2))                     &
    &    / (depth(3:nz) - depth(1:nz-2))                                     &
    &  )                                                                     &
    &  - ((alpha_u(2:nz-1) + alpha_l(2:nz-1)) * head(2:nz-1))                &
    &  - (-alpha_u(2:nz-1) * head(3:nz))                                     &
    &  - (-alpha_l(2:nz-1) * head(1:nz-2))                                   &
    )

    ! lower boundary condition (constant flux)
    IF ( outflux >= 0.0_dp ) THEN

      inter_cond = (                                  &
      &  (hydraulic_cond(2) + hydraulic_cond(1))      &
      &  / (2.0_dp * ((depth(2) - depth(1))**2.0_dp)) &
      )

      au(1) = (-inter_cond)

      am(1) = (((2.0_dp * specific_cap(1)) / dt) + inter_cond)

      b(1) = (                                                         &
      &  - (+inter_cond * head(1))                                     &
      &  - (-inter_cond * head(2))                                     &
      &  + ((2.0_dp * (n_water_content(1) - c_water_content(1))) / dt) &
      &  + (last_dwcdt(1))                                             &
      &  + (inter_cond * (depth(2) - depth(1)))                        &
      &  + (outflux / (depth(2) - depth(1)))                           &
      )

    ! free surface
    ELSE IF ( outflux < 0.0_dp ) THEN

      am(1) = 1.0_dp
      au(1) = -1.0_dp
      b(1) = (head(2) - head(1))

    ELSE
      PRINT *, '*************************************************************'
      PRINT *, '* Module: vgRichards                                        *'
      PRINT *, '* Function: solveRichardsEquation                           *'
      PRINT *, '*                                                           *'
      PRINT *, '* The variables "outflux" seems to show an irregular value! *'
      PRINT *, '* Stopping right here!                                      *'
      PRINT *, '*************************************************************'
      STOP
    END IF

    ! upper boundary condition (constant flux)
    IF ( influx >= 0.0_dp ) THEN

      inter_cond = (                                      &
      &  (hydraulic_cond(nz) + hydraulic_cond(nz-1))      &
      &  / (2.0_dp * ((depth(nz) - depth(nz-1))**2.0_dp)) &
      )

      al(nz) = (-inter_cond)

      am(nz) = ((2.0_dp * (specific_cap(nz) / dt)) + inter_cond)

      b(nz) = (                                                          &
      &  - (-inter_cond * head(nz-1))                                    &
      &  - (+inter_cond * head(nz))                                      &
      &  + ((2.0_dp * (n_water_content(nz) - c_water_content(nz))) / dt) &
      &  + (last_dwcdt(nz))                                              &
      &  - (inter_cond * (depth(nz) - depth(nz-1)))                      &
      &  - (-influx / (depth(nz) - depth(nz-1)))                         &
      )

    ! free surface
    ELSE IF ( influx < 0.0_dp ) THEN

      am(nz) = 1.0_dp
      al(nz) = -1.0_dp
      b(nz) = (head(nz-1) - head(1))

    ELSE
      PRINT *, ''
      PRINT *, '***********************************************************'
      PRINT *, '* Module: vgRichards                                      *'
      PRINT *, '* Function: solveRichardsEquation                         *'
      PRINT *, '*                                                         *'
      PRINT *, '* The variable "influx" seems to show an irregular value! *'
      PRINT *, '* Stopping right here!                                    *'
      PRINT *, '***********************************************************'
      STOP
    END IF

    ! TDMA
    DO m = 2, nz, 1
      w = al(m) / am(m-1)
      am(m) = am(m) - (w * au(m-1))
      b(m)  = b(m)  - (w * b(m-1))
    END DO

    d_head(nz) = (b(nz) / am(nz))

    DO m = (nz - 1), 1, -1
      d_head(m) = (b(m) - (au(m) * d_head(m+1))) / am(m)
    END DO
  END FUNCTION solveRichardsEquation


  FUNCTION vgRichardsAdvanceTimeStep( &
  &  nz,                              &
  &  dt,                              &
  &  depth,                           &
  &  density,                         &
  &  grain_radius,                    &
  &  water_content,                   &
  &  liquid_accumulation,             &
  &  van_genuchten_model              &
  ) RESULT(n_water_content)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      density,                             &
      grain_radius,                        &
      water_content

    REAL(dp), INTENT(IN) :: liquid_accumulation

    PROCEDURE(van_genuchten_inter), POINTER :: van_genuchten_model

    REAL(dp), DIMENSION(nz) :: n_water_content

    INTEGER  :: iter
    REAL(dp) ::     &
      backsteps,    &
      local_dt,     &
      elapsed_time, &
      influx,       &
      outflux

    REAL(dp), dimension(2)  :: &
      picard_tolerance,        &
      step_tolerance

    TYPE(vanGenuchtenParameters) :: vg_params

    REAL(dp), DIMENSION(nz) :: &
      last_dhdt,               &
      residuum,                &
      saturation_wc,           &
      residual_wc,             &
      saturation_cond,         &
      rel_hydraulic_cond,      &
      hydraulic_cond,          &
      head,                    &
      last_head,               &
      d_head,                  &
      c_water_content,         &
      specific_cap,            &
      last_dwcdt,              &
      dry_layers,              &
      eff_saturation

    !---------------------------------------------------------------------------

    eff_saturation = 1.0E-3_dp
    saturation_wc = 0.9_dp * (1.0_dp - (density / ICE_DENSITY))

    WHERE ( water_content == 0.0_dp )
      dry_layers = (                                             &
      &  (eff_saturation * saturation_wc)                        &
      &  / (1.0_dp - 0.75_dp + (0.75_dp * eff_saturation)) &
      )
    ELSE WHERE
      dry_layers = 0.0_dp
    END WHERE

    ! defintions
    step_tolerance = [ &
    &  0.0_dp,         &
    &  1.0E-2_dp       &
    ]
    picard_tolerance = [            &
    &  0.01_dp * step_tolerance(1), &
    &  0.01_dp * step_tolerance(2)  &
    ]

    backsteps = 0.0_dp

    n_water_content = water_content + dry_layers
    influx = liquid_accumulation
    outflux = -1

    ! van Genuchten parameters
    CALL vgAllocateParams( &
    &  nz=nz,              &
    &  vg_params=vg_params &
    )
    CALL van_genuchten_model(     &
    &  nz=nz,                     &
    &  density=density,           &
    &  grain_radius=grain_radius, &
    &  vg_params=vg_params        &
    )

    ! props
    residual_wc = min(0.02_dp, (0.75_dp * n_water_content), saturation_wc)

    saturation_cond = vgSaturationCondShimizu1970( &
    &  nz=nz,                                      &
    &  density=density,                            &
    &  grain_radius=grain_radius                   &
    )
    !saturation_cond = vgSaturationCondCalonne2012(nz, density, grain_radius)
    !saturation_cond = vgSaturationSelfConsistent(nz, density, grain_radius)
    !saturation_cond = vgSaturationCarmanKozeny(nz, density, grain_radius)

    CALL initKavetski2001(                &
    &  nz=nz,                             &
    &  water_content=n_water_content,     &
    &  saturation_wc=saturation_wc,       &
    &  residual_wc=residual_wc,           &
    &  saturation_cond=saturation_cond,   &
    &  vg_params=vg_params,               &
    &  dt=dt,                             &
    &  truncerr_tolerance=step_tolerance, &
    &  head=head,                         &
    &  dhdt=last_dhdt,                    &
    &  dwcdt=last_dwcdt,                  &
    &  local_dt=local_dt                  &
    )
    local_dt = 1.0E-6_dp
    last_head = head

    ! time loop
    elapsed_time = 0.0_dp
    DO WHILE ( elapsed_time /= dt )

      residuum = 999999.9_dp
      iter = 0

      ! Picard loop
      DO WHILE ( (maxval(residuum) >= 0.0_dp) .AND. (iter < 100) )

        rel_hydraulic_cond = vgRelativeHydraulicCond( &
        &  nz=nz,                                     &
        &  head=head,                                 &
        &  vg_params=vg_params                        &
        )
        hydraulic_cond = (rel_hydraulic_cond * saturation_cond)

        specific_cap = vgSpecificMoistureCapNum( &
        &  nz=nz,                                &
        &  head=head,                            &
        &  residual_wc=residual_wc,              &
        &  saturation_wc=saturation_wc,          &
        &  vg_params=vg_params                   &
        )

        d_head = solveRichardsEquation(     &
        &  nz=nz,                           &
        &  dt=local_dt,                     &
        &  depth=depth,                     &
        &  head=head,                       &
        &  n_water_content=n_water_content, &
        &  c_water_content=c_water_content, &
        &  last_dwcdt=last_dwcdt,           &
        &  hydraulic_cond=hydraulic_cond,   &
        &  specific_cap=specific_cap,       &
        &  influx=influx,                   &
        &  outflux=outflux                  &
        )

        head = (head + d_head)

        c_water_content = vgWaterContent( &
        &  nz=nz,                         &
        &  head=head,                     &
        &  saturation_wc=saturation_wc,   &
        &  residual_wc=residual_wc,       &
        &  vg_params=vg_params            &
        )

        ! exception head
        IF ( ( any(isnan(head)) ) .OR. ( any(abs(head) > huge(abs(head))) ) ) THEN
          PRINT *, ''
          PRINT *, '**************************************************'
          PRINT *, '* Module: tfm_liquid                             *'
          PRINT *, '* Function: vgRichardsAdvanceTimeStep            *'
          PRINT *, '*                                                *'
          PRINT *, '* The variable "head" shows one or more NaN or   *'
          PRINT *, '* Inifinity values! This should not be the case! *'
          PRINT *, '* Stopping right here!                           *'
          PRINT *, '**************************************************'
          STOP
        END IF

        iter = (iter + 1)
        residuum = (                           &
        &  abs(d_head)                         &
        &  - (picard_tolerance(2) * abs(head)) &
        &  - (picard_tolerance(1))             &
        )
      END DO

      ! time and time step control
      CALL timeStepKavetski2001(           &
      &  nz=nz,                            &
      &  local_dt=local_dt,                &
      &  elapsed_time=elapsed_time,        &
      &  head=head,                        &
      &  last_head=last_head,              &
      &  last_dhdt=last_dhdt,              &
      &  n_water_content=n_water_content,  &
      &  c_water_content=c_water_content,  &
      &  last_dwcdt=last_dwcdt,            &
      &  truncerr_tolerance=step_tolerance &
      )

      IF ( (elapsed_time + local_dt) > dt ) THEN
        local_dt = (dt - elapsed_time)
      END IF
    END DO

    n_water_content = (n_water_content - dry_layers)
    WHERE ( (n_water_content <= 1.0D-4) )
      n_water_content = 0.0_dp
    END WHERE

    ! deallocation of the van Genuchten parameters
    CALL vgDeallocateParams(vg_params)
  END FUNCTION vgRichardsAdvanceTimeStep


  SUBROUTINE tfm_liquid_RichardsEquation( &
  &  nz,                                  &
  &  dt,                                  &
  &  depth,                               &
  &  density,                             &
  &  temperature,                         &
  &  grain_radius,                        &
  &  water_content,                       &
  &  liquid_accumulation,                 &
  &  runoff,                              &
  &  van_genuchten_model                  &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt
    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      grain_radius

    REAL(dp), INTENT(IN) :: liquid_accumulation

    PROCEDURE(van_genuchten_inter), POINTER :: van_genuchten_model

    REAL(dp), DIMENSION(nz), INTENT(INOUT) :: &
      density,                                &
      temperature,                            &
      water_content

    REAL(dp), intent(inout) :: runoff

    REAL(dp), DIMENSION(nz) :: &
      dz,                      &
      water_mass,              &
      refreeze_cap,            &
      ice_cap,                 &
      storage,                 &
      d_temperature,           &
      d_density

    !---------------------------------------------------------------------------

    IF ( (liquid_accumulation > 0.0_dp) .OR. any(water_content > 0.0_dp) ) THEN
      water_content = vgRichardsAdvanceTimeStep(  &
      &  nz=nz,                                   &
      &  dt=dt,                                   &
      &  depth=depth,                             &
      &  density=density,                         &
      &  grain_radius=grain_radius,               &
      &  water_content=water_content,             &
      &  liquid_accumulation=liquid_accumulation, &
      &  van_genuchten_model=van_genuchten_model  &
      )
    END IF

    dz(1) = 1.0_dp
    dz(2:nz) = (depth(2:nz) - depth(1:nz-1))

    water_mass = (water_content * WATER_DENSITY * dz)

    IF ( any(water_mass > 0.0_dp) ) THEN

      ! potential mass per square meter that might refreeze
      refreeze_cap = (                                                  &
      &  (SPECIFIC_HEAT_ICE * density * dz * (MELT_TEMP - temperature)) &
      &  / (LATENT_HEAT)                                                &
      )

      ! pore space in kg ice equivalent per square meter
      ice_cap = ((1.0_dp - (density / ICE_DENSITY)) * ICE_DENSITY * dz)

      storage = min(refreeze_cap, ice_cap, water_mass)

      ! temperature change
      d_temperature = (LATENT_HEAT / (SPECIFIC_HEAT_ICE * density * dz)) * storage
      !d_temperature(1) = 0.0_dp
      temperature = (temperature + d_temperature)

      ! density change
      d_density = (storage / dz)
      d_density(1) = 0.0_dp
      density = (density + d_density)

      water_content = water_content - (storage / WATER_DENSITY / dz)
      WHERE ( water_content < 1.0E-10_dp )
        water_content = 0.0_dp
      END WHERE
    END IF
  END SUBROUTINE tfm_liquid_RichardsEquation


  SUBROUTINE tfm_liquid_bucket( &
  &  nz,                        &
  &  dt,                        &
  &  depth,                     &
  &  density,                   &
  &  temperature,               &
  &  grain_radius,              &
  &  liquid_water,              &
  &  infiltration_rate,         &
  &  runoff,                    &
  &  van_genuchten_model        &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      grain_radius

    REAL(dp), INTENT(IN) :: infiltration_rate

    PROCEDURE(van_genuchten_inter), POINTER :: van_genuchten_model

    REAL(dp), DIMENSION(nz), INTENT(INOUT) :: &
      density,                                &
      temperature,                            &
      liquid_water

    REAL(dp), INTENT(INOUT) :: runoff

    INTEGER :: n
    REAL(dp) ::          &
      water,             &
      dz,                &
      irr_water_content, &
      ice_cap,           &
      refreeze_cap,      &
      imp_cap,           &
      storage

    !---------------------------------------------------------------------------

    CALL tfm_essentials_do_nothing(nz, grain_radius)

    ! infiltrated water mass
    IF (infiltration_rate <= 0.0_dp) RETURN
    water = dt * infiltration_rate * WATER_DENSITY

    DO n = nz - 1, 1, -1

      dz = (depth(n+1) - depth(n))

      ! potential mass per square meter that might be frozen
      refreeze_cap = (                                                        &
      &  (SPECIFIC_HEAT_ICE * density(n) * dz * (MELT_TEMP - temperature(n))) &
      &  / LATENT_HEAT                                                        &
      )

      ! pore space in kg ice equivalent per square meter
      ice_cap = (1.0_dp - (density(n) / ICE_DENSITY)) * ICE_DENSITY * dz

      ! pore space available until a impermeable layer is formed
      imp_cap = max(0.0_dp, (IMP_DENSITY - density(n)) * dz)

      ! maximum amount of water to be refrozen
      storage = min(refreeze_cap, ice_cap, water, imp_cap)

      IF ( storage < 0.0_dp ) THEN
        PRINT *, 'module: tfm_liquid                                   '
        PRINT *, 'subroutine: tfm_liquid_bucket                        '
        PRINT *, 'Variable storage became negative.                    '
        PRINT *, 'strorage: ', storage, 'kg m-2'
        PRINT *, 'Phyiscally this is not possible. Most certainly there'
        PRINT *, 'are temperatures above the melting point occuring.   '
        PRINT *, 'max(temperature): ', maxval(temperature), 'K'
        PRINT *, 'Stopping the simulation right here.                  '
        STOP
      END IF

      ! temperature change due to refreezing
      temperature(n) = (                                                   &
      &  temperature(n)                                                    &
      &  + (LATENT_HEAT / (SPECIFIC_HEAT_ICE * density(n) * dz)) * storage &
      )

      ! remaining water
      water = water - storage

      ! density change
      density(n) = density(n) + (storage / dz)
      IF ( density(n) >= IMP_DENSITY .OR. water <= 0.0_dp ) EXIT

      ! irreducable water content
      CALL tfm_liquid_Coleou1998(            &
      &  density=density(n),                 &
      &  irr_water_content=irr_water_content &
      )

      ! liquid water remaining in the layer
      storage = min(water, (irr_water_content * WATER_DENSITY * dz))
      liquid_water(n) = liquid_water(n) + storage

      ! remaining water
      water = water - storage
      IF ( water <= 0.0_dp ) EXIT
    END DO

    runoff = water ! total runoff water kg
  END SUBROUTINE tfm_liquid_bucket


  SUBROUTINE tfm_liquid_Coleou1998( &
  &  density,                       &
  &  irr_water_content              &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)

    REAL(dp), INTENT(IN)    :: density
    REAL(dp), INTENT(INOUT) :: irr_water_content

    irr_water_content = 0.017_dp + 0.057_dp * ((ICE_DENSITY - density) / density)
  END SUBROUTINE tfm_liquid_Coleou1998
END MODULE tfm_liquid

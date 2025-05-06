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
MODULE tfm_grain
!-----------------------------------------------------------------------
! Module: tfm_grain
!
! Dependencies: tfm_essentials, tfm_constants
!
! Functions:
!  tfm_grain_arthern2010: Arthern et al. (2010)
!  tfm_grain_zwally2002: Zwally & Li (2002)
!  tfm_grain_brun1989: Brun (1989)
!  tfm_grain_tusima1978: Tusima (1978)
!  tfm_grain_katsushima2009: Katsushima et al. (2009)
!-----------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! modules used
  !-----------------------------------------------------------------------------
  USE, INTRINSIC :: iso_fortran_env, ONLY: dp => real64

  USE tfm_essentials, ONLY : &
    tfm_essentials_do_nothing

  USE tfm_constants, ONLY : &
    WATER_DENSITY,          &
    GAS_CONST,              &
    PI,                     &
    SECONDS_YEAR

  !-----------------------------------------------------------------------------
  ! declarations
  !-----------------------------------------------------------------------------
  IMPLICIT NONE (TYPE, EXTERNAL)
  PRIVATE

  !-----------------------------------------------------------------------------
  ! public routines
  !-----------------------------------------------------------------------------
  PUBLIC ::                &
    tfm_grain_arthern2010, &
    tfm_grain_zwally2002,  &
    tfm_grain_brun1989,    &
    tfm_grain_tusima1978,  &
    tfm_grain_katsushima2009

  !-----------------------------------------------------------------------------
  ! public routines of this module
  !-----------------------------------------------------------------------------
  CONTAINS


  SUBROUTINE tfm_grain_Area2RadiusChangeExplicit( &
  &  nz,                                          &
  &  dt,                                          &
  &  grain_radius,                                &
  &  area_growth_rate,                            &
  &  d_grain_radius                               &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      grain_radius,                        &
      area_growth_rate

    REAL(dp), DIMENSION(nz), INTENT(INOUT) :: d_grain_radius

    !---------------------------------------------------------------------------

    d_grain_radius = (                &
    &  ((2.0_dp * PI)**(-1.0_dp)) &
    &  * (grain_radius**(-1.0_dp))  &
    &  * area_growth_rate             &
    &  * dt                           &
    )
  END SUBROUTINE tfm_grain_Area2RadiusChangeExplicit


  SUBROUTINE tfm_grain_Area2RadiusChangeImplicit( &
  &  nz,                                          &
  &  dt,                                          &
  &  grain_radius,                                &
  &  area_growth_rate,                            &
  &  d_grain_radius                               &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      grain_radius,                        &
      area_growth_rate

    REAL(dp), DIMENSION(nz), INTENT(INOUT) :: d_grain_radius

    !---------------------------------------------------------------------------

    d_grain_radius = ((                    &
    &  (grain_radius / 2.0_dp)             &
    &  + ((                                &
    &    ((grain_radius**2.0_dp) / 4.0_dp) &
    &    + (                               &
    &      ((2.0_dp * PI)**(-1.0_dp))      &
    &      * area_growth_rate              &
    &      * dt                            &
    &    )                                 &
    &  )**0.5_dp)                          &
    ) - grain_radius)
  END SUBROUTINE tfm_grain_Area2RadiusChangeImplicit


  SUBROUTINE tfm_grain_Volume2RadiusChangeExplicit( &
  &  nz,                                            &
  &  dt,                                            &
  &  grain_radius,                                  &
  &  volume_growth_rate,                            &
  &  d_grain_radius                                 &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      grain_radius,                        &
      volume_growth_rate

    REAL(dp), DIMENSION(nz), INTENT(INOUT) :: d_grain_radius

    !---------------------------------------------------------------------------

    d_grain_radius = (             &
    &  ((4.0_dp * PI)**(-1.0_dp))  &
    &  * (grain_radius**(-2.0_dp)) &
    &  * volume_growth_rate        &
    &  * dt                        &
    )
  END SUBROUTINE tfm_grain_Volume2RadiusChangeExplicit


  SUBROUTINE tfm_grain_Volume2RadiusChangeImplicit( &
  &  nz,                                            &
  &  dt,                                            &
  &  grain_radius,                                  &
  &  volume_growth_rate,                            &
  &  d_grain_radius                                 &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt

    REAL(dp), DIMENSION(nz), intent(in) :: &
      grain_radius,                        &
      volume_growth_rate

    REAL(dp), DIMENSION(nz), INTENT(INOUT) :: d_grain_radius

    real(dp), dimension(nz) :: &
      p,                       &
      q,                       &
      delta,                   &
      u,                       &
      v

    !---------------------------------------------------------------------------

    p = ((-1.0_dp / 3.0_dp) * (-grain_radius**2.0_dp))
    q = (                                                      &
    &  ((2.0_dp / 27.0_dp) * (-grain_radius**3.0))             &
    &  + ((-1.0_dp / (4.0_dp * PI)) * volume_growth_rate * dt) &
    )

    delta = (((q**2.0_dp) / 4.0_dp) + ((p**3.0_dp) / 27.0_dp))

    u = (((-q / 2.0_dp) + (delta**0.5_dp))**(1.0_dp / 3.0_dp))
    v = (((-q / 2.0_dp) - (delta**0.5_dp))**(1.0_dp / 3.0_dp))

    d_grain_radius = ((grain_radius / 3.0_dp) + u + v - grain_radius)
  END SUBROUTINE tfm_grain_Volume2RadiusChangeImplicit


  SUBROUTINE tfm_grain_Area2RadiusChange( &
  &  nz,                                  &
  &  dt,                                  &
  &  grain_radius,                        &
  &  area_growth_rate,                    &
  &  d_grain_radius,                      &
  &  solving_method                       &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt

    REAL(dp), DIMENSION(nz), INTENT(IN)    :: &
      grain_radius,                           &
      area_growth_rate

    REAL(dp), DIMENSION(nz), INTENT(INOUT) :: d_grain_radius
    CHARACTER(LEN=*), INTENT(IN)           :: solving_method

    !---------------------------------------------------------------------------

    IF ( solving_method == 'explicit' ) THEN
      CALL tfm_grain_Area2RadiusChangeExplicit( &
      &  nz=nz,                                 &
      &  dt=dt,                                 &
      &  grain_radius=grain_radius,             &
      &  area_growth_rate=area_growth_rate,     &
      &  d_grain_radius=d_grain_radius          &
      )

    ELSE IF ( solving_method == 'implicit' ) THEN
      CALL tfm_grain_Area2RadiusChangeImplicit( &
      &  nz=nz,                                 &
      &  dt=dt,                                 &
      &  grain_radius=grain_radius,             &
      &  area_growth_rate=area_growth_rate,     &
      &  d_grain_radius=d_grain_radius          &
      )

    ELSE
      PRINT *, '*******************************************************'
      PRINT *, '* Module: tfm_grain                                   *'
      PRINT *, '* Subroutine: tfm_grain_Area2RadiusChange             *'
      PRINT *, '*                                                     *'
      PRINT *, '* Argument "solving_method" has eiter to be           *'
      PRINT *, '* "explicit" or "implicit". This seems not to be the  *'
      PRINT *, '* case.                                               *'
      PRINT *, '* Stopping right here!                                *'
      PRINT *, '*******************************************************'
      STOP
    END IF
  END SUBROUTINE tfm_grain_Area2RadiusChange


  SUBROUTINE tfm_grain_Volume2RadiusChange( &
  &  nz,                                    &
  &  dt,                                    &
  &  grain_radius,                          &
  &  volume_growth_rate,                    &
  &  d_grain_radius,                        &
  &  solving_method                         &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      grain_radius,                        &
      volume_growth_rate

    REAL(dp), DIMENSION(nz), INTENT(INOUT) :: d_grain_radius
    CHARACTER(LEN=*), INTENT(IN)           :: solving_method

    !---------------------------------------------------------------------------

    IF ( solving_method == 'explicit' ) THEN
      CALL tfm_grain_Volume2RadiusChangeExplicit( &
      &  nz=nz,                                   &
      &  dt=dt,                                   &
      &  grain_radius=grain_radius,               &
      &  volume_growth_rate=volume_growth_rate,   &
      &  d_grain_radius=d_grain_radius            &
      )

    ELSE IF ( solving_method == 'implicit' ) THEN
      CALL tfm_grain_Volume2RadiusChangeImplicit( &
      &  nz=nz,                                   &
      &  dt=dt,                                   &
      &  grain_radius=grain_radius,               &
      &  volume_growth_rate=volume_growth_rate,   &
      &  d_grain_radius=d_grain_radius            &
      )

    ELSE
      PRINT *, '*******************************************************'
      PRINT *, '* Module: tfm_grain                                   *'
      PRINT *, '* Subroutine: tfm_grain_Volume2RadiusChange           *'
      PRINT *, '*                                                     *'
      PRINT *, '* Argument "solving_method" has eiter to be           *'
      PRINT *, '* "explicit" or "implicit". This seems not to be the  *'
      PRINT *, '* case.                                               *'
      PRINT *, '* Stopping right here!                                *'
      PRINT *, '*******************************************************'
      STOP
    END IF
  END SUBROUTINE tfm_grain_Volume2RadiusChange


  FUNCTION tfm_grain_arthern2010( &
  &  nz,                          &
  &  dt,                          &
  &  temperature,                 &
  &  density,                     &
  &  liquid_water,                &
  &  grain_radius                 &
  ) RESULT(d_grain_radius)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_grain_arthern2010
    !
    ! Arthern, R. J., Vaughan, D. G., Rankin, A. M., Mulvaney, R., and
    ! Thomas, E. R. In situ measurements of Antarctic snow compaction
    ! compared with predicitions of models. Journal of Geophysical Research:
    ! Earth Surface, 115 (F3), (2010). https://doi.org/10.1029/2009JF001306
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of the variables "temperature", "density", and
    !       "d_grain_radius".
    !   dt: Time step (s).
    !   temperature: Temperature along the firn profile (K).
    !   density: Density along the firn profile (kg / m**3).
    !   liquid_water: Volumetric liquid water content along the profile (1).
    !   grain_radius: Grain radius alon ght eprofile (m).
    !
    ! Result:
    !   d_grain_radius: Grain radius change along the firn profile (m).
    !---------------------------------------------------------------------------

    ! Parameters from Arthern & Wingham (2010) following Paterson (1994)
    REAL(dp), PARAMETER ::    &
      PRE_FACTOR = 1.3e-7_dp, & ! m**2 s**-1
      ACTIVATION_ENERGY = 42400.0_dp !J mol**-1

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      temperature,                         &
      density,                             &
      liquid_water,                        &
      grain_radius

    REAL(dp), DIMENSION(nz) :: &
      d_grain_radius,          &
      area_growth_rate

    !---------------------------------------------------------------------------

    CALL tfm_essentials_do_nothing(nz, density)
    CALL tfm_essentials_do_nothing(nz, liquid_water)

    ! growth rate (m**2 s**-1)
    area_growth_rate = (                                     &
    &  PRE_FACTOR                                            &
    &  * exp(-ACTIVATION_ENERGY / (GAS_CONST * temperature)) &
    )

    CALL tfm_grain_Area2RadiusChange(     &
    &  nz=nz,                             &
    &  dt=dt,                             &
    &  grain_radius=grain_radius,         &
    &  area_growth_rate=area_growth_rate, &
    &  d_grain_radius=d_grain_radius,     &
    &  solving_method='implicit'          &
    )
  END FUNCTION tfm_grain_arthern2010


  FUNCTION tfm_grain_zwally2002( &
  &  nz,                         &
  &  dt,                         &
  &  temperature,                &
  &  density,                    &
  &  liquid_water,               &
  &  grain_radius                &
  ) RESULT(d_grain_radius)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function: tfm_grain_li2002
    !
    ! Grain growth following Zwally & Li (2002) as described in
    ! Li & Zwally (2011).
    !
    ! Zwally, H. J. and Li, J. Seasonal and interannual variations of firn
    ! densification and ice-sheet elevation at the Greenland summit.
    ! Journal of Glaciology, 48 (171), (2002).
    ! https://doi.org/10.3189/172756502781831403
    !
    ! Li, J. and Zwally, H. J. Modeling of firn compaction for estimating
    ! ice-sheet mass change from observed ice-sheet elevation change. Annals
    ! of Glaciology, 52 (59), (2011).
    ! https://doi.org/10.3189/172756411799096321
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of the variables "temperature", "density", and
    !       "d_grain_radius".
    !   dt: Time step (s).
    !   temperature: Temperature along the firn profile (K).
    !   density: Density along the firn profile (kg / m**3).
    !   liquid_water: Volumetric liquid water content along the profile (1).
    !   grain_radius: Grain radius alon ght eprofile (m).
    !
    ! Result:
    !   d_grain_radius: Grain radius change along the firn profile (m).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      temperature,                         &
      density,                             &
      liquid_water,                        &
      grain_radius

    REAL(dp), DIMENSION(nz) :: &
      area_growth_rate,        &
      d_grain_radius

    !---------------------------------------------------------------------------

    CALL tfm_essentials_do_nothing(nz, density)
    CALL tfm_essentials_do_nothing(nz, liquid_water)

    ! growth rate (mm**2 a**-1)
    area_growth_rate = 8.36_dp * ((273.2_dp - temperature)**(-2.061_dp))
    area_growth_rate = area_growth_rate * 1.0E-6_dp / SECONDS_YEAR ! (m**2 s**-1)

    CALL tfm_grain_Area2RadiusChange(     &
    &  nz=nz,                             &
    &  dt=dt,                             &
    &  grain_radius=grain_radius,         &
    &  area_growth_rate=area_growth_rate, &
    &  d_grain_radius=d_grain_radius,     &
    &  solving_method='implicit'          &
    )
  END FUNCTION tfm_grain_zwally2002


  FUNCTION tfm_grain_brun1989( &
  &  nz,                       &
  &  dt,                       &
  &  temperature,              &
  &  density,                  &
  &  liquid_water,             &
  &  grain_radius              &
  ) RESULT(d_grain_radius)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Fucntion: tfm_grain_brun1989
    !
    ! Grain growth following Brun 1989.
    !
    ! Brun, E. (1989). Investigation On Wet-Snow Metamorphism in Respect of
    ! Liquid-Water Content. Annals of Glaciology, 13, pp. 22-26,
    ! https://doi.org/10.3189/S0260305500007576
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of the variables "temperature", "density", and
    !       "d_grain_radius".
    !   dt: Time step (s).
    !   temperature: Temperature along the firn profile (K).
    !   density: Density along the firn profile (kg / m**3).
    !   liquid_water: Volumetric liquid water content along the profile (1).
    !   grain_radius: Grain radius alon ght eprofile (m).
    !
    ! Result:
    !   d_grain_radius: Grain radius change along the firn profile (m).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      temperature,                         &
      density,                             &
      liquid_water,                        &
      grain_radius

    REAL(dp), DIMENSION(nz) :: &
      d_grain_radius,          &
      gravimetric_wc,          &
      volume_growth_rate

    !---------------------------------------------------------------------------

    CALL tfm_essentials_do_nothing(nz, temperature)

    gravimetric_wc = (100.0_dp * (                  &
    &  (liquid_water * WATER_DENSITY)               &
    &  / (density + (liquid_water * WATER_DENSITY)) &
    ))

    ! growth rate (mm**3 s**-1)
    volume_growth_rate = (1.28E-8_dp + (4.22E-10_dp * (gravimetric_wc**3.0_dp)))
    volume_growth_rate = (volume_growth_rate * 1.0E-9_dp) ! (m**3 s**-1)

    CALL tfm_grain_Volume2RadiusChange(       &
    &  nz=nz,                                 &
    &  dt=dt,                                 &
    &  grain_radius=grain_radius,             &
    &  volume_growth_rate=volume_growth_rate, &
    &  d_grain_radius=d_grain_radius,         &
    &  solving_method='explicit'              &
    )
  END FUNCTION tfm_grain_brun1989


  FUNCTION tfm_grain_tusima1978( &
  &  nz,                         &
  &  dt,                         &
  &  temperature,                &
  &  density,                    &
  &  liquid_water,               &
  &  grain_radius                &
  ) RESULT(d_grain_radius)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function tfm_grain_tusima1978
    !
    ! Grain growth following Tusima (1978) as described by Katsushima et al.
    ! (2009). The original paper is written in Japanese.
    !
    ! Tusima, K. (1978). Grain coarsening of ice particles immersed in pure
    ! water. Seppyo, 40 (4), pp. 155-165.
    !
    ! Katsushima, T., Kamakura, T., and Takeuchi, Y. (2009). A multiple
    ! layer model including a parametrization of vertical water channel
    ! process in snowpack. Cold Regions Science and Technology, 59, pp.
    ! 143-151, https://doi.org/10.1016/j.coldregions.2009.09.002
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of the variables "temperature", "density", and
    !       "d_grain_radius".
    !   dt: Time step (s).
    !   temperature: Temperature along the firn profile (K).
    !   density: Density along the firn profile (kg / m**3).
    !   liquid_water: Volumetric liquid water content along the profile (1).
    !   grain_radius: Grain radius alon ght eprofile (m).
    !
    ! Result:
    !   d_grain_radius: Grain radius change along the firn profile (m).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      temperature,                         &
      density,                             &
      liquid_water,                        &
      grain_radius

    REAL(dp), dimension(nz) :: &
      d_grain_radius,          &
      radius_growth_rate

    !---------------------------------------------------------------------------

    CALL tfm_essentials_do_nothing(nz, temperature)
    CALL tfm_essentials_do_nothing(nz, density)
    CALL tfm_essentials_do_nothing(nz, liquid_water)

    ! growth rate (mm s**-1)
    radius_growth_rate = (                               &
    &  (2.5E-4_dp / ((grain_radius * 1.0E3_dp)**2.0_dp)) &
    &  / 3600.0_dp                                       &
    )
    radius_growth_rate = (radius_growth_rate * 1.0E-3_dp) ! (m s**-1)

    d_grain_radius = (dt * radius_growth_rate)
  END FUNCTION tfm_grain_tusima1978


  FUNCTION tfm_grain_katsushima2009( &
  &  nz,                             &
  &  dt,                             &
  &  temperature,                    &
  &  density,                        &
  &  liquid_water,                   &
  &  grain_radius                    &
  ) RESULT(d_grain_radius)
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Function tfm_grain_Katushima2009
    !
    ! Grain growth following Katsushima et al. (2009) combining the
    ! parametrizations of Brun (1989) and Tusima (1978) depending on the
    ! liquid water content of the firn.
    !
    ! Katsushima, T., Kamakura, T., and Takeuchi, Y. (2009). A multiple
    ! layer model including a parametrization of vertical water channel
    ! process in snowpack. Cold Regions Science and Technology, 59, pp.
    ! 143-151, https://doi.org/10.1016/j.coldregions.2009.09.002
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of the variables "temperature", "density", and
    !       "d_grain_radius".
    !   dt: Time step (s).
    !   temperature: Temperature along the firn profile (K).
    !   density: Density along the firn profile (kg / m**3).
    !   liquid_water: Volumetric liquid water content along the profile (1).
    !   grain_radius: Grain radius alon ght eprofile (m).
    !
    ! Result:
    !   d_grain_radius: Grain radius change along the firn profile (m).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), INTENT(IN) :: dt

    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      temperature,                         &
      density,                             &
      liquid_water,                        &
      grain_radius

    REAL(dp), DIMENSION(nz) :: &
      gravimetric_wc,          &
      d_grain_radius_brun,     &
      d_grain_radius_tusima,   &
      d_grain_radius

    !---------------------------------------------------------------------------

    CALL tfm_essentials_do_nothing(nz, temperature)

    gravimetric_wc = (100.0_dp * (                  &
    &  (liquid_water * WATER_DENSITY)               &
    &  / (density + (liquid_water * WATER_DENSITY)) &
    ))

    d_grain_radius_brun = tfm_grain_brun1989( &
    &  nz=nz,                                 &
    &  dt=dt,                                 &
    &  temperature=temperature,               &
    &  density=density,                       &
    &  liquid_water=liquid_water,             &
    &  grain_radius=grain_radius              &
    )

    d_grain_radius_tusima = tfm_grain_tusima1978( &
    &  nz=nz,                                     &
    &  dt=dt,                                     &
    &  temperature=temperature,                   &
    &  density=density,                           &
    &  liquid_water=liquid_water,                 &
    &  grain_radius=grain_radius                  &
    )

    WHERE ( gravimetric_wc <= 10.0_dp )
      d_grain_radius = d_grain_radius_brun
    ELSE WHERE
      d_grain_radius = min(d_grain_radius_brun, d_grain_radius_tusima)
    END WHERE
  END FUNCTION tfm_grain_katsushima2009
END MODULE tfm_grain

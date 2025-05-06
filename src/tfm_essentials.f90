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
MODULE tfm_essentials
  !-----------------------------------------------------------------------------
  ! modules used
  !-----------------------------------------------------------------------------
  USE, INTRINSIC :: iso_fortran_env, ONLY: dp => real64

  USE tfm_constants, ONLY : &
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
  PUBLIC ::                    &
    tfm_essentials_mean_acc,   &
    tfm_essentials_do_nothing, &
    tfm_essentials_indicate_tstep

  !-----------------------------------------------------------------------------
  ! routines of this module
  !-----------------------------------------------------------------------------
  CONTAINS


  SUBROUTINE tfm_essentials_do_nothing( &
  &  nz,                                &
  &  variable                           &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Subroutne: tfm_essentials_do_nothing
    !
    ! Routine to do "nothing" with a given variable. Used to avoid warning
    ! at compile time.
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of the given variable.
    !   variable: Variable to do nothing about.
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN)                 :: nz
    REAL(dp), DIMENSION(nz), INTENT(IN) :: variable
    REAL(dp), DIMENSION(nz)             :: nothing

    !---------------------------------------------------------------------------

    nothing = variable
  END SUBROUTINE tfm_essentials_do_nothing


  SUBROUTINE tfm_essentials_mean_acc( &
  &  nz,                              &
  &  depth,                           &
  &  density,                         &
  &  age,                             &
  &  mean_acc                         &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)
    !---------------------------------------------------------------------------
    ! Subroutine : tfm_essentials_mean_acc
    !
    ! The subroutine computes the mean accumulation rate from the age of
    ! the firn profile. The concept follow the idea of calculating the
    ! mean accumulation rate over the life time of a firn parcel.
    !
    ! See for example:
    ! Stevens, C. M., Verjans, V., Luding, J. M. D., Kahle, E. C.,
    ! Horlings, A. N., Horlings, B. I., and Waddington, E. D. The Community
    ! Firn Model (CFM) v1.0. Geosci. Model. Dev., 13, 4355-4377, (2020).
    ! https://doi.org/10.5194/gmd-13-4355-2020
    !
    ! Author: Timm Schultz
    !
    ! Arguments:
    !   nz: Dimension of variables "depth", "density", "age", "mean_acc".
    !   depth: Depth of the firn profile (m).
    !   density: Density of the firn profile (kg m**-3).
    !   age: Age of the firn profile (s).
    !   mean_acc - on input: Variable to store the mean accumulation rate.
    !
    ! Result:
    !   mean_acc - on output: Mean accumulation rate (m weq. a**-1).
    !---------------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nz
    REAL(dp), DIMENSION(nz), INTENT(IN) :: &
      depth,                               &
      density,                             &
      age

    REAL(dp), DIMENSION(nz), INTENT(INOUT) :: mean_acc

    INTEGER :: n

    !---------------------------------------------------------------------------

    mean_acc(nz) = 0.0_dp
    DO n = nz - 1, 1, -1
      mean_acc(n) = (depth(n+1) - depth(n)) * (density(n) / WATER_DENSITY)
      mean_acc(n) = mean_acc(n) + mean_acc(n+1)
    END DO

    mean_acc(1:nz-1) = mean_acc(1:nz-1) / (age(1:nz-1) / SECONDS_YEAR)
  END SUBROUTINE tfm_essentials_mean_acc


  SUBROUTINE tfm_essentials_indicate_tstep( &
  &  nt,                                    &
  &  t                                      &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)

    INTEGER, INTENT(IN) :: &
      nt,                  &
      t

    REAL :: perc

    INTEGER :: &
      n,       &
      n_perc

    ! some computations
    perc = (real(t) / real(nt)) * 100.0
    n_perc = floor(perc / 5.0)

    ! number of time steps
    WRITE(*, '(a)', ADVANCE='no') '\rtimestep: '
    WRITE(*, '(i6,a,i6)', ADVANCE='no') t, ' of ', nt

    ! bar
    WRITE(*, '(a)', ADVANCE='no') ' |'
    DO n = 1, 20, 1
       IF ( n < n_perc ) THEN
         WRITE(*, '(a)', ADVANCE='no') '='
       ELSE IF (n == n_perc) THEN
         WRITE(*, '(a)', ADVANCE='no') '>'
       ELSE
         WRITE(*, '(a)', ADVANCE='no') ' '
       END IF
    END DO
    WRITE(*, '(a)', ADVANCE='no') '|'

    ! percentenge
    WRITE(*, '(f6.2,a)', ADVANCE='no') perc, ' %'

    ! done
    IF ( t == nt ) THEN
      WRITE(*, '(a)') ''
      WRITE(*, '(a)') 'Done!'
    END IF
  END SUBROUTINE tfm_essentials_indicate_tstep
END MODULE tfm_essentials

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
PROGRAM tfm_example
  !-----------------------------------------------------------------------------
  ! modules used
  !-----------------------------------------------------------------------------
  USE, INTRINSIC :: iso_fortran_env, ONLY: &
    dp => real64

  USE tfm_num, ONLY:           &
    tfm_num_modelinit,         &
    sim_models,                &
    tfm_num_trimProfileLength, &
    tfm_num_trimProfileAge,    &
    tfm_num_surface,           &
    tfm_num_step

  USE tfm_essentials, ONLY: &
    tfm_essentials_indicate_tstep

  USE tfm_io, ONLY : &
    simpleOutput,    &
    tfm_file_length, &
    tfm_read_init,   &
    tfm_read_csv

  USE tfm_llStructure, ONLY: &
    llProps,                 &
    llPropsFree

  USE tfm_constants, ONLY: &
    SECONDS_YEAR

  !-----------------------------------------------------------------------------
  ! declarrations
  !-----------------------------------------------------------------------------
  IMPLICIT NONE (TYPE, EXTERNAL)

  ! Parameters
  CHARACTER(LEN=*), PARAMETER :: CONFIGURATION_FILE = './tfm.conf'

  ! input variables read from configuration file
  CHARACTER(len=100) ::                    &
    solve_density,                         &
    solve_temperature,                     &
    solve_heat_capacity,                   &
    solve_thermal_conductivity,            &
    solve_liquid_thermal_conductivity,     &
    solve_saturation_thermal_conductivity, &
    solve_liquid,                          &
    solve_van_genuchten,                   &
    solve_grain_growth

  character(len=100) :: &
    forcing_input_file, &
    initial_input_file

  REAL(KIND=dp) ::      &
    time_step,          &
    spinup,             &
    max_profile_length, &
    max_profile_age

  TYPE(sim_models)                           :: models
  TYPE(llProps)                              :: props
  REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: forcing
  REAL(KIND=dp), DIMENSION(:), ALLOCATABLE   :: runoff
  REAL(KIND=dp), DIMENSION(7)                :: spinup_forcing

  INTEGER :: &
    nt = 0,  &
    snt = 0, &
    tic,     &
    toc,     &
    t
  REAL(KIND=dp) :: rate

  !-----------------------------------------------------------------------------

  ! configuration file list
  NAMELIST /config/                         &
  &  solve_density,                         &
  &  solve_temperature,                     &
  &  solve_heat_capacity,                   &
  &  solve_thermal_conductivity,            &
  &  solve_liquid_thermal_conductivity,     &
  &  solve_saturation_thermal_conductivity, &
  &  solve_liquid,                          &
  &  solve_van_genuchten,                   &
  &  solve_grain_growth,                    &
  &  forcing_input_file,                    &
  &  initial_input_file,                    &
  &  time_step,                             &
  &  spinup,                                &
  &  max_profile_length,                    &
  &  max_profile_age

  ! read configuration
  OPEN(111, FILE=CONFIGURATION_FILE, ACTION='read')
    READ(UNIT=111, NML=config)
  CLOSE(111)

  ! model initialization
  CALL tfm_num_modelinit(                                                               &
  &  solve_density=trim(solve_density),                                                 &
  &  solve_temperature=trim(solve_temperature),                                         &
  &  solve_heat_capacity=trim(solve_heat_capacity),                                     &
  &  solve_thermal_conductivity=trim(solve_thermal_conductivity),                       &
  &  solve_liquid_thermal_conductivity=trim(solve_liquid_thermal_conductivity),         &
  &  solve_saturation_thermal_conductivity=trim(solve_saturation_thermal_conductivity), &
  &  solve_liquid=trim(solve_liquid),                                                   &
  &  solve_van_genuchten=trim(solve_van_genuchten),                                     &
  &  solve_grain_growth=trim(solve_grain_growth),                                       &
  &  models=models                                                                      &
  )

  ! sinup configuration
  IF ( spinup > 0.0 ) THEN
    snt = int((spinup * SECONDS_YEAR) / time_step)
  ELSE
    snt = 0
  END IF

  ! forcing and init import
  CALL tfm_file_length(trim(forcing_input_file), nt)
  ALLOCATE(forcing(7,nt), runoff(nt))
  CALL tfm_read_csv(trim(forcing_input_file), 7, nt, forcing)
  CALL tfm_read_init(trim(initial_input_file), props, models)

  ! freedback
  PRINT '(a,a)', 'Model Definitions'
  PRINT '(a,a)', '================='
  PRINT '(a,a)', 'density:                         ', trim(solve_density)
  PRINT '(a,a)', 'temperature:                     ', trim(solve_temperature)
  PRINT '(a,a)', 'heat capacity:                   ', trim(solve_heat_capacity)
  PRINT '(a,a)', 'thermal conductivity:            ', trim(solve_thermal_conductivity)
  PRINT '(a,a)', 'liquid thermal conductivity:     ', trim(solve_liquid_thermal_conductivity)
  PRINT '(a,a)', 'saturation thermal conductivity: ', trim(solve_saturation_thermal_conductivity)
  PRINT '(a,a)', 'liquid:                          ', trim(solve_liquid)
  PRINT '(a,a)', 'van Genuchten:                   ', trim(solve_van_genuchten)
  PRINT '(a,a)', 'grain growth:                    ', trim(solve_grain_growth)
  PRINT '(a,a)', ''
  PRINT '(a,a)', 'read forcing from:         ', trim(forcing_input_file)
  PRINT '(a,a)', 'read initial profile from: ', trim(initial_input_file)
  PRINT '(a,a)', ''
  PRINT '(a,f13.2)', 'time step (s):            ', time_step
  PRINT '(a,i10)',   'number of time steps:     ', nt
  PRINT '(a,i10)',   'points in inital profile: ', props%depth%length
  PRINT *, ''

  ! optional spinup using mean values for forcing
  CALL system_clock(tic, rate)
  IF ( spinup > 0.0 ) THEN

    PRINT '(a)', 'Spin-Up:'

    ! mean forcing for spinupp
    spinup_forcing = [           &
    &  0.0_dp,                   & ! -> world time
    &  0.0_dp,                   & ! -> model time
    &  (sum(forcing(3,:)) / nt), & ! -> surface temperature
    &  (sum(forcing(4,:)) / nt), & ! -> surface densiy
    &  (sum(forcing(5,:)) / nt), & ! -> surface accumulation (solid)
    &  0.0_dp,                   & ! -> surface accumulation (liquid)
    &  (sum(forcing(7,:)) / nt)  & ! -> surface grain radius
    ]

    DO t = 1, snt, 1

      CALL tfm_essentials_indicate_tstep(snt, t)

      IF (max_profile_length > 0.0) then
        CALL tfm_num_trimProfileLength(props, max_profile_length)
      END IF

      IF (max_profile_age > 0.0) then
        CALL tfm_num_trimProfileAge(props, max_profile_age)
      END IF


      CALL tfm_num_surface( &
      &  time_step,         &
      &  spinup_forcing,    &
      &  models,            &
      &  props              &
      )

      CALL tfm_num_step(              &
      &  time_step,                   &
      &  models=models,               &
      &  props=props,                 &
      &  runoff=runoff(1),            &
      &  liquid_acc=spinup_forcing(6) &
      )
    END DO

    PRINT '(a)', ''
  END IF

  CALL simpleOutput(0, props)

  ! time loop
  PRINT '(a)', 'Simulation run:'
  DO t = 1, nt, 1

    CALL tfm_essentials_indicate_tstep(nt, t)

    IF (max_profile_length > 0.0) THEN
      CALL tfm_num_trimProfileLength(props, max_profile_length)
    END IF

    IF (max_profile_age > 0.0) THEN
      CALL tfm_num_trimProfileAge(props, max_profile_age)
    END IF

    CALL tfm_num_surface( &
    &  time_step,         &
    &  forcing(:,t),      &
    &  models,            &
    &  props              &
    )

    CALL tfm_num_step(         &
    &  time_step,              &
    &  models=models,          &
    &  props=props,            &
    &  runoff=runoff(t),       &
    &  liquid_acc=forcing(6,t) &
    )

    CALL simpleOutput(t, props)
  END DO

  ! feedback
  CALL system_clock(toc, rate)
  PRINt *, ''
  WRITE(*, '(a,f10.2,a)') 'time elapsed: ', real(toc - tic) / real(rate), ' s'

  ! memoray deallocation
  CALL llPropsFree(props)
  DEALLOCATE(forcing, runoff)
END PROGRAM tfm_example

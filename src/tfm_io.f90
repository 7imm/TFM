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
MODULE tfm_io
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
  PUBLIC ::          &
    tfm_file_length, &
    tfm_read_csv,    &
    tfm_read_init,   &
    simpleOutput

  !-----------------------------------------------------------------------------
  ! routines of this module
  !-----------------------------------------------------------------------------
  CONTAINS


  SUBROUTINE tfm_file_length( &
  &  input_file,              &
  &  length                   &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)

    CHARACTER(LEN=*), INTENT(IN) :: input_file
    INTEGER, INTENT(INOUT)       :: length

    INTEGER :: stat

    length = 0

    OPEN(111, FILE=input_file, ACTION='read')
      DO
        READ(111, FMT='(a)', IOSTAT=stat)
        IF ( stat /= 0 ) EXIT
        length = length + 1
      END DO
    CLOSE(111)
  END SUBROUTINE tfm_file_length


  SUBROUTINE tfm_read_csv( &
  &  input_file,           &
  &  width,                &
  &  length,               &
  &  arr                   &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)

    CHARACTER(LEN=*), INTENT(IN) :: input_file
    INTEGER, INTENT(IN) :: &
      width,               &
      length
    REAL(dp), DIMENSION(width,length), INTENT(INOUT) :: arr

    INTEGER :: n

    OPEN(111, FILE=input_file, ACTION='read')
      DO n = 1, length, 1
        READ(111,*) arr(:,n)
      END DO
    CLOSE(111)
  END SUBROUTINE tfm_read_csv


  SUBROUTINE tfm_read_init( &
  &  input_file,            &
  &  props,                 &
  &  models                 &
  )
    USE tfm_num, ONLY : sim_models
    USE tfm_llStructure, ONLY : &
      llProps,                  &
      llAppendData,             &
      llGetDAta
    IMPLICIT NONE (TYPE, EXTERNAL)

    CHARACTER(LEN=*), INTENT(IN) :: input_file
    TYPE(llProps), INTENT(INOUT) :: props
    TYPE(sim_models), INTENT(IN) :: models

    INTEGER :: &
      stat,    &
      nz
    REAL(dp), DIMENSION(6) :: init

    REAL(dp), DIMENSION(:), ALLOCATABLE :: &
      heatcap,                             &
      thermcond

    OPEN(111, FILE=input_file, ACTION='read')
      DO
        READ(111, *, IOSTAT=stat) init
        IF ( stat /= 0 ) EXIT

        CALL llAppendData(props%depth,        1, init(1))
        CALL llAppendData(props%density,      1, init(2))
        CALL llAppendData(props%temperature,  1, init(3))
        CALL llAppendData(props%grain_radius, 1, init(4))
        CALL llAppendData(props%liquidwater,  1, init(5))
        CALL llAppendData(props%age,          1, init(6))
      END DO
    CLOSE(111)

    nz = props%depth%length
    ALLOCATE(heatcap(nz), thermcond(nz))

    heatcap = models%heatcap_model(  &
    &  nz,                           &
    &  llGetData(props%density),     &
    &  llGetData(props%temperature), &
    &  llGetData(props%liquidwater)  &
    )
    thermcond = models%thermcond_model( &
    &  nz,                              &
    &  llGetData(props%density),        &
    &  llGetData(props%temperature)     &
    )

    CALL llAppendData(props%heatcap, nz, heatcap)
    CALL llAppendData(props%thermcond, nz, thermcond)

    DEALLOCATE(heatcap, thermcond)
  END SUBROUTINE tfm_read_init


  SUBROUTINE simpleOutput( &
  &  nt,                   &
  &  props                 &
  )
    USE tfm_llStructure, ONLY: &
      llProps,                 &
      llGetDAta
    IMPLICIT NONE (TYPE, EXTERNAL)

    INTEGER, INTENT(IN)       :: nt
    TYPE(llProps), INTENT(IN) :: props
    CHARACTER(LEN=6)          :: out_number
    CHARACTER(LEN=24)         :: out_name

    INTEGER :: n

    REAL(dp), DIMENSION(props%depth%length) :: &
      depth,                                   &
      density,                                 &
      temperature,                             &
      grain_radius,                            &
      age,                                     &
      liquidwater

    depth        = llGetData(props%depth)
    density      = llGetData(props%density)
    temperature  = llGetData(props%temperature)
    grain_radius = llGetData(props%grain_radius)
    age          = llGetData(props%age)
    liquidwater  = llGetData(props%liquidwater)

    WRITE (out_number, '(I0.6)') nt
    out_name = 'tfm_output/tfm'//out_number//'.out'

    OPEN(333, FILE=out_name, STATUS='replace', ACTION='write')
      DO n = props%depth%length, 1, -1
        WRITE(333,*) depth(n), density(n), temperature(n), grain_radius(n), age(n), liquidwater(n)
      END DO
    CLOSE(333)
  END SUBROUTINE simpleOutput
END MODULE tfm_io

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
MODULE tfm_llStructure
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
    llProps,         &
    llGetLast,       &
    llGetData,       &
    llUpdateList,    &
    llDropData,      &
    llAppendData,    &
    llPropsDropData, &
    llPropsFree

  !-----------------------------------------------------------------------------
  ! types
  !-----------------------------------------------------------------------------
  TYPE llNode
    REAL(dp), DIMENSION(:), ALLOCATABLE :: data
    TYPE(llNode), POINTER               :: next
  END TYPE llNode


  TYPE linkedList
    INTEGER ::          &
      node_size = 1000, &
      hind = 1,         &
      tind = 1001,      &
      length = 0

    TYPE(llNode), POINTER :: &
      head => null(),        &
      tail => null()
  END TYPE linkedList


  TYPE llProps
    TYPE(linkedList) :: &
      depth,            &
      density,          &
      temperature,      &
      heatcap,          &
      thermcond,        &
      liquidwater,      &
      age,              &
      grain_radius
  END TYPE llProps

  !-----------------------------------------------------------------------------
  ! routines of the module
  !-----------------------------------------------------------------------------
  CONTAINS


  SUBROUTINE llSetNodeSize( &
  &  self,                  &
  &  node_size              &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)

    TYPE(linkedList), INTENT(INOUT) :: self
    INTEGER, INTENT(IN)             :: node_size

    self%node_size = node_size
    self%tind = (self%tind + 1)
  END SUBROUTINE llSetNodeSize


  SUBROUTINE llAppendNode(self)
    IMPLICIT NONE (TYPE, EXTERNAL)

    TYPE(linkedList), INTENT(INOUT) :: self
    TYPE(llNode), POINTER           :: new_node

    ! generate new node
    ALLOCATE(new_node)
    ALLOCATE(new_node%data(self%node_size))
    new_node%next => null()

    ! if the appended node is the first node, head and tail are the new node
    IF ( .NOT. associated(self%head) ) THEN
      self%head => new_node
      self%tail => new_node

    ! otherwise the new node bcomes the tail
    ELSE
      self%tail%next => new_node
      self%tail      => new_node

    END IF

    ! change linked list attributes accordingly
    self%tind = 1
  END SUBROUTINE llAppendNode


  SUBROUTINE llAppendData( &
  &  self,                 &
  &  n,                    &
  &  data                  &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)

    TYPE(linkedList), INTENT(INOUT)    :: self
    INTEGER, INTENT(IN)                :: n
    REAL(dp), DIMENSION(n), INTENT(IN) :: data

    INTEGER ::   &
      remainder, &
      free,      &
      k,         &
      l,         &
      i,         &
      j

    remainder = n
    i = 1

    DO WHILE ( remainder > 0 )

      ! if necessary append another node
      IF ( self%tind > self%node_size ) THEN
        CALL llAppendNode(self)
      END IF

      ! node and data indices
      free = (self%node_size - self%tind + 1)
      k = self%tind
      l = min((k + free - 1), (k + remainder - 1))
      j = (i + (l - k))

      ! store data
      self%tail%data(k:l) = data(i:j)

      ! change counting variables accordingly
      self%tind = (self%tind + (l - k) + 1)
      i = (i + (l - k) + 1)
      remainder = (remainder - free)
    END DO

    self%length = (self%length + n)
  END SUBROUTINE llAppendData


  SUBROUTINE llDropHeadNode(self)
    IMPLICIT NONE(TYPE, EXTERNAL)

    TYPE(linkedList), INTENT(INOUT) :: self
    TYPE(llNode), POINTER           :: former_head_node

    former_head_node => self%head
    self%head => self%head%next

    DEALLOCATE(former_head_node%data)
    NULLIFY(former_head_node)
  END SUBROUTINE llDropHeadNode


  SUBROUTINE llDropTailNode(self)
    IMPLICIT NONE (TYPE, EXTERNAL)

    TYPE(linkedList), INTENT(INOUT) :: self
    TYPE(llNode), POINTER           :: new_tail_node

    new_tail_node => self%head
    DO WHILE ( associated(new_tail_node%next%next) )
      new_tail_node => new_tail_node%next
    END DO

    DEALLOCATE(self%tail%data)
    self%tail => new_tail_node
    NULLIFY(self%tail%next)
  END SUBROUTINE llDropTailNode


  SUBROUTINE llDropData( &
  &  self,               &
  &  n                   &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)

    TYPE(linkedList), INTENT(INOUT) :: self
    INTEGER, INTENT(IN)             :: n

    IF ( n > 0 ) THEN
      self%hind = (self%hind + n)

      ! if the head iis empty delete the current head node
      DO WHILE ( self%hind >= self%node_size )
        CALL llDropHeadNode(self)
        self%hind = (self%hind - self%node_size)
      END DO

    ELSE IF ( n < 0 ) THEN
      self%tind = (self%tind + n)

      ! if the tail is empty delete the current tail node
      DO WHILE ( self%tind < 1 )
        CALL llDropTailNode(self)
        self%tind = (self%tind + self%node_size)
      END DO

    ELSE
      CONTINUE
    END IF

    self%length = (self%length - abs(n))
  END SUBROUTINE llDropData


  FUNCTION llGetData(self) RESULT(output)
    IMPLICIT NONE (TYPE, EXTERNAL)

    TYPE(linkedList), INTENT(IN)     :: self
    REAL(dp), DIMENSION(self%length) :: output

    TYPE(llNode), POINTER :: curr_node

    INTEGER :: &
      k,       &
      l,       &
      i,       &
      j

    ! data from the head
    i = 1
    j = min((self%node_size - self%hind + 1), self%length)
    k = self%hind
    l = min(self%node_size, (self%hind + self%length - 1))
    output(i:j) = self%head%data(k:l)

    IF ( associated(self%head, self%tail) ) RETURN

    ! data from the middle part
    curr_node => self%head%next
    DO WHILE ( associated(curr_node%next) )
      i = (j + 1)
      j = (i + self%node_size - 1)
      output(i:j) = curr_node%data(1:l)

      curr_node => curr_node%next
    END DO

    ! data from the tail
    i = (j + 1)
    j = (i + (self%tind - 1) - 1)
    l = (self%tind - 1)
    output(i:j) = self%tail%data(1:l)
  END FUNCTION llGetData


  SUBROUTINE llUpdateList( &
  &  self,                 &
  &  data                  &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)

    TYPE(linkedList), INTENT(INOUT)              :: self
    REAL(dp), DIMENSION(self%length), INTENT(IN) :: data

    INTEGER :: &
      k,       &
      l,       &
      i,       &
      j

    TYPE(llNode), POINTER :: curr_node

    ! head
    i = 1
    j = min((self%node_size - self%hind + 1), self%length)
    k = self%hind
    l = min(self%node_size, (self%hind + self%length - 1))
    self%head%data(k:l) = data(i:j)

    IF ( associated(self%head, self%tail) ) RETURN

    ! middle part
    curr_node => self%head%next
    DO WHILE ( associated(curr_node%next) )
      i = (j + 1)
      j = (i + self%node_size - 1)
      curr_node%data(:) = data(i:j)

      curr_node => curr_node%next
    END DO

    ! tail
    i = (j + 1)
    j = (i + (self%tind - 1) - 1)
    k = 1
    l = (self%tind - 1)
    self%tail%data(k:l) = data(i:j)
  END SUBROUTINE llUpdateList


  FUNCTION llGetFirst(self) RESULT(first)
    IMPLICIT NONE (TYPE, EXTERNAL)

    TYPE(linkedList), INTENT(IN) :: self
    REAL(dp)                     :: first

    first = self%head%data(self%hind)
  END FUNCTION llGetFirst


  FUNCTION llGetLast(self) RESULT(last)
    IMPLICIT NONE (TYPE, EXTERNAL)

    TYPE(linkedList), INTENT(IN)     :: self
    REAL(dp)                         :: last
    REAL(dp), DIMENSION(self%length) :: all_data

    IF ( (self%tind - 1) /= 0 ) THEN
      last = self%tail%data((self%tind - 1))

    ! if the last index is the first of the tail
    ! get the last entrace of all data
    ELSE
      all_data = llGetData(self)
      last = all_data(self%length)

    END IF
  END FUNCTION llGetLast


  SUBROUTINE llFreeList(self)
    IMPLICIT NONE (TYPE, EXTERNAL)

    TYPE(linkedList), INTENT(INOUT) :: self

    TYPE(llNode), POINTER :: &
      curr_node,             &
      next_node

    curr_node => self%head
    DO WHILE ( associated(curr_node) )
      next_node => curr_node%next

      DEALLOCATE(curr_node%data)
      NULLIFY(curr_node)

      curr_node => next_node
    END DO
  END SUBROUTINE llFreeList


  SUBROUTINE llPropsFree(self)
    IMPLICIT NONE (TYPE, EXTERNAL)

    TYPE(llProps), INTENT(INOUT) :: self

    CALL llFreeList(self%depth)
    CALL llFreeList(self%density)
    CALL llFreeList(self%grain_radius)
    CALL llFreeList(self%temperature)
    CALL llFreeList(self%heatcap)
    CALL llFreeList(self%thermcond)
    CALL llFreeList(self%liquidwater)
    CALL llFreeList(self%age)
  END SUBROUTINE llPropsFree


  SUBROUTINE llPropsDropData( &
  &  self,                    &
  &  n                        &
  )
    IMPLICIT NONE (TYPE, EXTERNAL)

    TYPE(llProps), INTENT(INOUT) :: self
    INTEGER, INTENT(IN)          :: n

    CALL llDropData(self%depth,        n)
    CALL llDropData(self%density,      n)
    CALL llDropData(self%grain_radius, n)
    CALL llDropData(self%temperature,  n)
    CALL llDropData(self%heatcap,      n)
    CALL llDropData(self%thermcond,    n)
    CALL llDropData(self%liquidwater,  n)
    CALL llDropData(self%age,          n)
  END SUBROUTINE llPropsDropData
END MODULE tfm_llStructure

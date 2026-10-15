!!****ih* source/Grid/GridMain/UG/Grid_iterator
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distributed under the License is distributed on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
!!
!!
!!
!!****

#include "Simulation.h"
#include "constants.h"
!! defines IMPURE_ELEMENTAL:
#include "FortranLangFeatures.fh"

#ifdef DEBUG_ALL
#define DEBUG_GRID_ITER
#endif

module Grid_iterator

    implicit none

    private

    public :: build_iterator, destroy_iterator

    !!****ic* Grid_iterator/Grid_iterator_t
    !!
    !! NAME
    !!  Grid_iterator_t
    !!
    !!****
    type, public :: Grid_iterator_t
        integer :: curBlk = 2
        integer :: lev = INVALID_LEVEL
        integer, private :: tileSize(MDIM) = (/1,1,1/)
        integer, private :: nxt=1,nyt=1,nzt=1
        logical, private :: tiling         = .FALSE.
        integer, private :: curTile        = 0 ! local tileNo in the block, 0-based
        integer, private :: ntiles         = 0 ! number of local tiles in the block
#ifdef DEBUG_GRID_ITER
        integer, private :: index_         = 0 ! 1-based index of the 'current item'
#endif
        integer, private :: nthreads       = 1 ! number of threads to participate in iterations
#ifdef _OPENMP
        integer, private :: threadOffs     = 0
#endif
    contains
        procedure, public :: isValid
        procedure         :: next_ ! simple next
        procedure         :: nextN ! next by n
        procedure, public :: next
        procedure, public :: currentTile
#ifdef DEBUG_GRID_ITER
        procedure, public :: itemNumber
#endif
    end type Grid_iterator_t

contains

    !!****im* Grid_iterator_t/build_iterator
    !!
    !! NAME
    !!  build_iterator
    !!
    !! SYNOPOSIS
    !!  build_iterator(Grid_iterator_t(OUT)  :: itor,
    !!                 integer(IN)           :: nodetype,
    !!                 integer(IN), optional :: level,
    !!                 logical(IN), optional :: tiling,
    !!                 integer(IN), optional :: nthreads)
    !!
    !! DESCRIPTION
    !!  Construct an iterator for walking across a specific subset of blocks or
    !!  tiles within the current paramesh octree structure.  The iterator is already
    !!  set to the first matching block/tile.
    !!  In the UG implementation, there is exactly one block on each MPI rank.
    !!  
    !! ARGUMENTS
    !!  itor     - the iterator
    !!  nodetype - the class of blocks to iterate over (e.g. LEAF, ACTIVE_BLKS).
    !!             in the current UG implementation, passing an unsupported
    !!             value results in program abort.
    !!  level    - iterate only over leaf blocks/tiles located at this level of
    !!             refinement.
    !!             A level value of UNSPEC_LEVEL is equivalent to omitting
    !!             this optional argument.
    !!             In this UG implementation, there is only one block at level 1.
    !!             Therefore this argument may be ignored in the current UG implementation.
    !!  tiling   - an optional optimization hint.  If TRUE, proper tiling is
    !!             requested.
    !!  nthreads - an optional argument that may affect division of the iteration
    !!             space in an active OpenMP parallel region.
    !!             Can usually be omitted since the default behavior is reasonable.
    !!
    !! SEE ALSO
    !!  constants.h
    !!****
    subroutine build_iterator(itor, nodetype, level, tiling, tileSize, nthreads)
      use Grid_data, ONLY: gr_enableTiling, gr_useTiling, gr_tileSize
      use Grid_data, ONLY:           gr_ilo, gr_ihi, &
                                     gr_jlo, gr_jhi, &
                                     gr_klo, gr_khi
      !$ use omp_lib
        type(Grid_iterator_t), intent(OUT)          :: itor
        integer,               intent(IN)           :: nodetype
        integer,               intent(IN), optional :: level
        logical,               intent(IN), optional :: tiling
        integer,               intent(IN), optional :: tileSize(1:MDIM)
        integer,               intent(IN), optional :: nthreads

        if (present(nthreads)) then
           itor%nthreads = nthreads
        else
#ifdef _OPENMP
           itor%nthreads = omp_get_num_threads()
#else
           itor%nthreads = 1
#endif
        end if
#ifdef _OPENMP
        itor%threadOffs = omp_get_thread_num()
        if (itor%threadOffs .GE. itor%nthreads) itor%threadOffs = 0
#endif

        itor%curBlk = 1
        itor%lev = 1

        if (present(level)) then
            itor%lev = level
            if      (itor%lev == UNSPEC_LEVEL) then
                itor%lev = 1    
            else if (itor%lev /= 1) then
                itor%curBlk = 2
            end if
        end if

        if (gr_enableTiling) then ! else keep default value of .FALSE.
           if (present(tiling)) then
              itor%tiling = tiling
           else
              itor%tiling = gr_useTiling
           end if
        end if
        if (present(tileSize)) then
           itor%tileSize(1:NDIM) = tileSize(1:NDIM)
        else
           itor%tileSize(:) = gr_tileSize(:)
        end if

        if ((nodetype /= LEAF) .AND. (nodetype /= ALL_BLKS) .AND. (nodetype /= ACTIVE_BLKS)) then
            call Driver_abort("[Grid_iterator]: Unsupported nodetype")
        end if
        if (itor%tiling .AND. itor%curBlk == 1) then
           itor%nxt = max(1,(gr_ihi-gr_ilo+itor%tileSize(IAXIS))/itor%tileSize(IAXIS))
           itor%ntiles = itor%nxt
#if NDIM > 1
           itor%nyt = max(1,(gr_jhi-gr_jlo+itor%tileSize(JAXIS))/itor%tileSize(JAXIS))
           itor%ntiles = itor%ntiles * itor%nyt
#endif
#if NDIM >= 3
           itor%nzt = max(1,(gr_khi-gr_klo+itor%tileSize(KAXIS))/itor%tileSize(KAXIS))
           itor%ntiles = itor%ntiles * itor%nzt
#endif
        else
          itor% ntiles = 1
        end if
        itor%curTile = 0
#ifdef DEBUG_GRID_ITER
        itor%index_ = 1
#endif
#ifdef _OPENMP
        if (itor%nthreads > 1 .AND. itor%threadOffs > 0) then
           itor%curTile = itor%curTile + itor%threadOffs
           if (itor%curTile .GE. itor%ntiles) itor%curBlk = 2
        end if
#endif
    end subroutine build_iterator

    !!****im* Grid_iterator_t/destroy_iterator
    !!
    !! NAME
    !!  destroy_iterator
    !!
    !! SYNPOSIS
    !!  Destroy given iterator.
    !!
    !! DESCRIPTION
    !!  Clean-up block interator object at destruction
    !!
    !!****
    IMPURE_ELEMENTAL subroutine destroy_iterator(itor)
        type(Grid_iterator_t), intent(INOUT) :: itor

        itor%curBlk = 2
        itor%lev = INVALID_LEVEL
    end subroutine destroy_iterator

    !!****m* Grid_iterator_t/isValid
    !!
    !! NAME
    !!  isValid
    !!
    !! SYNPOSIS
    !!  logical valid = itor%isValid()
    !!
    !! DESCRIPTION
    !!  Determine if the iterator is currently set to a valid block.
    !!
    !! RETURN VALUE 
    !!  True if iterator is currently set to a valid block
    !!
    !!****
    logical function isValid(this)
        class(Grid_iterator_t), intent(IN) :: this

        isValid = (this%curBlk == 1)
    end function isValid

    !!****m* Grid_iterator_t/next_
    !!
    !! NAME
    !!  next_
    !!
    !! SYNPOSIS
    !!  call itor%next_()
    !!
    !! DESCRIPTION
    !!  Advance the iterator to the next tile managed by the process that meets
    !!  the iterator constraints given at instantiation.
    !!
    !!****
    subroutine next_(this)
        class(Grid_iterator_t), intent(INOUT) :: this

        if (this%curBlk == 1 .AND. this%curTile < this%ntiles-1) then
           this%curTile = this%curTile + 1
        else
           this%curBlk = 2
        end if
    end subroutine next_

    !!****m* Grid_iterator_t/nextN
    !!
    !! NAME
    !!  nextN
    !!
    !! SYNPOSIS
    !!  call itor%nextN(n)
    !!
    !! DESCRIPTION
    !!  Advance the iterator up to n times, as if next() was invoked
    !!  until either n invocations were done or isValid becomes false,
    !!  whichever occurs first.
    !!
    !!****
    subroutine nextN(this, n)
        class(Grid_iterator_t), intent(INOUT) :: this
        integer,                intent(IN)    :: n

        integer :: i

        do i = 1,n
           if (this%curBlk .NE. 1) EXIT
           call this%next_()
        end do

    end subroutine nextN

    !!****m* Grid_iterator_t/next
    !!
    !! NAME
    !!  next
    !!
    !! SYNPOSIS
    !!  call itor%next()
    !!
    !! DESCRIPTION
    !!  Advance the iterator to the next tile managed by the current process and thread
    !! that meets the iterator constraints given at instantiation.
    !!
    !!****
    subroutine next(this)
        class(Grid_iterator_t), intent(INOUT) :: this

#ifdef _OPENMP
        call this%nextN(this%nthreads)
#else
        call this%next_()
#endif

#ifdef DEBUG_GRID_ITER
        this%index_ = this%index_ + 1
#endif
    end subroutine next

#ifdef DEBUG_GRID_ITER
    !! NAME
    !!  itemNumber
    !!
    !! SYNPOSIS
    !!  iterno = itor%itemNumber()
    !!
    !! DESCRIPTION
    !!  Return the (1-based) number of the current item, i.e.,
    !!  the item to which an following call of the currentTile
    !!  method would refer, if the current state of the iterator
    !!  is valid (points to a valid item).
    !!  If the iterator has proceeded beyond the last item,
    !!  returns 1 plus the number of iterations performed.
    elemental integer function itemNumber(this)
        class(Grid_iterator_t), intent(IN) :: this

        itemNumber = this%index_
    end function itemNumber
#endif

    !!****m* Grid_iterator_t/currentTile
    !!
    !! NAME
    !!  currentTile 
    !!
    !! SYNPOSIS
    !!  call itor%currentTile(Grid_tile_t(OUT) : tileDesc)
    !!
    !! DESCRIPTION
    !!  Obtain meta data that characterizes the block currently set in the
    !!  iterator.
    !!
    !!****
    subroutine currentTile(this, tileDesc)
        use Grid_tile,        ONLY : Grid_tile_t
        use Grid_data,        ONLY : gr_blkCornerID, &
                                     gr_ilo, gr_ihi, &
                                     gr_jlo, gr_jhi, &
                                     gr_klo, gr_khi, &
                                     gr_iloGc, gr_ihiGc, &
                                     gr_jloGc, gr_jhiGc, &
                                     gr_kloGc, gr_khiGc
        use Driver_interface, ONLY : Driver_abort

        class(Grid_iterator_t), intent(IN)  :: this
        type(Grid_tile_t),      intent(OUT) :: tileDesc

        integer :: loBnd(1:MDIM), hiBnd(1:MDIM)
        integer :: tx,ty,tz
        if (.NOT. this%isValid()) then
            call Driver_abort("[currentTile] No current tile")
        end if

        tileDesc%id = this%curBlk
        tileDesc%level = this%lev
        tileDesc%cid = gr_blkCornerID
        tileDesc%stride = 1

        associate(lo      => tileDesc%limits(LOW, :), &
                  hi      => tileDesc%limits(HIGH, :), &
                  loGrown => tileDesc%grownLimits(LOW, :), &
                  hiGrown => tileDesc%grownLimits(HIGH, :), &
                  loBlkGC => tileDesc%blkLimitsGC(LOW, :), &
                  hiBlkGC => tileDesc%blkLimitsGC(HIGH, :))
            loBlkGC(1:MDIM) = [gr_iloGc, gr_jloGc, gr_kloGc]
            hiBlkGC(1:MDIM) = [gr_ihiGc, gr_jhiGc, gr_khiGc]

            lo(1:MDIM)   = [gr_ilo,   gr_jlo,   gr_klo]
            hi(1:MDIM)   = [gr_ihi,   gr_jhi,   gr_khi]

            ! If there is no tiling with UG, the grown tile
            ! is just the block + GC halo
            loGrown(:) = loBlkGC(:)
            hiGrown(:) = hiBlkGC(:)

            if (this%tiling) then
               call tile2txtytz(this%curTile,tx,ty,tz)
               loBnd(:) = max(0, 1 - (/tx,ty,tz/) ) ! indicates whether range includes first cell
!!$            hiBnd(:) = 1 - min(1, (/this%nxt-1,this%nyt-1,this%nzt-1/) - (/tx,ty,tz/) )
               hiBnd(:) = max(0, (/tx,ty,tz/) - (/this%nxt-2,this%nyt-2,this%nzt-2/) )
               loBnd(:) = loBnd(:) * (lo(:) - loBlkGC(:))
               hiBnd(:) = hiBnd(:) * (hiBlkGC(:) - hi(:))

               hi(IAXIS) = min(hi(IAXIS), lo(IAXIS)+(gr_ihi-gr_ilo+1)*(tx+1)/this%nxt-1)
               lo(IAXIS) =     lo(IAXIS) + (gr_ihi-gr_ilo+1)*tx/this%nxt
#if (NDIM > 1)
               hi(JAXIS) = min(hi(JAXIS), lo(JAXIS)+(gr_jhi-gr_jlo+1)*(ty+1)/this%nyt-1)
               lo(JAXIS) =     lo(JAXIS) + (gr_jhi-gr_jlo+1)*ty/this%nyt
#endif
#if (NDIM > 2)
               hi(KAXIS) = min(hi(KAXIS), lo(KAXIS)+(gr_khi-gr_klo+1)*(tz+1)/this%nzt-1)
               lo(KAXIS) =     lo(KAXIS) + (gr_khi-gr_klo+1)*tz/this%nzt
#endif
               loGrown(1:NDIM) = max(loBlkGC(1:NDIM), lo(1:NDIM) - loBnd(1:NDIM))
               hiGrown(1:NDIM) = min(hiBlkGC(1:NDIM), hi(1:NDIM) + hiBnd(1:NDIM))
            end if
        end associate

      contains
        subroutine tile2txtytz(tileID,tx,ty,tz)
          integer,intent(IN) :: tileID
          integer,intent(OUT) :: tx,ty,tz

          integer :: tn
          tn = tileID

          associate(nxt      => this%nxt, &
                    nyt      => this%nyt, &
                    nzt      => this%nzt)
            tx = mod(tn,nxt)
            ty = mod(tn/nxt,nyt)
            tz = tn/(nxt*nyt)
          end associate
        end subroutine tile2txtytz

    end subroutine currentTile

end module Grid_iterator


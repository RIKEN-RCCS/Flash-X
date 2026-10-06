! Provider-facing contract. No Grid implementation or mutable global state.
module hydro_grid_contract
  use newimpl_hydro, only: cell_field,flux_field
  implicit none
  integer, parameter :: hydro_kind=selected_real_kind(15,307)
  type :: hydro_block_t
    integer :: block=1,level=1,rank=0
    integer :: lo(3)=1,hi(3)=1,guardLo(3)=1,guardHi(3)=1
    real(hydro_kind) :: spacing(3)=1._hydro_kind,faceArea(3)=1._hydro_kind,volume=1._hydro_kind
  end type
  type, abstract :: hydro_grid_context_t
  end type
  abstract interface
    subroutine grid_guard_interface(context,state,blocks,nVars,status)
      import :: hydro_grid_context_t,cell_field,hydro_block_t
      class(hydro_grid_context_t), intent(inout) :: context
      type(cell_field), intent(inout) :: state
      type(hydro_block_t), intent(in) :: blocks(:)
      integer, intent(in) :: nVars
      integer, intent(out) :: status
    end subroutine
    subroutine grid_flux_interface(context,flux,blocks,fluxMap,nDim,status)
      import :: hydro_grid_context_t,flux_field,hydro_block_t
      class(hydro_grid_context_t), intent(inout) :: context
      type(flux_field), intent(inout) :: flux
      type(hydro_block_t), intent(in) :: blocks(:)
      integer, intent(in) :: fluxMap(5),nDim
      integer, intent(out) :: status
    end subroutine
    subroutine grid_reduce_interface(context,localDt,localSpeed,localLoc,globalDt,globalSpeed,globalLoc,status)
      import :: hydro_grid_context_t,hydro_kind
      class(hydro_grid_context_t), intent(inout) :: context
      real(hydro_kind), intent(in) :: localDt,localSpeed
      integer, intent(in) :: localLoc(5)
      real(hydro_kind), intent(out) :: globalDt,globalSpeed
      integer, intent(out) :: globalLoc(5),status
    end subroutine
    subroutine hydro_eos_interface(context,state,blocks,inputMap,status)
      import :: hydro_grid_context_t,cell_field,hydro_block_t
      class(hydro_grid_context_t), intent(inout) :: context
      type(cell_field), intent(inout) :: state
      type(hydro_block_t), intent(in) :: blocks(:)
      integer, intent(in) :: inputMap(7)
      integer, intent(out) :: status
    end subroutine
  end interface
end module

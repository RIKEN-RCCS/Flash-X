! Single-stage hydrodynamics. Mesh references and storage are backend macros.
module hy_gpu_amr_types
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer, parameter :: rk=kind(1.0)
  integer, parameter, public :: QRHO=1, QUX=2, QUY=3, QUZ=4, QP=5, QGAMMA=6, QEINT=7
  integer, parameter, public :: URHO=1, UMX=2, UMY=3, UMZ=4, UENER=5
  integer, parameter, public :: HY_OK=0, HY_BAD_INPUT=1, HY_BAD_STATE=2, HY_DEGENERATE=3
  type :: cell_field
  real(rk), allocatable :: values(:,:,:,:,:)
end type
type :: flux_field
  real(rk), allocatable :: values(:,:,:,:,:,:)
end type
  type, public :: hydro_parameters
    logical :: flattening=.true., hybrid=.false., shock_marker_present=.false.
    real(rk) :: density_floor=1.e-12_rk, pressure_floor=1.e-12_rk
    real(rk) :: tiny_state=1.e-30_rk, viscosity=0._rk
    real(rk) :: weno_epsilon=1.e-36_rk, flattening_eta=real(max(1.e-99_8,real(tiny(1._rk),8)),rk)
  end type
  type, public :: hydro_workspace
    type(cell_field) :: flat, minus, plus
    integer :: lo(3)=0, hi(3)=-1, blocks=0
  end type
  public :: cell_field, flux_field
end module

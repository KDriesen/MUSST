!< author: Noël Brunetière<br/>&emsp;Arthur Francisco
!  version: 1.0.0
!  date: july, 12 2018
!
!  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!     **MUSST general parameters**
!  </span>
!@note
!        The type *NUM_PAR* is not used in **MSOLV**, it is a type related to the iterative process of **MUSST**
!@endnote
module num_param
use iso_fortran_env, only : input_unit, output_unit
use data_arch, only : I4, R8
implicit none

public

! codes for message
integer(kind=I4), parameter :: NO_MESS    = 0    !! *code for no message on screen during problem solving*
integer(kind=I4), parameter :: PRINT_MESS = 1    !! *code for printing message during problem solving*

integer(kind=I4)   :: VERBOSE             !! *Output detail control*
integer(kind=I4)   :: OPU = output_unit   !! *Output unit*
character(len=128) :: OUTPUT_FILE         !! *When needed, output file*

! num_par definition
type NUM_PAR
   real(kind=R8)    :: relax            !! *relaxation parameter*
   real(kind=R8)    :: eps              !! *error for convergence*
   integer(kind=I4) :: it_max           !! *maximal number of iterations*
   integer(kind=I4) :: mess             !! *message*
endtype NUM_PAR

private :: input_unit, output_unit

endmodule num_param

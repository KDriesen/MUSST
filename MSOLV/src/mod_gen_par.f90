!< author: Noël Brunetière<br/>&emsp;Arthur Francisco
!  version: 1.0.0
!  date: july, 12 2018
!
!  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!     **MUSST general parameters**
!  </span>
module gen_param
use data_arch, only : I4, R8
implicit none

! codes for message
integer(kind=I4), parameter :: NO_MESS    = 0    !! *code for no message on screen during problem solving*
integer(kind=I4), parameter :: PRINT_MESS = 1    !! *code for printing message during problem solving*

integer(kind=I4)   :: SOLV_MESS = NO_MESS        !! *Solver output detail control*

integer(kind=I4)   :: VERBOSE             !! *Output detail control*
character(len=128) :: OUTPUT_FILE         !! *When needed, output file*

endmodule gen_param

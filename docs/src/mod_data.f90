
!! author: Noël Brunetière<br/>&emsp;Arthur Francisco
!! version: 1.0.0
!! date: April,20 2017
!! summary: Data for a fluid film by FEM

!< <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<       **Definition of the data for a film solution problem**
!< </span>

module data_film_hd
use data_arch, only : I4, R8
use fluid_law
implicit none

private

type DATA_FILM
   real(kind=R8)    :: h_0       !! *\(h_0\) : ambient pressure*
   real(kind=R8)    :: h_g       !! *\(h_g\) : cavitation pressure*
   real(kind=R8)    :: V_x       !! *\(V_x\) : surface velocity along \(x\)*
   real(kind=R8)    :: V_y       !! *\(V_y\) : surface velocity along \(y\)*
   type(FLUID)      :: fl        !! *fluid rheological properties*
   integer(kind=I4) :: pb_type   !! *problem type, e.g. HD (hydrodynamic)*
endtype DATA_FILM


integer(kind=I4), parameter :: HD  = 0 !! *code for problem type, 0 : hydrodynamic problem*

public :: DATA_FILM, HD

endmodule data_film_hd

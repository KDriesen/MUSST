
!! author: Noël Brunetière<br/>&emsp;Arthur Francisco
!! version: 1.0.0
!! date: May, 3 2017

!< <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<       **Rheological and thermodynamic behavior of fluid**
!< </span>

!<# Description of the fluid module
!< The fluid module module_fluid_law proposes different rheological and thermodynamic laws for the fluids.
!
   !<##  Fluid type
      !< In the current version, three fluid types are proposed:
      !<
      !<- incompressible liquid
      !<- perfect gas
      !<- mixture of liquid and gaz
      !<
      !< The fluid parameters are:
      !<
      !<- \(\rho_0\),  liquid density;
      !<- \(r_g\),     the perfect gas constant of the considered gas;
      !<- \(\lambda\), the gas mass fraction.
      !<
      !
   !<## Density and compressibility
      !< The module is limited to the density \(\rho\) [[FLUID:rho]] and compressibility
      !< \(\frac{\partial \rho}{\partial p}\) [[FLUID:drhodp]] computation based on the pressure \(p\) and absolute temperature \(t\)
      !
      !<### Incompressible [[fluid_law:INC]]
         !<- \(\rho=\rho_0\)
         !<- \(\frac{\partial \rho}{\partial p}=0\)
      !<### Perfect gas [[fluid_law:GP]]
         !<- \(\rho=\frac{p}{r_g T}\)
         !<- \(\frac{\partial \rho}{\partial p}=\frac{1}{r_g T}\)
      !<### Mixture [[fluid_law:MIXT]]
         !<- \(\rho=\frac{1}{\frac{1-\lambda}{\rho_0}+\frac{\lambda r_g T}{p}}\)
         !<- \(\frac{\partial \rho}{\partial p}=\frac{\lambda \rho_O^2 r_g T}{\left(p(1-\lambda) +\lambda r_g T \rho_0\right)^2}\)

module fluid_law
use data_arch, only : I4, R8
implicit none

private

type FLUID
   integer(kind=I4)             :: fluid_type   !! *constant for the fluid type: INC, GP, MIXT*
   real(kind=R8)                :: p_0          !! *reference pressure*
   real(kind=R8)                :: rho_0        !! *reference liquid density*
   real(kind=R8)                :: mu_0         !! *reference dynamic viscosity*
   real(kind=R8)                :: T_0          !! *reference temperature*
   real(kind=R8)                :: rg           !! *perfect gas constant*
   real(kind=R8)                :: lambda       !! *gas mass fraction*
   real(kind=R8), dimension(20) :: cst          !! *table of parameters for the thermodynamic and rheological laws*
   contains
      procedure :: rho          !! *density*
      procedure :: drhodp       !! *compressibility*
      procedure :: pres         !! *pressure (perfect gas)*
endtype FLUID

! Codes for basic fluid type
integer(kind=I4), parameter :: INC  = 0 !! *incompressible fluid*
integer(kind=I4), parameter :: GP   = 1 !! *perfect gas*
integer(kind=I4), parameter :: MIXT = 2 !! *mixture of liquid and gas*

public :: FLUID, INC, GP, MIXT

contains

   real(kind=R8) function rho(fl, p, t)
   !=========================================================================================
   !< @note Function to calculate the density of the fluid based on *t* and *p*
   !
   !+ if fl is set to [[fluid_law:INC]]  (incompressible) \(\rho=\rho_0\)          where \(\rho_0\) is a member of fl [[FLUID]]
   !+ if fl is set to [[fluid_law:GP]]   (perfect gas)    \(\rho=\frac{p}{r_g T}\) where \(r_g\)    is a member of fl [[FLUID]]
   !+ if fl is set to [[fluid_law:MIXT]] (mixture)        \(\rho=\frac{1}{\frac{1-\lambda}{\rho_0}+\frac{\lambda r_g T}{p}}\) where \(\lambda \text{ and } r_g\) are members of fl [[FLUID]]
   !
   !  @endnote
   !------------------------------------------------------------------------------------------
   implicit none
   class(FLUID),  intent(in) :: fl        !! *fluid type*
   real(kind=R8), intent(in) :: p         !! *pressure*
   real(kind=R8), intent(in) :: t         !! *absolute temperature*
      select case (fl%fluid_type)
         case(INC)    ; rho = fl%rho_0
         case(GP)     ; rho = p / (fl%rg * t)
         case(MIXT)   ; rho = 1._R8 / ((1._R8 - fl%lambda) / fl%rho_0 + (fl%lambda * fl%rg * t) / p)
         case default ; stop 'Bad case number rho'
      endselect
   return
   endfunction rho

   real(kind=R8) function drhodp(fl, p, t)
   !=========================================================================================
   !< @note Function to calculate compressibility (drdp)T based on T and p
   !
   !+ if fl is set to [[fluid_law:INC]]  (incompressible) \(\rho=0\)
   !+ if fl is set to [[fluid_law:GP]]   (perfect gas)    \(\rho=\frac{1}{r_g T}\) where \(r_g\)
   !                                                            is a member of fl [[FLUID]]
   !+ if fl is set to [[fluid_law:MIXT]] (mixture)        \(\rho=\frac{\lambda \rho_0^2 r_g T}{(p(1-\lambda) +
   !                  \lambda r_g T \rho_0)^2}\) where \(\lambda \' , \' \rho_0 \text{ and } r_g\) are members of fl [[FLUID]]
   !  @endnote
   !-----------------------------------------------------------------------------------------
   implicit none
   class(FLUID), intent(in)  :: fl           !! *fluid type*
   real(kind=R8), intent(in) :: p            !! *pressure*
   real(kind=R8), intent(in) :: t            !! *absolute temperature*
      select case (fl%fluid_type)
         case(INC)    ; drhodp = 0._R8
         case(GP)     ; drhodp = 1._R8 / (fl%rg * t)
         case(MIXT)   ; drhodp = fl%lambda * (fl%rho_0**2) * fl%rg * t / (((1._R8 - fl%lambda) * p + fl%lambda * fl%rho_0 * fl%rg * t) ** 2)
         case default ; stop 'Bad case number drhodp'
      endselect
   return
   endfunction drhodp

   real(kind=R8) function pres(fl, rho, t)
   !=========================================================================================
   !< @note Function to calculate pressure (pres) based on T and density
   !
   !+ if fl is set to [[fluid_law:GP]]   (perfect gas)    \(pres=\rho r_g T\) where \(r_g\)
   !                                                            is a member of fl [[FLUID]]
   !+ if fl is set to [[fluid_law:MIXT]] (mixture)        \(pres=\frac{\lambda r_g T}{\frac{1}{\rho} -
   !                                           \frac{1-\lambda}{\rho_0}}\) are members of fl [[FLUID]]
   !  @endnote
   !-----------------------------------------------------------------------------------------
   implicit none
   class(FLUID), intent(in)  :: fl           !! *fluid type*
   real(kind=R8), intent(in) :: rho          !! *density*
   real(kind=R8), intent(in) :: t            !! *absolute temperature*
      select case (fl%fluid_type)
         case(GP)     ; pres = rho * fl%rg * T
         case(MIXT)   ; pres = fl%lambda * fl%rg * T / (1.0_R8 / rho - (1.0_R8 - fl%lambda) / fl%rho_0)
         case default ; stop 'Bad case number pres'
      endselect
   endfunction pres

endmodule fluid_law

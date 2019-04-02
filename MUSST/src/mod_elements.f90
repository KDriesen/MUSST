
!! author: Noël Brunetière<br/>&emsp;Arthur Francisco
!! version: 1.0.0
!! date: March, 26 2019

!< <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<       **Subroutines about the finite elements that are used**
!< </span>

module elements
use data_arch, only : I4, R8, HIG_R8
implicit none

private

public :: ni4, ni4_up_2d, dj4, calc_ni4_xy_derivatives

contains


   !=========================================================================================
   !< @note Standard QU4 element
   !  @endnote
   !-----------------------------------------------------------------------------------------
   real(kind=R8) function ni4(numn, ksi, eta)
   implicit none
   real(kind=R8),    intent(in) :: ksi
   real(kind=R8),    intent(in) :: eta
   integer(kind=I4), intent(in) :: numn   !! *local node number: 1->SW, 2-> SE, 3->NE, 4->NW*
      ni4 = HIG_R8
      select case(numn)
         case(1); ni4 = (1._R8 -ksi)*(1._R8 -eta)/4._R8
         case(2); ni4 = (1._R8 +ksi)*(1._R8 -eta)/4._R8
         case(3); ni4 = (1._R8 +ksi)*(1._R8 +eta)/4._R8
         case(4); ni4 = (1._R8 -ksi)*(1._R8 +eta)/4._R8
         case default; stop 'ni4 bad case'
      endselect
   return
   endfunction ni4


   !=========================================================================================
   !< @note Up-Downwind 1D element for perfect gaz cases
   !
   !   Instead of the standard linear 1D element, a parameter *pe* is introduced to control the
   !   slope. The higher *pe* the steepest curve => it can even lead to a quasi-discontinuous function.
   !   If *pe* is zero, then the shape functions are the standard ones.
   !
   !   *se* is linked to the flow direction: se=1 for upwind, se=-1 for downwind.
   !
   !   <p style="text-align:center;"><img src="../media/up_downwind.png" alt="Upwinding coeffs" width="450px"/></p>
   !
   ! @endnote
   !-----------------------------------------------------------------------------------------
   real(kind=R8) function ni4_up_1d(numn, ksi, pe, se)
   implicit none
   integer(kind=I4), intent(in) :: numn   !! *local node number*
   real(kind=R8),    intent(in) :: ksi
   real(kind=R8),    intent(in) :: pe     !! *parameter controlling the shape function slope*
   real(kind=R8),    intent(in) :: se     !! *+1 -> upwind, -1 -> downwind*
      real(kind=R8) :: p, t
      p = (pe +1._R8)/2
      t = p*(ksi -se)
      select case(numn)
         case(1); ni4_up_1d = -t + (1._R8 -se)/2
         case(2); ni4_up_1d = +t + (1._R8 +se)/2
         case default; stop 'ni4_up_1d bad case'
      endselect
      ni4_up_1d = max( ni4_up_1d, 0._R8)
      ni4_up_1d = min( ni4_up_1d, 1._R8)
   return
   endfunction ni4_up_1d


   !=========================================================================================
   !< @note Up-Downwind 2D element for perfect gaz cases
   !
   !  @endnote
   !-----------------------------------------------------------------------------------------
   real(kind=R8) function ni4_up_2d(numn, ksi, eta, pe, se)
   implicit none
   integer(kind=I4), intent(in)               :: numn !! *local node number*
   real(kind=R8),    intent(in)               :: ksi
   real(kind=R8),    intent(in)               :: eta
   real(kind=R8),    intent(in), dimension(2) :: pe   !! *parameter controling the shape function slope for each direction*
   real(kind=R8),    intent(in), dimension(2) :: se   !! *+1 -> upwind, -1 -> downwind, for each direction*
      select case(numn)
         case(1); ni4_up_2d = ni4_up_1d(1, ksi, pe(1), se(1)) * ni4_up_1d(1, eta, pe(2), se(2))
         case(2); ni4_up_2d = ni4_up_1d(2, ksi, pe(1), se(1)) * ni4_up_1d(1, eta, pe(2), se(2))
         case(3); ni4_up_2d = ni4_up_1d(2, ksi, pe(1), se(1)) * ni4_up_1d(2, eta, pe(2), se(2))
         case(4); ni4_up_2d = ni4_up_1d(1, ksi, pe(1), se(1)) * ni4_up_1d(2, eta, pe(2), se(2))
         case default; stop 'ni4_up_2d bad case'
      endselect
   return
   endfunction ni4_up_2d


   !=========================================================================================
   !< @note Up-Downwind 1D element for perfect gaz cases, smooth functions
   !  @endnote
   !  @warning
   !    NOT USED
   !  @endwarning
   !-----------------------------------------------------------------------------------------
   real(kind=R8) function ni_tanh_1d(numn, ksi, pe, se)
   implicit none
   integer(kind=I4), intent(in) :: numn
   real(kind=R8),    intent(in) :: ksi
   real(kind=R8),    intent(in) :: pe
   real(kind=R8),    intent(in) :: se
      real(kind=R8) :: tksi, x
      if (se>0.) then
         x = 1._r8 -ksi
      else
         x = 1._r8 +ksi
      endif
      tksi = tanh(pe*x)/tanh(2*pe)

      select case(numn)
         case(1); ni_tanh_1d = ((1._R8 - se)/2 +se*tksi)
         case(2); ni_tanh_1d = ((1._R8 + se)/2 -se*tksi)
         case default; stop 'ni_tanh_1d bad case'
      endselect
   return
   endfunction ni_tanh_1d


   !=========================================================================================
   !< @note Up-Downwind 1D element for perfect gaz cases, smooth functions
   !  @endnote
   !  @warning
   !    NOT USED
   !  @endwarning
   !-----------------------------------------------------------------------------------------
   real(kind=R8) function ni_tanh_1d_der(numn, ksi, pe, se)
   implicit none
   integer(kind=I4), intent(in) :: numn
   real(kind=R8),    intent(in) :: ksi
   real(kind=R8),    intent(in) :: pe
   real(kind=R8),    intent(in) :: se
      real(kind=R8) :: tksi_der, tksi1, tksi2, x
      if (se>0.) then
         x = 1._r8 -ksi
      else
         x = 1._r8 +ksi
      endif
      tksi1 = tanh(pe*x)
      tksi2 = tanh(2*pe)
      tksi_der = (pe/tksi2)*( 1._R8 -tksi1**2 )

      select case(numn)
         case(1); ni_tanh_1d_der = ((1._R8 - se)/2 +se*tksi_der)
         case(2); ni_tanh_1d_der = ((1._R8 + se)/2 -se*tksi_der)
         case default; stop 'ni_tanh_1d_der bad case'
      endselect
   return
   endfunction ni_tanh_1d_der


   !=========================================================================================
   !< @note Up-Downwind 2D element for perfect gaz cases, smooth functions
   !  @endnote
   !  @warning
   !    NOT USED
   !  @endwarning
   !-----------------------------------------------------------------------------------------
   real(kind=R8) function ni_tanh_2d(numn, ksi, eta, pe, se)
   implicit none
   integer(kind=I4), intent(in)               :: numn
   real(kind=R8),    intent(in)               :: ksi
   real(kind=R8),    intent(in)               :: eta
   real(kind=R8),    intent(in), dimension(2) :: pe
   real(kind=R8),    intent(in), dimension(2) :: se
      select case(numn)
         case(1); ni_tanh_2d = ni_tanh_1d(1, ksi, pe(1), se(1))*ni_tanh_1d(1, eta, pe(2), se(2))
         case(2); ni_tanh_2d = ni_tanh_1d(2, ksi, pe(1), se(1))*ni_tanh_1d(1, eta, pe(2), se(2))
         case(3); ni_tanh_2d = ni_tanh_1d(2, ksi, pe(1), se(1))*ni_tanh_1d(2, eta, pe(2), se(2))
         case(4); ni_tanh_2d = ni_tanh_1d(1, ksi, pe(1), se(1))*ni_tanh_1d(2, eta, pe(2), se(2))
         case default; stop 'ni_tanh_2d bad case'
      endselect
   return
   endfunction ni_tanh_2d


   !=========================================================================================
   !< @note Up-Downwind 2D element for perfect gaz cases, smooth functions
   !  @endnote
   !  @warning
   !    NOT USED
   !  @endwarning
   !-----------------------------------------------------------------------------------------
   real(kind=R8) function ni_tanh_2d_ksi(numn, ksi, eta, pe, se)
   implicit none
   integer(kind=I4), intent(in)               :: numn
   real(kind=R8),    intent(in)               :: ksi
   real(kind=R8),    intent(in)               :: eta
   real(kind=R8),    intent(in), dimension(2) :: pe
   real(kind=R8),    intent(in), dimension(2) :: se
      select case(numn)
         case(1); ni_tanh_2d_ksi = ni_tanh_1d_der(1, ksi, pe(1), se(1))*ni_tanh_1d(1, eta, pe(2), se(2))
         case(2); ni_tanh_2d_ksi = ni_tanh_1d_der(2, ksi, pe(1), se(1))*ni_tanh_1d(1, eta, pe(2), se(2))
         case(3); ni_tanh_2d_ksi = ni_tanh_1d_der(2, ksi, pe(1), se(1))*ni_tanh_1d(2, eta, pe(2), se(2))
         case(4); ni_tanh_2d_ksi = ni_tanh_1d_der(1, ksi, pe(1), se(1))*ni_tanh_1d(2, eta, pe(2), se(2))
         case default; stop 'ni_tanh_2d bad case'
      endselect
   return
   endfunction ni_tanh_2d_ksi


   !=========================================================================================
   !< @note Up-Downwind 2D element for perfect gaz cases, smooth functions
   !  @endnote
   !  @warning
   !    NOT USED
   !  @endwarning
   !-----------------------------------------------------------------------------------------
   real(kind=R8) function ni_tanh_2d_eta(numn, ksi, eta, pe, se)
   implicit none
   integer(kind=I4), intent(in)               :: numn
   real(kind=R8),    intent(in)               :: ksi
   real(kind=R8),    intent(in)               :: eta
   real(kind=R8),    intent(in), dimension(2) :: pe
   real(kind=R8),    intent(in), dimension(2) :: se
      select case(numn)
         case(1); ni_tanh_2d_eta = ni_tanh_1d(1, ksi, pe(1), se(1))*ni_tanh_1d_der(1, eta, pe(2), se(2))
         case(2); ni_tanh_2d_eta = ni_tanh_1d(2, ksi, pe(1), se(1))*ni_tanh_1d_der(1, eta, pe(2), se(2))
         case(3); ni_tanh_2d_eta = ni_tanh_1d(2, ksi, pe(1), se(1))*ni_tanh_1d_der(2, eta, pe(2), se(2))
         case(4); ni_tanh_2d_eta = ni_tanh_1d(1, ksi, pe(1), se(1))*ni_tanh_1d_der(2, eta, pe(2), se(2))
         case default; stop 'ni_tanh_2d_eta bad case'
      endselect
   return
   endfunction ni_tanh_2d_eta


   !=========================================================================================
   !< @note Standard 2D element function derivatives (with respect to \(\xi\))
   !  @endnote
   !-----------------------------------------------------------------------------------------
   real(kind=R8) function ni4ksi(numn, eta)
   implicit none
   integer(kind=I4), intent(in) :: numn
   real(kind=R8),    intent(in) :: eta
      select case(numn)
         case(1); ni4ksi = -(1._R8 -eta)/4._R8
         case(2); ni4ksi = +(1._R8 -eta)/4._R8
         case(3); ni4ksi = +(1._R8 +eta)/4._R8
         case(4); ni4ksi = -(1._R8 +eta)/4._R8
         case default; stop 'ni4ksi bad case'
      endselect
   return
   endfunction ni4ksi


   !=========================================================================================
   !< @note Standard 2D element function derivatives (with respect to \(\eta\))
   !  @endnote
   !-----------------------------------------------------------------------------------------
   real(kind=R8) function ni4eta(numn, ksi)
   implicit none
   integer(kind=I4), intent(in) :: numn
   real(kind=R8),    intent(in) :: ksi
      ni4eta = HIG_R8
      select case(numn)
         case(1); ni4eta = -(1._R8 -ksi)/4._R8
         case(2); ni4eta = -(1._R8 +ksi)/4._R8
         case(3); ni4eta = +(1._R8 +ksi)/4._R8
         case(4); ni4eta = +(1._R8 -ksi)/4._R8
         case default; stop 'ni4eta bad case'
      endselect
   return
   endfunction ni4eta


   !=========================================================================================
   !< @note Jacobian intermediate calculus
   !
   !   *Example*
   !   if numn=1 \(j_4=\sum_{i=1}^{4} \frac{\partial N_i}{\partial \xi} x_i\)
   !  @endnote
   !-----------------------------------------------------------------------------------------
   real(kind=R8) function j4(numn, ksi, eta, x, y)
   implicit none
   integer(kind=I4), intent(in)               :: numn
   real(kind=R8),    intent(in)               :: ksi
   real(kind=R8),    intent(in)               :: eta
   real(kind=R8),    intent(in), dimension(4) :: x
   real(kind=R8),    intent(in), dimension(4) :: y
      real(kind=R8)    :: valint
      integer(kind=I4) :: i
      valint = 0._R8
      select case(numn)
         case(1)
            do i = 1, 4
               valint = valint +ni4ksi(i, eta)*x(i)
            enddo
         case(2)
            do i = 1, 4
               valint = valint +ni4ksi(i, eta)*y(i)
            enddo
         case(3)
            do i = 1, 4
               valint = valint +ni4eta(i, ksi)*x(i)
            enddo
         case(4)
            do i = 1, 4
               valint = valint +ni4eta(i, ksi)*y(i)
            enddo
         case default
            stop 'j4 bad case'
      endselect
      j4 = valint
   return
   endfunction j4


   !=========================================================================================
   !< @note Jacobian calculated without calling any function, for computing speed reasons.
   !
   !  @endnote
   !< @note **Maxima**, **Wxmaxima** script
   !
   !   ```ni4_1(xi, eta):= (1. -xi)*(1. -eta)/4.$```                                                                                <br/>
   !   ```ni4_2(xi, eta):= (1. +xi)*(1. -eta)/4.$```                                                                                <br/>
   !   ```ni4_3(xi, eta):= (1. +xi)*(1. +eta)/4.$```                                                                                <br/>
   !   ```ni4_4(xi, eta):= (1. -xi)*(1. +eta)/4.$```                                                                                <br/>
   !   ```ni4xi_1(xi, eta):=diff(ni4_1(xi, eta),xi,1)$```                                                                           <br/>
   !   ```ni4xi_2(xi, eta):=diff(ni4_2(xi, eta),xi,1)$```                                                                           <br/>
   !   ```ni4xi_3(xi, eta):=diff(ni4_3(xi, eta),xi,1)$```                                                                           <br/>
   !   ```ni4xi_4(xi, eta):=diff(ni4_4(xi, eta),xi,1)$```                                                                           <br/>
   !   ```ni4eta_1(xi, eta):=diff(ni4_1(xi, eta),eta,1)$```                                                                         <br/>
   !   ```ni4eta_2(xi, eta):=diff(ni4_2(xi, eta),eta,1)$```                                                                         <br/>
   !   ```ni4eta_3(xi, eta):=diff(ni4_3(xi, eta),eta,1)$```                                                                         <br/>
   !   ```ni4eta_4(xi, eta):=diff(ni4_4(xi, eta),eta,1)$```                                                                         <br/>
   !   ```j4_1(xi, eta):=factorsum(ni4xi_1(xi, eta)*x1+ni4xi_2(xi, eta)*x2+ni4xi_3(xi, eta)*x3+ni4xi_4(xi, eta)*x4)$```             <br/>
   !   ```j4_1(xi, eta);```                                                                                                         <br/>
   !   ```j4_2(xi, eta):=factorsum(ni4xi_1(xi, eta)*y1+ni4xi_2(xi, eta)*y2+ni4xi_3(xi, eta)*y3+ni4xi_4(xi, eta)*y4)$```             <br/>
   !   ```j4_2(xi, eta);```                                                                                                         <br/>
   !   ```j4_3(xi, eta):=factorsum(ni4eta_1(xi, eta)*x1+ni4eta_2(xi, eta)*x2+ni4eta_3(xi, eta)*x3+ni4eta_4(xi, eta)*x4)$```         <br/>
   !   ```j4_3(xi, eta);```                                                                                                         <br/>
   !   ```j4_4(xi, eta):=factorsum(ni4eta_1(xi, eta)*y1+ni4eta_2(xi, eta)*y2+ni4eta_3(xi, eta)*y3+ni4eta_4(xi, eta)*y4)$```         <br/>
   !   ```j4_4(xi, eta);```                                                                                                         <br/>
   !   ```dj4(xi, eta):=factorout(j4_1(xi, eta)*j4_4(xi, eta)-j4_2(xi, eta)*j4_3(xi, eta),xi,eta)$```                               <br/>
   !   ```dj4(xi, eta);```                                                                                                          <br/>
   !   ```ni4x_1(xi,eta):=factorout(j4_4(xi, eta)*ni4xi_1(xi, eta),eta,xi) - factorout(j4_3(xi, eta)*ni4eta_1(xi, eta),xi,eta)$```  <br/>
   !   ```ni4x_1(xi,eta);```                                                                                                        <br/>
   !   ```ni4x_2(xi,eta):=factorout(j4_4(xi, eta)*ni4xi_2(xi, eta),eta,xi) - factorout(j4_3(xi, eta)*ni4eta_2(xi, eta),eta,xi)$```  <br/>
   !   ```ni4x_2(xi,eta);```                                                                                                        <br/>
   !   ```ni4x_3(xi,eta):=factorout(j4_4(xi, eta)*ni4xi_3(xi, eta),eta,xi) - factorout(j4_3(xi, eta)*ni4eta_3(xi, eta),xi,eta)$```  <br/>
   !   ```ni4x_3(xi,eta);```                                                                                                        <br/>
   !   ```ni4x_4(xi,eta):=factorout(j4_4(xi, eta)*ni4xi_4(xi, eta),eta,xi) - factorout(j4_3(xi, eta)*ni4eta_4(xi, eta),eta,xi)$```  <br/>
   !   ```ni4x_4(xi,eta);```                                                                                                        <br/>
   !   ```ni4y_1(xi,eta):=factorout(-j4_2(xi, eta)*ni4xi_1(xi, eta),eta,xi) + factorout(j4_1(xi, eta)*ni4eta_1(xi, eta),xi,eta)$``` <br/>
   !   ```ni4y_1(xi,eta);```                                                                                                        <br/>
   !   ```ni4y_2(xi,eta):=factorout(-j4_2(xi, eta)*ni4xi_2(xi, eta),eta,xi) + factorout(j4_1(xi, eta)*ni4eta_2(xi, eta),eta,xi)$``` <br/>
   !   ```ni4y_2(xi,eta);```                                                                                                        <br/>
   !   ```ni4y_3(xi,eta):=factorout(-j4_2(xi, eta)*ni4xi_3(xi, eta),eta,xi) + factorout(j4_1(xi, eta)*ni4eta_3(xi, eta),xi,eta)$``` <br/>
   !   ```ni4y_3(xi,eta);```                                                                                                        <br/>
   !   ```ni4y_4(xi,eta):=factorout(-j4_2(xi, eta)*ni4xi_4(xi, eta),eta,xi) + factorout(j4_1(xi, eta)*ni4eta_4(xi, eta),eta,xi)$``` <br/>
   !   ```ni4y_4(xi,eta);```                                                                                                        <br/>
   !
   !  @endnote
   !-----------------------------------------------------------------------------------------
   real(kind=R8) function dj4(ksi, eta, x, y)
   implicit none
   real(kind=R8), intent(in)               :: ksi
   real(kind=R8), intent(in)               :: eta
   real(kind=R8), intent(in), dimension(4) :: x
   real(kind=R8), intent(in), dimension(4) :: y

   !~       standard calculation
   !~       dj4 = j4(1, ksi, eta, x, y) * j4(4, ksi, eta, x, y) - & !
   !~             j4(2, ksi, eta, x, y) * j4(3, ksi, eta, x, y)

      real(kind=R8) :: km1, kp1, em1, ep1, kme, kpe, x1, x2, x3, x4, y1, y2, y3, y4

      km1 = ksi -1._R8
      kp1 = ksi +1._R8
      em1 = eta -1._R8
      ep1 = eta +1._R8
      kme = ksi -eta
      kpe = ksi +eta
      x1 = x(1) ; x2 = x(2) ; x3 = x(3) ; x4 = x(4)
      y1 = y(1) ; y2 = y(2) ; y3 = y(3) ; y4 = y(4)

      dj4 = 0.125_R8*( +km1*(x1*y4 -x4*y1) & !
                       +kp1*(x2*y3 -x3*y2) & !
                       -em1*(x1*y2 -x2*y1) & !
                       +ep1*(x3*y4 -x4*y3) & !
                       -kpe*(x2*y4 -x4*y2) & !
                       -kme*(x1*y3 -x3*y1) )

   return
   endfunction dj4


   !=========================================================================================
   !< @note QU4 derivatives with respect to \(x\) and \(y\) without any function call, for speed reasons
   !<
   !-----------------------------------------------------------------------------------------
   subroutine calc_ni4_xy_derivatives(ni4x, ni4y, ksi, eta, x, y, dj)
   implicit none
   real(kind=R8),    intent(out), dimension(4) :: ni4x
   real(kind=R8),    intent(out), dimension(4) :: ni4y
   real(kind=R8),    intent(in )               :: ksi
   real(kind=R8),    intent(in )               :: eta
   real(kind=R8),    intent(in ), dimension(4) :: x
   real(kind=R8),    intent(in ), dimension(4) :: y
   real(kind=R8),    intent(in )               :: dj

   !~       standard calculation
   !~       ni4x = ( j4(4, ksi, eta, x, y)*ni4ksi(numn, eta) -&
   !~                j4(3, ksi, eta, x, y)*ni4eta(numn, ksi) )/dj4(ksi, eta, x, y)
   !~       ni4y = ( -j4(2, ksi, eta, x, y)*ni4ksi(numn, eta) +&
   !~                 j4(1, ksi, eta, x, y)*ni4eta(numn, ksi) )/dj4(ksi, eta, x, y)

      real(kind=R8)    :: km1, kp1, em1, ep1, kme, kpe, & !
                          x1, x2, x3, x4,               & !
                          y1, y2, y3, y4,               & !
                          kx, ky, ex, ey

      km1 = ksi -1._R8
      kp1 = ksi +1._R8
      em1 = eta -1._R8
      ep1 = eta +1._R8
      kme = ksi -eta
      kpe = ksi +eta

      x1 = x(1) ; x2 = x(2) ; x3 = x(3) ; x4 = x(4)
      y1 = y(1) ; y2 = y(2) ; y3 = y(3) ; y4 = y(4)

      kx = -km1*x1 +kp1*x2 -kp1*x3 +km1*x4
      ky = -km1*y1 +kp1*y2 -kp1*y3 +km1*y4
      ex = -em1*x1 +ep1*x2 -ep1*x3 +em1*x4
      ey = -em1*y1 +ep1*y2 -ep1*y3 +em1*y4

      ni4x(1) = 0.0625_R8*(+km1*kx -em1*ky)/dj
      ni4x(2) = 0.0625_R8*(-kp1*kx +em1*ky)/dj
      ni4x(3) = 0.0625_R8*(+kp1*kx -ep1*ky)/dj
      ni4x(4) = 0.0625_R8*(-km1*kx +ep1*ky)/dj

      kx = -km1*x1 +km1*x2 -kp1*x3 +kp1*x4
      ky = -km1*y1 +km1*y2 -kp1*y3 +kp1*y4
      ex = -em1*x1 +em1*x2 -ep1*x3 +ep1*x4
      ey = -em1*y1 +em1*y2 -ep1*y3 +ep1*y4

      ni4y(1) = 0.0625_R8*(-km1*ex +em1*ey)/dj
      ni4y(2) = 0.0625_R8*(+kp1*ex -em1*ey)/dj
      ni4y(3) = 0.0625_R8*(-kp1*ex +ep1*ey)/dj
      ni4y(4) = 0.0625_R8*(+km1*ex -ep1*ey)/dj

   return
   endsubroutine calc_ni4_xy_derivatives


endmodule elements

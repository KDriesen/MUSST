
!! author: Noël Brunetière<br/>&emsp;Arthur Francisco
!! version: 1.0.0
!! date: April,17 2017
!! summary: save subroutines

!< <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!         **Subroutines to save the data**
!  </span>

module inout_files
use VTK
use data_arch, only : I4, R8, get_unit
use film
use ms_film
use surfile
implicit none

private

integer(kind=I4)  :: s_vtk                                           ! flag for vtk save
integer(kind=I4), parameter :: NO_VTK = 0, TS_VTK = 1, BS_VTK = 2    ! parameters for s_vtk

public :: s_vtk, NO_VTK, TS_VTK, BS_VTK, save_fe_f_vtk,      &
                                         save_ms_fe_f_vtk,   &
                                         save_profile_x_comp_slider,       &
                                         save_profile_x_comp_air_pocket,   &
                                         save_profile_y_comp_air_pocket,   &
                                         save_profile_x_fe,                &
                                         save_profile_x_ms
contains


   !=========================================================================================
   !< @note Subroutine to save deterministic data following VTK model
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine save_fe_f_vtk(fe_f, nom_fic)
   implicit none
   type(FE_FILM),    intent(in) :: fe_f      !! [[FE_FILM]] *element to store*
   character(len=*), intent(in) :: nom_fic   !! *filename*

      integer(kind=I4) :: E_IO   ! parametre identification vtk
      integer(kind=I4) :: i, k

      real(kind=R8),    dimension(:), allocatable :: z            ! tableau supplementaire pour la sortie
      integer(kind=I4), dimension(:), allocatable :: connec, tipo

      allocate (z(fe_f%m%n), connec(fe_f%m%ne*(1+4)), tipo(fe_f%m%ne))
      z    = 0._R8
      tipo = 9       ! elements a quatre noeuds: 9 dans le formalisme vtk

      k = 0
      do i = 1, fe_f%m%ne
         connec(k+1) = 4
         connec(k+2) = fe_f%m%con(i, 1) -1
         connec(k+3) = fe_f%m%con(i, 2) -1
         connec(k+4) = fe_f%m%con(i, 3) -1
         connec(k+5) = fe_f%m%con(i, 4) -1
         k=k+5
      enddo
      ! ouverture du fichier
      E_IO = VTK_INI(output_format = 'binary', &
                          filename = nom_fic, &
                             title = 'Resultats fe_f fluide', &
                     mesh_topology = 'UNSTRUCTURED_GRID')
      ! ecriture des coordonnees
      E_IO = VTK_GEO(NN = fe_f%m%n, &
                      X = fe_f%m%x, &
                      Y = fe_f%m%y, &
                      Z = Z)
      ! definition des cellules
      E_IO = VTK_CON(NC      = fe_f%m%ne, &
                     connect = connec,    &
                   cell_type = tipo)
      ! choix du type de variable (elementaire)
      E_IO = VTK_DAT(NC_NN   = fe_f%m%ne, &
                var_location = 'cell')
      ! passage des variables cellules
      do i = 1, fe_f%n_vc
          E_IO = VTK_VAR(NC_NN = fe_f%m%ne,       &
                       varname = fe_f%vc_name(i), &
                           var = fe_f%vc(:,i))
      enddo
      ! choix du type de variable (nodale)
      E_IO = VTK_DAT(NC_NN   = fe_f%m%n, &
                var_location = 'node')
      ! passage des variables nodales
      do i = 1, fe_f%n_vn
          E_IO = VTK_VAR(NC_NN = fe_f%m%n,        &
                       varname = fe_f%vn_name(i), &
                           var = fe_f%vn(:, i))
      enddo
      ! fermeture du fichier vtk
      E_IO = VTK_END()
      ! liberation des tableaux
      deallocate(z, connec, tipo)
   return
   endsubroutine save_fe_f_vtk


   !=========================================================================================
   !< @note Subroutine to save the top scale data following VTK model
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine save_ms_fe_f_vtk(ms_fe_f, nom_fic)
   implicit none
   type(MS_FE_FILM), intent(in) :: ms_fe_f   !! [[MS_FE_FILM]] *element to store*
   character(len=*), intent(in) :: nom_fic   !! *filename*

      integer(kind=I4)   :: e, ne, i, l
      character(len=256) :: nom_fic_bs, nom_fic_ts
      character(len=128) :: suff
      nom_fic_bs = repeat(" ", len(nom_fic_bs))
      nom_fic_ts = repeat(" ", len(nom_fic_bs))
      suff       = repeat(" ", len(suff)      )
      i = index(trim(nom_fic), "/", back=.true.)
      l = len_trim(nom_fic)

      if (s_vtk > NO_VTK) then
         nom_fic_ts = nom_fic(1:i)//"ts_"//nom_fic(i+1:l)
         call save_fe_f_vtk(ms_fe_f%ts_fe_f, trim(nom_fic_ts))
         write(*,*) trim(nom_fic_ts)
      endif

      if (s_vtk == BS_VTK) then
         ne = ms_fe_f%ts_fe_f%m%ne
         do e = 1, ne
            write(suff,'(i5.5, 2a)') e, "_", nom_fic(i+1:l)
            nom_fic_bs = nom_fic(1:i)//trim(suff)
            call save_fe_f_vtk(ms_fe_f%bs_fe_f(e), trim(nom_fic_bs))
         enddo
      endif

   return
   endsubroutine save_ms_fe_f_vtk


   !=========================================================================================
   !< @note Subroutine to save the pressures along a line following the flow, at distance
   !        \(zy\) from a border.
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine save_profile_x_fe(fe_f, file_name, lx, zy)
   implicit none
   type(FE_FILM),    intent(in) :: fe_f         !! [[FE_FILM]] *element*
   character(len=*), intent(in) :: file_name    !! *filename*
   real(kind=R8),    intent(in) :: lx           !! *surface length*
   real(kind=R8),    intent(in) :: zy           !! *distance from a border*

      real(kind=R8)    :: dr
      integer(kind=I4) :: i, k

      dr = lx / fe_f%m%n
      dr = dr / 10
      call get_unit(k)
      open(k, file = file_name, status = 'unknown')
      do i = 1, fe_f%m%n
         if (((fe_f%m%y(i) - zy)**2) < (dr ** 2)) then
            write(k, *) fe_f%m%x(i), fe_f%vn(i, H_N), fe_f%vn(i, P_N)
         endif
      enddo
      close(k)

   return
   endsubroutine save_profile_x_fe


   !=========================================================================================
   !< @note Subroutine to save the pressures along a line following the flow, at distance
   !        \(zy\) from a border.
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine save_profile_x_ms(ms_fe_f, file_name, lx, zy)
   implicit none
   type(MS_FE_FILM), intent(in) :: ms_fe_f      !! [[MS_FE_FILM]] *element*
   character(len=*), intent(in) :: file_name    !! *filename*
   real(kind=R8),    intent(in) :: lx           !! *surface length*
   real(kind=R8),    intent(in) :: zy           !! *distance from a border*

      real(kind=R8)    :: dr
      integer(kind=I4) :: i, e, k

      dr = lx / ms_fe_f%bs_fe_f(1)%m%n
      dr = dr / 10
      call get_unit(k)
      open(k, file = file_name, status = 'unknown')
         do e = 1, ms_fe_f%ts_fe_f%m%ne
         do i = 1, ms_fe_f%bs_fe_f(e)%m%n
            if (((ms_fe_f%bs_fe_f(e)%m%y(i) - zy)**2) < (dr ** 2)) then
               write (k, *) ms_fe_f%bs_fe_f(e)%m%x(i), ms_fe_f%bs_fe_f(e)%vn(i, H_N), &
                                                       ms_fe_f%bs_fe_f(e)%vn(i, P_N)
            endif
         enddo
         enddo
      close(k)

   return
   endsubroutine save_profile_x_ms


   !=========================================================================================
   !< @note Subroutine to save the pressures along a line following the flow, at distance
   !        \(zy\) from a border. The theoretical values are also stored.
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine save_profile_x_comp_slider(fe_f, file_name, lx, zy)
   implicit none
   type(FE_FILM),    intent(in) :: fe_f         !! [[FE_FILM]] *element*
   character(len=*), intent(in) :: file_name    !! *filename*
   real(kind=R8),    intent(in) :: lx           !! *slider length*
   real(kind=R8),    intent(in) :: zy           !! *distance from a border*

      real(kind=R8)    :: pref, K, dr, hr, hb, p0
      integer(kind=I4) :: i, kk

      dr = lx / fe_f%m%n
      dr = dr / 10
      K  = maxval(fe_f%vn(:, H_N))/minval(fe_f%vn(:, H_N))
      hr = minval(fe_f%vn(:, H_N))
      pref = 6 * fe_f%data_f%fl%mu_0 * fe_f%data_f%V_x * lx / (hr ** 2)
      p0 = minval(fe_f%vn(:, P_N))
      call get_unit(kk)
      open(kk, file = file_name, status = 'unknown')
      do i = 1, fe_f%m%n
         if (((fe_f%m%y(i) - zy)**2) < (dr ** 2)) then
            hb = fe_f%vn(i,H_N)/hr
            write (kk, *) fe_f%m%x(i)/lx, hb, (fe_f%vn(i, P_N) - p0) / pref, &
                          1.0 / (K - 1) * (1.0 / hb - K / (K + 1) * (hb ** (-2) - K ** (-2)) - 1.0 / K)
         endif
      enddo
      close(kk)

   return
   endsubroutine save_profile_x_comp_slider


   !=========================================================================================
   !< @note Subroutine to save the pressures along a line following the flow, at distance
   !        \(zy\) from a border
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine save_profile_x_comp_air_pocket(fe_f, file_name, lx, zy, bc)
   implicit none
   type(FE_FILM),    intent(in)               :: fe_f          !! [[FE_FILM]] *element*
   character(len=*), intent(in)               :: file_name     !! *filename*
   real(kind=R8),    intent(in)               :: lx            !! *pocket length*
   real(kind=R8),    intent(in)               :: zy            !! *distance from a border*
   real(kind=R8),    intent(in), dimension(4) :: bc            !! *pressure boundaries*

      real(kind=R8)    :: dr, hr, hb, p0
      integer(kind=I4) :: i, k

      dr = lx / fe_f%m%n
      dr = dr / 10
      hr = minval(fe_f%vn(:, H_N))
      p0 = bc(1)
      call get_unit(k)
      open(k, file = file_name, status = 'unknown')
      do i = 1, fe_f%m%n
         if (((fe_f%m%y(i) - zy)**2) < (dr ** 2)) then
            hb = (fe_f%vn(i, H_N) + fe_f%vc(i, HG_C))/hr
            write (k, *) fe_f%m%x(i)/lx, hb, fe_f%vn(i, P_N) / p0, hb
         endif
      enddo
      close(k)

   return
   endsubroutine save_profile_x_comp_air_pocket


   !=========================================================================================
   !< @note Subroutine to save the pressures along a line perpendicular to the flow, at distance
   !        \(zx\)
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine save_profile_y_comp_air_pocket(fe_f, file_name, ly, zx, bc)
   implicit none
   type(FE_FILM),    intent(in)               :: fe_f          !! [[FE_FILM]] *element*
   character(len=*), intent(in)               :: file_name     !! *filename*
   real(kind=R8),    intent(in)               :: ly            !! *pocket width*
   real(kind=R8),    intent(in)               :: zx            !! *distance from the pocket entry*
   real(kind=R8),    intent(in), dimension(4) :: bc            !! *pressure boundaries*

      real(kind=R8)    :: dr, hr, hb, p0
      integer(kind=I4) :: i, ii, jj, k
      logical(kind=I4), allocatable, dimension(:) :: done

      dr = ly / fe_f%m%n
      dr = dr / 10
      hr = minval(fe_f%vn(:, H_N))
      p0 = bc(1)
      allocate(done(1:fe_f%m%n)) ; done = .false.
      call get_unit(k)
      open(k, file = file_name, status = 'unknown')
      do ii = 1, fe_f%m%ne
      do jj = 1, fe_f%m%el_t(ii)
         i = fe_f%m%con(ii, jj)
         if (done(i)) cycle
         done(i) = .true.
         if (((fe_f%m%x(i) - zx)**2) < (dr ** 2)) then
            hb = (fe_f%vn(i, H_N) + fe_f%vc(ii, HG_C))/hr
            write (k, *) fe_f%m%y(i)/ly, hb, fe_f%vn(i, P_N) / p0, hb
         endif
      enddo
      enddo
      close(k)
      deallocate(done)

   return
   endsubroutine save_profile_y_comp_air_pocket

endmodule inout_files

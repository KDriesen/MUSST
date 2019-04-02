
!! author: Noël Brunetière<br/>&emsp;Arthur Francisco
!! version: 1.0.0
!! date: April,17 2017
!! summary: run different tests

!< <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<       **Subroutines to read a batch file and run tests**
!< </span>

module test_musst
use data_arch, only : R4, PI_R8, get_unit
use data_film_hd
use ms_film
use film
use inout_files
use surfile
use num_param
use solver
implicit none

private

type(FE_FILM)    :: fe_f      !! [[FE_FILM]] *deterministic finite element*
type(MS_FE_FILM) :: ms_fe_f   !! [[MS_FE_FILM]] *multiscale finite element*
type(DATA_FILM)  :: data_f    !! [[DATA_FILM]] *parameter*
type(NUM_PAR)    :: num_pbs   !! [[NUM_PAR]] *parameter for the bottom-scale*
type(NUM_PAR)    :: num_pts   !! [[NUM_PAR]] *parameter for the top-scale*

type(MAT_SOLV)    :: mat      !! [[MAT_SOLV]] *solver type matrices*
type(MS_MAT_SOLV) :: ms_mat   !! [[MS_MAT_SOLV]] *solver type matrices*


integer(kind=I4)                            :: nx        !! *total number of nodes in \(x\) direction*
integer(kind=I4)                            :: ny        !! *total number of nodes in \(y\) direction*
integer(kind=I4)                            :: n_mac     !! *number of macro elements in a direction*
integer(kind=I4)                            :: n_mic     !! *number of nodes in \(x\) or \(y\) direction for the bottom scale*
real(kind=R8)                               :: lx        !! *domain size along \(x\)*
real(kind=R8)                               :: ly        !! *domain size along \(y\)*
real(kind=R8)                               :: sq        !! *roughness height*
real(kind=R8), dimension(:, :), allocatable :: tab_s     !! *roughness table*

real(kind=R8), dimension(4) :: bc, bf                 !! *boundary conditions*

real(kind=R4)     :: t1, t2                           !! *cpu time*
integer(kind=I4)  :: cend, cr, cinit                  !! *real time*
integer(kind=I4)  :: unit_num_res                     !! *file number*

character(len=256):: ms_vtk                           !! *output vtk file name*
character(len=256):: prof_ts, prof_bs                 !! *ts/bs mat profile name*
character(len=256):: res_file                         !! *result file name*
character(len=256):: surface_file                     !! *surface file name*
character(len= 15):: res_dir                          !! *"/out" subdirectory for results*

type(SCALE_SURF) :: scal_tmp                          !! *object [[SCALE_SURF]]*

integer(kind=I4) :: test_num                          !! *test number*
logical(kind=I4) :: save_PeK                          !! *save \(x\) Peclet field*
logical(kind=I4) :: save_PeE                          !! *save \(y\) Peclet field*

public :: run_test

contains


   !=========================================================================================
   !< @note Subroutine to read the data file ```.dat``` in the 'EXEC_MUSST' section and then
   !        to run the specified test
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine run_test(test, iunit, dir)
   implicit none
   integer(kind=I4),  intent(in) :: test   !! *number of the test to be performed*
   integer(kind=I4),  intent(in) :: iunit  !! *unit number of the data file*
   character(len=15), intent(in) :: dir    !! *output directory*

      test_num = test

      call read_data(iunit, dir)
      ! when ```run_test``` has been launched, *SOLVER_TS* has been determined
      mat%slv_t = SOLVER_TS

      call get_unit(unit_num_res) ; open(unit = unit_num_res, file = trim(res_file), status = 'unknown')

      select case(test)
         case( 1); call test_slider_fe
         case(11); call test_slider_ms

         case( 2); call test_bearing_x_fe
         case( 3); call test_bearing_y_fe

         case( 4); call test_rough_fe
         case(14); call test_rough_ms

         case( 5); call test_pocket_fe
         case default
      endselect

      close(unit_num_res)

   return
   endsubroutine run_test


   !=========================================================================================
   !< @note Subroutine to read the 'EXEC_MUSST' section
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine read_data(iunit, dir)
   implicit none
   integer(kind=I4),  intent(in) :: iunit
   character(len=15), intent(in) :: dir    !! *output directory*

      res_dir = dir

      read( iunit,*)

      read(iunit,*) data_f%h_0, data_f%h_g

      read(iunit,*) data_f%V_x, data_f%V_y

      read(iunit,*) data_f%pb_type
      read(iunit,*) data_f%fl%fluid_type
      read(iunit,*) data_f%fl%p_0
      read(iunit,*) data_f%fl%rho_0
      read(iunit,*) data_f%fl%mu_0
      read(iunit,*) data_f%fl%rg
      read(iunit,*) data_f%fl%lambda
      read(iunit,*) data_f%fl%T_0
      read(iunit,*)
      read(iunit,*) surface_file
      read(iunit,*) lx, ly

      read(iunit,*) nx, ny

      read(iunit,*) n_mac
      read(iunit,*) bc(1:4)

      read(iunit,*)
      read(iunit,*) num_pts

      read(iunit,*) num_pbs

      read(iunit,*) sq
      read(iunit,*) s_vtk, ms_vtk

      read(iunit,*) prof_ts
      read(iunit,*) prof_bs
      read(iunit,*) res_file

      ms_vtk   = "out/"//res_dir//"/"//trim(ms_vtk  )
      prof_ts  = "out/"//res_dir//"/"//trim(prof_ts )
      prof_bs  = "out/"//res_dir//"/"//trim(prof_bs )
      res_file = "out/"//res_dir//"/"//trim(res_file)

   return
   endsubroutine read_data


   !=========================================================================================
   !< @note Subroutine to define the slider geometry
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine modify_h_slider(fe_f, lx)
   implicit none
   type(FE_FILM), intent(inout) :: fe_f
   real(kind=R8), intent(in)    :: lx

      integer(kind=I4) :: i

      do i = 1, fe_f%m%n
         fe_f%vn(i, H2_N) = fe_f%data_f%h_0 + fe_f%data_f%h_0 * (lx - fe_f%m%x(i)) / lx
      enddo
      fe_f%vn(:, H_N) = fe_f%vn(:, H2_N) - fe_f%vn(:, H1_N)

   return
   endsubroutine modify_h_slider


   !=========================================================================================
   !< @note Subroutine to define the pocket geometry
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine modify_h_pocket(fe_f, lx)
   implicit none
   type(FE_FILM), intent(inout) :: fe_f
   real(kind=R8), intent(in)    :: lx

      integer(kind=I4) :: e, i
      logical(kind=I4) :: groove

      fe_f%vn(:,H2_N) = fe_f%data_f%h_0
      fe_f%vn(:,H1_N) = 0._R8
      do e = 1, fe_f%m%ne
         groove = .true.
         do i = 1, 4
            if (fe_f%m%x(fe_f%m%con(e,i)) > (lx / 2._R8)) groove = .false.
         enddo
         if (groove) then
            fe_f%vc(e, HG_C) = fe_f%data_f%h_g
         endif
      enddo
      fe_f%vn(:, H_N) = fe_f%vn(:, H2_N) - fe_f%vn(:, H1_N)

   return
   endsubroutine modify_h_pocket


   !=========================================================================================
   !< @note Subroutine to define the slider, in a multiscale problem
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine modify_h_slider_MS(ms_fe_f, lx)
   implicit none
   type(MS_FE_FILM), intent(inout) :: ms_fe_f
   real(kind=R8), intent(in)       :: lx

      integer(kind=I4) :: e

      call modify_h_slider(ms_fe_f%ts_fe_f, lx)
      do e = 1, ms_fe_f%ts_fe_f%m%ne
         call modify_h_slider(ms_fe_f%bs_fe_f(e), lx)
      enddo

   return
   endsubroutine modify_h_slider_MS


   !=========================================================================================
   !< @note Subroutine to define the bearing geometry
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine modify_h_bearing(fe_f, lx)
   implicit none
   type(FE_FILM), intent(inout) :: fe_f
   real(kind=R8), intent(in)    :: lx

      integer(kind=I4) :: i

      do i = 1, fe_f%m%n
         fe_f%vn(i, H2_N) = fe_f%data_f%h_0 + 0.5_R8 * fe_f%data_f%h_0 * cos(2 * PI_R8 * fe_f%m%x(i) / lx)
      enddo
      fe_f%vn(:, H_N) = fe_f%vn(:, H2_N) - fe_f%vn(:, H1_N)

   return
   endsubroutine modify_h_bearing


   !=========================================================================================
   !< @note Subroutine to define the bearing geometry
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine modify_h_bearing_y(fe_f, ly)
   implicit none
   type(FE_FILM), intent(inout) :: fe_f
   real(kind=R8), intent(in)    :: ly

      integer(kind=I4) :: i

      do i = 1, fe_f%m%n
         fe_f%vn(i, H2_N) = fe_f%data_f%h_0 + 0.5_R8 * fe_f%data_f%h_0 * cos(2 * PI_R8 * fe_f%m%y(i) / ly)
      enddo
      fe_f%vn(:, H_N) = fe_f%vn(:, H2_N) - fe_f%vn(:, H1_N)

   return
   endsubroutine modify_h_bearing_y


   !=========================================================================================
   !< @note Subroutine to apply a roughness table to a surface
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine apply_roughness(fe_f, tab_s)
   implicit none
   type(FE_FILM), intent(inout)                 :: fe_f
   real(kind=R8), intent(in   ), dimension(:,:) :: tab_s

      real(kind=R8)    :: xp, yp
      integer(kind=I4) :: k, i, j

      do k = 1, fe_f%m%n
         xp = fe_f%m%x(k)
         yp = fe_f%m%y(k)
         i = int((fe_f%m%nx - 1) * xp / fe_f%m%lx) + 1
         j = int((fe_f%m%ny - 1) * yp / fe_f%m%ly) + 1
         fe_f%vn(k, H2_N) = fe_f%data_f%h_0 - tab_s(i, j)
      enddo
      fe_f%vn(:, H_N) = fe_f%vn(:, H2_N) - fe_f%vn(:, H1_N)

   return
   endsubroutine apply_roughness


   !=========================================================================================
   !< @note Subroutine to apply a roughness table to a surface for a multiscale problem
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine apply_roughness_MS(ms_fe_f, tab_s)
   implicit none
   type(MS_FE_FILM), intent(inout)                 :: ms_fe_f
   real(kind=R8),    intent(in   ), dimension(:,:) :: tab_s

      real(kind=R8)    :: xp, yp
      integer(kind=I4) :: k, i, j, e

      do k = 1, ms_fe_f%ts_fe_f%m%n
         xp = ms_fe_f%ts_fe_f%m%x(k)
         yp = ms_fe_f%ts_fe_f%m%y(k)
         i = int((ms_fe_f%ts_fe_f%m%nx - 1) * xp / ms_fe_f%ts_fe_f%m%lx) + 1
         j = int((ms_fe_f%ts_fe_f%m%ny - 1) * yp / ms_fe_f%ts_fe_f%m%ly) + 1
         ms_fe_f%ts_fe_f%vn(k, H2_N) = ms_fe_f%ts_fe_f%data_f%h_0 - tab_s(i, j)
      enddo
      ms_fe_f%ts_fe_f%vn(:, H_N) = ms_fe_f%ts_fe_f%vn(:, H2_N) - ms_fe_f%ts_fe_f%vn(:, H1_N)

      do e = 1, ms_fe_f%ts_fe_f%m%ne
         call apply_roughness(ms_fe_f%bs_fe_f(e), tab_s)
      enddo

   return
   endsubroutine apply_roughness_ms


   !=========================================================================================
   !< @note Subroutine to run the slider test
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine test_slider_fe
   implicit none
      call init_fe_prob
      call modify_h_slider(fe_f, lx)
      call solve_fe_prob
   return
   endsubroutine test_slider_fe


   !=========================================================================================
   !< @note Subroutine to run the bearing test
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine test_bearing_x_fe
   implicit none
      call init_fe_prob
      call modify_h_bearing(fe_f, lx)
      call solve_fe_prob
   return
   endsubroutine test_bearing_x_fe


   !=========================================================================================
   !< @note Subroutine to run the bearing test
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine test_bearing_y_fe
   implicit none
      call init_fe_prob
      call modify_h_bearing_y(fe_f, ly)
      call solve_fe_prob
   return
   endsubroutine test_bearing_y_fe


   !=========================================================================================
   !< @note Subroutine to run the deterministic rough surface test
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine test_rough_fe
   implicit none
      call init_rough_prob
      call init_fe_prob
      call apply_roughness(fe_f, tab_s)
      call solve_fe_prob
   return
   endsubroutine test_rough_fe


   !=========================================================================================
   !< @note Subroutine to run the pocket test
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine test_pocket_fe
   implicit none
      call init_fe_prob
      call modify_h_pocket(fe_f, lx)
      call solve_fe_prob
   return
   endsubroutine test_pocket_fe


   !=========================================================================================
   !< @note Subroutine to run the multiscale slider test
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine test_slider_ms
   implicit none
      call init_ms_prob
      call modify_h_slider_MS(ms_fe_f, lx)
      call solve_ms_prob
   return
   endsubroutine test_slider_ms


   !=========================================================================================
   !< @note Subroutine to run the rough surface multiscale problem
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine test_rough_ms
   implicit none
      call init_rough_prob
      call init_ms_prob
      call apply_roughness_MS(ms_fe_f, tab_s)
      call solve_ms_prob
   return
   endsubroutine test_rough_ms


   !=========================================================================================
   !< @note Subroutine to initialize a deterministic 'smooth' problem
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine init_fe_prob
   implicit none
      if (mod(nx, 2) /= 1) then
         write(OPU,*) 'It is better to use an odd number of nodes in the x dir to have profile for comparison with analytical sol.'
      endif
      fe_f%m%nx = nx
      fe_f%m%ny = ny
      fe_f%m%lx = lx
      fe_f%m%ly = ly
      fe_f%m%zx = 0.
      fe_f%m%zy = 0.
      call create_rect_FE_film( data_f = data_f,   &
                                 num_p = num_pts,  &
                                  fe_f = fe_f)
      write(OPU,*) 'film created'

   return
   endsubroutine init_fe_prob


   !=========================================================================================
   !< @note Subroutine to initialize a deterministic 'rough' problem
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine init_rough_prob
   implicit none
      allocate(tab_s(nx, ny))
      tab_s = 0._R8
      call read_surf(trim(surface_file), sq, tab_s, scal_tmp)
   return
   endsubroutine init_rough_prob


   !=========================================================================================
   !< @note Subroutine to initialize a multiscale 'smooth' problem
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine init_ms_prob
   implicit none
      n_mic = (nx - 1) / n_mac +1
      ms_fe_f%ts_fe_f%m%nx = n_mac +1
      ms_fe_f%ts_fe_f%m%ny = n_mac +1
      ms_fe_f%ts_fe_f%m%lx = lx
      ms_fe_f%ts_fe_f%m%ly = ly
      ms_fe_f%ts_fe_f%m%zx = 0.
      ms_fe_f%ts_fe_f%m%zy = 0.

      call multi_scale_create_rect_fe_film(data_f = data_f,    &
                                            bs_nx = n_mic,     &
                                            bs_ny = n_mic,     &
                                          num_pts = num_pts,   &
                                          num_pbs = num_pbs,   &
                                          ms_fe_f = ms_fe_f)
      write(OPU,*) 'film created'
   return
   endsubroutine init_ms_prob


   !=========================================================================================
   !< @note Subroutine to solve a deterministic 'smooth' problem
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine solve_fe_prob
   implicit none

      call system_clock(count=cinit)
      call cpu_time(t1)
         call solve_FE_film(fe_f, mat, bc)
         call compute_corner_fluxes(fe_f, mat, bf)
      call cpu_time(t2)
      call system_clock(count=cend,count_rate=cr)

      write(OPU,*)          'FE cpu time (s):',          char(9), t2 - t1
      write(unit_num_res,*) 'FE_cpu_time_(s):',          char(9), t2 - t1
      write(OPU,*)          'FE real comp time (s):',    char(9), real(cend - cinit) / real(cr)
      write(unit_num_res,*) 'FE_real_comp_time_(s):',    char(9), real(cend - cinit) / real(cr)

      save_PeK = .false.
      save_PeE = .false.

      if (s_vtk /= NO_VTK) call save_fe_f_vtk(fe_f, trim(ms_vtk))

      select case(test_num)
         case (1)
            call save_profile_x_comp_slider(fe_f, trim(prof_ts), lx, ly / 2)
            call execute_command_line("python3 bin/pyt/filetoplot.py "// trim(prof_ts) //" '$x/L$' '$ph^2/6 \mu V L$' ")

         case (5)
            call save_profile_y_comp_air_pocket(fe_f, trim(prof_ts), ly, lx/2, bc)
            call execute_command_line("python3 bin/pyt/filetoplot.py "// trim(prof_ts) //" '$x/L$' '$ph^2/6 \mu V L$' ")

            call save_profile_x_comp_air_pocket(fe_f, trim(prof_ts)//'2', lx, ly / 2, bc)
            call execute_command_line("python3 bin/pyt/filetoplot.py "// trim(prof_ts)//'2' //" '$x/L$' '$ph^2/6 \mu V L$' ")

         case default
            call save_profile_x_fe(fe_f, trim(prof_ts), lx, ly / 2)

      endselect

      call save_fe_field(fe_f = fe_f,                                   & !
                    file_name = "out/"//res_dir//"/"//"pressure.sur",   & !
                         code = P_N,                                    & !
                        nodal = .true.)

      if (save_PeK) & !
      call save_fe_field(fe_f = fe_f,                                   & !
                    file_name = "out/"//res_dir//"/"//"PeK.sur",        & !
                         code = PEK_C,                                  & !
                        nodal = .false.)

      if (save_PeE) & !
      call save_fe_field(fe_f = fe_f,                                   & !
                    file_name = "out/"//res_dir//"/"//"PeE.sur",        & !
                         code = PEE_C,                                  & !
                        nodal = .false.)

      write(OPU,*)          'load FE (N):',  char(9), fe_f%fz()
      write(unit_num_res,*) 'load_FE_(N):',  char(9), fe_f%fz()
      write(OPU,*)          'fric/x FE (N)', char(9), fe_f%fx()
      write(unit_num_res,*) 'fric/x_FE_(N)', char(9), fe_f%fx()
      write(OPU,*)          'fric/y FE (N)', char(9), fe_f%fy()
      write(unit_num_res,*) 'fric/y_FE_(N)', char(9), fe_f%fy()
   return
   endsubroutine solve_fe_prob


   !=========================================================================================
   !< @note Subroutine to solve a 'smooth' multiscale problem
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine solve_ms_prob
   implicit none
      call system_clock(count=cinit)
      call cpu_time(t1)
         call multi_scale_solve_fe_film(ms_fe_f, ms_mat, bc)
      call cpu_time(t2)
      call system_clock(count=cend,count_rate=cr)

      write(OPU,*)          'MS cpu time (s):',          char(9), t2 - t1
      write(unit_num_res,*) 'MS_cpu_time_(s):',          char(9), t2 - t1
      write(OPU,*)          'MS real comp time (s):',    char(9), real(cend - cinit) / real(cr)
      write(unit_num_res,*) 'MS_real_comp_time_(s):',    char(9), real(cend - cinit) / real(cr)

      write(OPU,*)          'load MS_FE (N):',           char(9), ms_fe_f%ms_fz()
      write(unit_num_res,*) 'load_MS_FE_(N):',           char(9), ms_fe_f%ms_fz()
      write(OPU,*)          'fric/x MS FE (N):',         char(9), ms_fe_f%ms_fx()
      write(unit_num_res,*) 'fric/x_MS_FE_(N):',         char(9), ms_fe_f%ms_fx()
      write(OPU,*)          'fric/y MS FE (N):',         char(9), ms_fe_f%ms_fy()
      write(unit_num_res,*) 'fric/y_MS_FE_(N):',         char(9), ms_fe_f%ms_fy()

      call save_ms_fe_f_vtk(ms_fe_f, trim(ms_vtk))

      call save_profile_x_fe(ms_fe_f%ts_fe_f, trim(prof_ts), lx, ly/2)
      call save_profile_x_ms(ms_fe_f,         trim(prof_bs), lx, ly/2)

      call save_ms_field(ms_fe_f = ms_fe_f,                                   & !
                       file_name = "out/"//res_dir//"/"//"ms_pressure.sur",   & !
                            code = P_N,                                       & !
                           nodal = .true.                                     )
   return
   endsubroutine solve_ms_prob

endmodule test_musst

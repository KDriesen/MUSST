!! author: Noël Brunetiere, Arthur Francisco
!! version: 1.0.0
!! date: February, 15 2017
!! summary: Creation and resolution of multi-scale FE film

!< <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<       **MultiScale FE solution of the Reynolds equation**
!< </span>

!<# Description of the film module
!< This module can be used for a two-scale solution of the lubrication problems (more particularly for rough surface problems)
!
   !<## Definition of MS_FE_film
      !< [[MS_FE_FILM]] is a data structure containing a [[FE_FILM]] which is the Top Scale or macro-scale of the problem and a table of [[FE_FILM]] which is the Bottom Scale or micro-scale
   !<## Solution procedure
      !< This module can be used to create a [[MS_FE_FILM]], assemble the [[MS_FE_FILM]], and solve it. Some post-treatements like fluxes, forces are available.

module ms_film
use data_arch, only : I4, R8, NB_THREADS_MAX
use film
use mesh, only : MAX_NNC, MAX_NBS
use solver
use omp_lib
use bspline
use surfile
use data_film_hd
use num_param
use fluid_law
implicit none

private

! multi scale fe fluid film type
type MS_FE_FILM
!! <span style="color:green">
!!   *MS_FE_FILM* is the top-scale [[FE_FILM]] plus all of the bottom-scale [[FE_FILM]]
!! </span>
   type(FE_FILM) :: ts_fe_f                              !! *top-scale fe_film*
   type(FE_FILM), dimension(:), allocatable :: bs_fe_f   !! *bottom-scale fe_film*

   contains
      procedure :: ms_fx !! *force computation along \(x\)*
      procedure :: ms_fy !! *force computation along \(y\)*
      procedure :: ms_fz !! *force computation along \(z\)*
endtype ms_fe_film

type(SCALE_SURF) :: scal_tmp !! *object [[SCALE_SURF]]*

public :: MS_FE_FILM, multi_scale_create_rect_fe_film, multi_scale_solve_fe_film, save_ms_field

contains

   !=========================================================================================
   !< @note Subroutine to create a [[MS_FE_FILM]]
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine multi_scale_create_rect_fe_film(data_f, bs_nx, bs_ny, num_pts, num_pbs, ms_fe_f)
   implicit none
   type(DATA_FILM),  intent(inout) :: data_f     !! *data of the film*
   integer(kind=I4), intent(in   ) :: bs_nx      !! *number of nodes in \(x\) direction for the bottom scale*
   integer(kind=I4), intent(in   ) :: bs_ny      !! *number of nodes in \(y\) direction for the bottom scale*
   type(NUM_PAR),    intent(in   ) :: num_pts    !! *numerical parameters for iterative solution, top scale*
   type(NUM_PAR),    intent(in   ) :: num_pbs    !! *numerical parameters for iterative solution, bottom scale*
   type(MS_FE_FILM), intent(inout) :: ms_fe_f    !! *MS FE film*

      integer(kind=I4) :: e, ne, i1, i3

      ! creation of the top scale fe film
      call create_rect_FE_film(data_f = data_f,          &
                                num_p = num_pts,         &
                                 fe_f = ms_fe_f%ts_fe_f)

      ! allocation of the bottom scale table
      ne = ms_fe_f%ts_fe_f%m%ne
      allocate(ms_fe_f%bs_fe_f(ne))

      ! creation of the bottom scale fe_film
      do e = 1, ne
         i1    = ms_fe_f%ts_fe_f%m%con(e, 1)
         i3    = ms_fe_f%ts_fe_f%m%con(e, 3)

         ms_fe_f%bs_fe_f(e)%m%nx = bs_nx
         ms_fe_f%bs_fe_f(e)%m%ny = bs_ny
         ms_fe_f%bs_fe_f(e)%m%zx = ms_fe_f%ts_fe_f%m%x(i1)
         ms_fe_f%bs_fe_f(e)%m%zy = ms_fe_f%ts_fe_f%m%y(i1)
         ms_fe_f%bs_fe_f(e)%m%lx = ms_fe_f%ts_fe_f%m%x(i3) - ms_fe_f%bs_fe_f(e)%m%zx
         ms_fe_f%bs_fe_f(e)%m%ly = ms_fe_f%ts_fe_f%m%y(i3) - ms_fe_f%bs_fe_f(e)%m%zy
         call create_rect_FE_film(data_f = data_f,                &
                                   num_p = num_pbs,               &
                                    fe_f = ms_fe_f%bs_fe_f(e))
      enddo
   return
   endsubroutine multi_scale_create_rect_fe_film


   !=========================================================================================
   !< @note Subroutine to solve a [[MS_FE_FILM]]
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine multi_scale_solve_fe_film(ms_fe_f, ms_mat, bc)
   implicit none
   type(MS_FE_FILM),  intent(inout)                     :: ms_fe_f   !! *multi-scale FE film*
   type(MS_MAT_SOLV), intent(inout)                     :: ms_mat    !! *multi-scale solver matrices*
   real(kind=R8),     intent(in   ), dimension(MAX_NNC) :: bc        !! *boundary conditions at the corners*

      logical(kind=I4) :: decomp
      integer(kind=I4) :: i, e, it, fluid
      logical(kind=I4) :: conv, ts_init
      real(kind=R8)    :: error, relax
      integer(kind=I4) :: ass_c
      integer(kind=I4), dimension(2) :: compt

      relax = ms_fe_f%ts_fe_f%num_p%relax
      fluid = ms_fe_f%ts_fe_f%data_f%fl%fluid_type

      VERBOSE = nint(VERBOSE/10.)

      ! open mp instructions (parallele computation)
      if (VERBOSE >= 1) write(OPU,*) 'nb_threads_max', NB_THREADS_MAX
      if (NB_THREADS_MAX <= 0) NB_THREADS_MAX = 1

      !$ call omp_set_num_threads(NB_THREADS_MAX)
      if (VERBOSE >= 1) write(OPU,*) 'nb_threads_used', NB_THREADS_MAX


      ! check of matrices allocation
      allocate(ms_mat%bs_mat(ms_fe_f%ts_fe_f%m%ne))
      ms_mat%bs_mat(:)%slv_t = SOLVER_BS
      ms_mat%bs_mat(:)%first = .true.

      ! solution parameters
      ass_c  = ASS
      decomp = (ass_c == ASS)

      ! update fluid properties
      do i = 1, ms_fe_f%ts_fe_f%m%n
         ms_fe_f%ts_fe_f%vn(i,RHO_N) = ms_fe_f%ts_fe_f%data_f%fl%rho( ms_fe_f%ts_fe_f%vn(i, P_N), &
                                                                      ms_fe_f%ts_fe_f%vn(i, T_N) )

         ms_fe_f%ts_fe_f%vn(i,DRHODP_N) = ms_fe_f%ts_fe_f%data_f%fl%drhodp( ms_fe_f%ts_fe_f%vn(i, P_N), &
                                                                            ms_fe_f%ts_fe_f%vn(i, T_N) )
      enddo
      if (VERBOSE >= 2) write(OPU,*) 'fluid properties updated'

      ms_mat%ts_mat%slv_t = SOLVER_TS
      ms_mat%ts_mat%first = .true.
      ms_mat%ts_mat%nn = ms_fe_f%ts_fe_f%m%n
      ms_mat%ts_mat%ne = ms_fe_f%ts_fe_f%m%ne

      ! matrices allocation
      compt(:) = 0
      do e = 1, ms_fe_f%ts_fe_f%m%ne
         compt(1) = compt(1) +ms_fe_f%ts_fe_f%m%el_n(e)       ! for each element, the number of lines is added
         do i = 1, ms_fe_f%ts_fe_f%m%el_t(e)
            compt(2) = compt(2) + ms_fe_f%ts_fe_f%m%el_t(e)   ! " ", for each node, the number of contributions is added
         enddo
      enddo
      ms_mat%ts_mat%nvar = compt(1)
      ms_mat%ts_mat%nt   = compt(2)
      allocate( ms_mat%ts_mat%eltvar(ms_mat%ts_mat%nvar  ) )
      allocate( ms_mat%ts_mat%a_elt( ms_mat%ts_mat%nt    ) )
      allocate( ms_mat%ts_mat%eltptr(ms_mat%ts_mat%ne +1 ) )

      ! check of precomputed tables allocation
      if (.not.allocated(ms_fe_f%ts_fe_f%prc%vcal)) call init_prc_tab(ms_fe_f%ts_fe_f)

      ! convergence is false at zero iteration
      conv = .false.
      it   = 0

      ts_init = .false.
      if (ts_init) then
         call solve_FE_film(fe_f = ms_fe_f%ts_fe_f,   &
                             mat = ms_mat%ts_mat,     &
                              bc = bc,                &
                        flag_ass = .true.)
      else
         ! apply boundary conditions
         call apply_bc_FE_film_simple(fe_f = ms_fe_f%ts_fe_f,   &
                                        bc = bc)
         if (VERBOSE >= 2) write(OPU,*) 'bc applied'
         call solve_syst(mat = ms_mat%ts_mat, &
                        step = 'ini')
         if (VERBOSE >= 2) write(OPU,*) 'TS solver initialized'
      endif

      if ( sum(ms_fe_f%bs_fe_f(1)%m%ed(:)%n) > MAX_NBS ) stop 'MAX_NBS under estimated'

      ! solution loop
      do
         if (VERBOSE >= 3) write(OPU,*) "loop MS *******************", it

         if (conv) then !ù
            do e = 1, ms_fe_f%ts_fe_f%m%ne
               call solve_syst(mat = ms_mat%bs_mat(e), & !
                              step = 'end')
               if (VERBOSE >= 3) write(OPU,*) '   Matrix BS released, thread ', omp_get_thread_num()!ù
            enddo
            exit
         endif

         if (SMOOTH_MS) call smooth_ms_fe_f(ms_fe_f, code=P_N, nodal=.true.)

         if (BC_SPLINE) call interp_ts_bs_splin(ms_fe_f)

         ! assembly of the system
         call multi_scale_assembly_fe_film_reynolds(ms_fe_f = ms_fe_f, &
                                                    ms_mat  = ms_mat,  &
                                                      ass_c = ass_c)
         if (VERBOSE >= 3) write(OPU,*) 'ms reynolds assembled'
         if (VERBOSE >= 3) write(OPU,*) 'ass_c', ass_c, 'first', ms_mat%ts_mat%first

!~          if (ass_c == ASS) then
            ! some stuff can be saved here, provided the reloading of jptr, irow, ... (instead of convert_matrice_format)
            call convert_matrice_format(mat = ms_mat%ts_mat)
            if (VERBOSE >= 3) write(OPU,*) 'Matrix TS formated, thread ', omp_get_thread_num()
!~          endif

         if (ms_mat%ts_mat%first) then
            call solve_syst(mat = ms_mat%ts_mat, &
                           step = 'ana')
            ms_mat%ts_mat%first = .false. !ù
            if (VERBOSE >= 3) write(OPU,*) 'Matrix TS analyzed, thread ', omp_get_thread_num()
         endif

         ! solution of the system
!~          if (ass_c == ASS) then
            call solve_syst(mat = ms_mat%ts_mat, &
                           step = 'fac')
            if (VERBOSE >= 3) write(OPU,*) 'Matrix TS factorized, thread ', omp_get_thread_num()
!~          endif
         call solve_syst(mat = ms_mat%ts_mat, &
                        step = 'sol')
         if (VERBOSE >= 3) write(OPU,*) 'System TS solved, thread ', omp_get_thread_num()

         !if (ass_c == ASS) then
            call solve_syst(mat = ms_mat%ts_mat, &
                           step = 'fre')
            if (VERBOSE >= 3) write(OPU,*) 'Matrix factors freed, thread ', omp_get_thread_num()
         !endif

         ! error computation
         error = (sum(ms_mat%ts_mat%x ** 2) / sum(ms_fe_f%ts_fe_f%vn(:, P_N) ** 2)) ** (0.5)
         it = it + 1
         if (VERBOSE >= 1) write(OPU,*) 'Iteration ', it, 'error TS', error

         ! convergence check
         if (error <= ms_fe_f%ts_fe_f%num_p%eps) conv = .true.

         ! update of variables
         if (fluid == MIXT) then
            do i = 1, ms_fe_f%ts_fe_f%m%n
               if (ms_mat%ts_mat%x(i) < 0.) then
                  ms_fe_f%ts_fe_f%vn(i, RHO_N) = ms_fe_f%ts_fe_f%vn(i, RHO_N) + ms_fe_f%ts_fe_f%vn(i, DRHODP_N) * ms_mat%ts_mat%x(i) * relax
                  if (ms_fe_f%ts_fe_f%vn(i, RHO_N) < 0.) then
                      ms_fe_f%ts_fe_f%vn(i,   P_N) = ms_fe_f%ts_fe_f%data_f%fl%p_0/100
                  else
                      ms_fe_f%ts_fe_f%vn(i,   P_N) = ms_fe_f%ts_fe_f%data_f%fl%pres( ms_fe_f%ts_fe_f%vn(i, RHO_N), & !
                                                                                     ms_fe_f%ts_fe_f%vn(i,   T_N) )
                  endif
               else
                  ms_fe_f%ts_fe_f%vn(i, P_N) = ms_fe_f%ts_fe_f%vn(i, P_N) + ms_mat%ts_mat%x(i) * relax
               endif
            enddo
         else
            ms_fe_f%ts_fe_f%vn(:, P_N) = ms_fe_f%ts_fe_f%vn(:, P_N) + ms_mat%ts_mat%x * relax
         endif

         ! check pressure
         if ( fluid == GP ) then
            if (minval(ms_fe_f%ts_fe_f%vn(:, P_N)) < 0._R8) write(OPU,*) 'MS p negative'
            where (ms_fe_f%ts_fe_f%vn(:, P_N) < 0._R8) ms_fe_f%ts_fe_f%vn(:, P_N) = ms_fe_f%ts_fe_f%data_f%fl%p_0 / 1.e2_R8
         endif

         ! update fluid properties
         do i = 1, ms_fe_f%ts_fe_f%m%n
            ms_fe_f%ts_fe_f%vn(i, RHO_N) = ms_fe_f%ts_fe_f%data_f%fl%rho( ms_fe_f%ts_fe_f%vn(i, P_N), &
                                                                          ms_fe_f%ts_fe_f%vn(i, T_N) )

            ms_fe_f%ts_fe_f%vn(i, DRHODP_N) = ms_fe_f%ts_fe_f%data_f%fl%drhodp( ms_fe_f%ts_fe_f%vn(i, P_N), &
                                                                                ms_fe_f%ts_fe_f%vn(i, T_N) )
         enddo

         if (it >= ms_fe_f%ts_fe_f%num_p%it_max) then
            conv = .true.
            write(OPU,*) 'maximum number of iteration reached before convergence'
         endif


      enddo

      call solve_syst(mat = ms_mat%ts_mat, & !
                     step = 'end')
      if (VERBOSE >= 2) write(OPU,*) 'Matrix TS released, thread ', omp_get_thread_num()

   return
   endsubroutine multi_scale_solve_fe_film


   !=========================================================================================
   !< @note Subroutine to interpolate the top-scale nodes for the bottom-scale boundaries
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine interp_ts_bs_splin(ms_fe_f)
   implicit none
   type(MS_FE_FILM), intent(inout) :: ms_fe_f

      integer(kind=I4) :: i, j, inbvx, inbvy, iloy, iflag, ts_nx, ts_ny, e, ed, bs_nx, bs_ny, n_edg, n_nod, nn, node_number
      real(kind=R8) :: x_n, y_n, val
      real(kind=R8), dimension(:,:), allocatable :: log_press

      real(kind=R8), dimension(:,:), allocatable :: coeff
      real(kind=R8), dimension(:),   allocatable :: tx ! x knots
      real(kind=R8), dimension(:),   allocatable :: ty ! y knots
      real(kind=R8), dimension(:),   allocatable ::  x
      real(kind=R8), dimension(:),   allocatable ::  y

      integer(kind=I4), parameter :: deg = 2
      logical(kind=I4) :: mixture

      mixture = (ms_fe_f%ts_fe_f%data_f%fl%fluid_type == MIXT)

      ts_nx = ms_fe_f%ts_fe_f%m%nx
      ts_ny = ms_fe_f%ts_fe_f%m%ny

      allocate( coeff(1:ts_nx, 1:ts_ny) )
      allocate( tx(1:(ts_nx +deg)),  &
                ty(1:(ts_ny +deg)) )
      allocate(  x(1: ts_nx    ),  &
                 y(1: ts_ny    ) )

      allocate( log_press(1:ts_nx, 1:ts_ny) )
      log_press = reshape( ms_fe_f%ts_fe_f%vn(:, P_N), (/ts_nx, ts_ny/) )

      if (mixture) log_press = log(log_press)

      x(1:ts_nx) = ms_fe_f%ts_fe_f%m%x(1:ts_nx)
      i = 0
      do j = 1, ts_nx*(ts_ny -1) +1, ts_nx
         i = i +1
         y(i) = ms_fe_f%ts_fe_f%m%y(j)
      enddo

      iflag = 0
      call db2ink(   x = x(1:ts_nx),                    & ! Array of x abcissae. Must be strictly increasing.
                    nx = ts_nx,                         & ! Number of x abcissae
                     y = y(1:ts_ny),                    & ! Array of y abcissae. Must be strictly increasing.
                    ny = ts_ny,                         & ! Number of y abcissae
                   fcn = log_press(1:ts_nx, 1:ts_ny),   & ! Array of function values to interpolate. fcn(i,j) should
                                                          !    contain the function value at the point (x(i),y(j))
                    kx = deg,                           & ! The order of spline pieces in x (>= 2, < nx). (order = polynomial degree + 1)
                    ky = deg,                           & ! The order of spline pieces in y (>= 2, < ny). (order = polynomial degree + 1)
                    tx = tx(1:(ts_nx +deg)),            & ! The knots in the x direction for the spline interpolant.
                                                          !    If iflag=0 these are chosen by [[db2ink]].
                                                          !    If iflag=1 these are specified by the user.
                                                          !    Must be non-decreasing.
                    ty = ty(1:(ts_ny +deg)),            & ! The knots in the y direction for the spline interpolant.
                                                          !    If iflag=0 these are chosen by [[db2ink]].
                                                          !    If iflag=1 these are specified by the user.
                                                          !    Must be non-decreasing.
                 bcoef = coeff(1:ts_nx, 1:ts_ny),       & ! Array of coefficients of the b-spline interpolant.
                 iflag = iflag)                           ! **on input:**  0 = knot sequence chosen by [[db2ink]].
                                                          !                1 = knot sequence chosen by user.
                                                          ! **on output:** 1 = successful execution.
                                                          !                2 = iflag out of range.
                                                          !                3 = nx out of range.
                                                          !                4 = kx out of range.
                                                          !                5 = x not strictly increasing.
                                                          !                6 = tx not non-decreasing.
                                                          !                7 = ny out of range.
                                                          !                8 = ky out of range.
                                                          !                9 = y not strictly increasing.
                                                          !               10 = ty not non-decreasing.
      if (iflag/=1) error stop 'error calling db2ink'

      inbvx = 1
      inbvy = 1
      iloy  = 1
      do e = 1, ms_fe_f%ts_fe_f%m%ne
         bs_nx = ms_fe_f%bs_fe_f(e)%m%nx
         bs_ny = ms_fe_f%bs_fe_f(e)%m%ny
         n_edg = ms_fe_f%bs_fe_f(e)%m%ned
         do ed = 1, n_edg
            n_nod = ms_fe_f%bs_fe_f(e)%m%ed(ed)%n
            do nn = 1, n_nod
               node_number = ms_fe_f%bs_fe_f(e)%m%ed(ed)%nm(nn)
               x_n = ms_fe_f%bs_fe_f(e)%m%x(node_number)
               y_n = ms_fe_f%bs_fe_f(e)%m%y(node_number)

               call db2val(xval = x_n,                      &  ! xval     !! x coordinate of evaluation point.
                           yval = y_n,                      &  ! yval     !! y coordinate of evaluation point.
                            idx = 0,                        &  ! idx      !! x derivative of piecewise polynomial to evaluate.
                            idy = 0,                        &  ! idy      !! y derivative of piecewise polynomial to evaluate.
                             tx = tx(1:(ts_nx +deg)),       &  ! tx       !! sequence of knots defining the piecewise polynomial in the x direction. (same as in last call to [[db2ink]])
                             ty = ty(1:(ts_ny +deg)),       &  ! ty       !! sequence of knots defining the piecewise polynomial in the y direction. (same as in last call to [[db2ink]])
                             nx = ts_nx,                    &  ! nx       !! the number of interpolation points in x. (same as in last call to [[db2ink]])
                             ny = ts_ny,                    &  ! ny       !! the number of interpolation points in y. (same as in last call to [[db2ink]])
                             kx = deg,                      &  ! kx       !! order of polynomial pieces in x. (same as in last call to [[db2ink]])
                             ky = deg,                      &  ! ky       !! order of polynomial pieces in y. (same as in last call to [[db2ink]])
                          bcoef = coeff(1:ts_nx, 1:ts_ny),  &  ! bcoef    !! the b-spline coefficients computed by [[db2ink]].
                              f = val,                      &  ! f        !! interpolated value &
                          iflag = iflag,                    &  ! iflag    !! status flag: 0 : no errors, /=0 : error
                          inbvx = inbvx,                    &  ! inbvx    !! initialization parameter which must be set to 1 the first time this routine is called, and must not be changed by the user.
                          inbvy = inbvy,                    &  ! inbvy    !! initialization parameter which must be set to 1 the first time this routine is called, and must not be changed by the user.
                           iloy = iloy)                        ! iloy     !! initialization parameter which must be set to 1 the first time this routine is called, and must not be changed by the user.
               if (iflag/=0) error stop 'error'
               if (mixture) val = exp(val)

               ms_fe_f%bs_fe_f(e)%vn(node_number, P_N) = val
               ms_fe_f%bs_fe_f(e)%bc(node_number, REY) = 0
            enddo
         enddo
      enddo

      deallocate( log_press, coeff, tx, ty, x, y )

   return
   endsubroutine interp_ts_bs_splin


   !=========================================================================================
   !< @note Subroutine to save the information controlled by *code* of the whole mesh.
   !
   !        The genertaed file is a ```.sur``` file.
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine save_ms_field(ms_fe_f, file_name, code, nodal)
   implicit none
   type(MS_FE_FILM), intent(in) :: ms_fe_f
   character(len=*), intent(in) :: file_name    !! *filename like "./out/pressure.sur"*
   integer(kind=I4), intent(in) :: code         !! *saved information like P_N*
   logical(kind=I4), intent(in) :: nodal        !! *if false : cell value, if true : nodal value*

      integer(kind=I4) :: nnx, nny
      character(len=8) :: unit_z, string
      integer(kind=I4), dimension(1) :: i1, i2
      real(kind=R8), allocatable, dimension(:,:) :: mat

      call ms_fe_f_2_mat(ms_fe_f, code, nodal, mat)

      nnx = ubound(mat, 1)
      nny = ubound(mat, 2)

      call empty(unit_z)

      if (nodal) then
         string = trim(ms_fe_f%ts_fe_f%vn_name(code))
      else
         string = trim(ms_fe_f%ts_fe_f%vc_name(code))
      endif

      i1 = index(string, '(') +1
      i2 = index(string, ')') -1
      unit_z = string(i1(1):i2(1))

      call init_scal(scal = scal_tmp,                   & !
                       nx = nnx,                        & !
                       ny = nny,                        & !
                       lx = ms_fe_f%ts_fe_f%m%lx,       & ! default unit : m
                       ly = ms_fe_f%ts_fe_f%m%ly,       & !
                   unit_z = unit_z                      ) !

      call write_surf(nom_fic = file_name,  & !
                        tab_s = mat,        & !
                         scal = scal_tmp    )

      deallocate(mat)

   return
   endsubroutine save_ms_field


   !=========================================================================================
   !< @note Subroutine to transform a MS_FE_FILM information into a matrix
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine ms_fe_f_2_mat(ms_fe_f, code, nodal, mat)
   implicit none
   type(MS_FE_FILM), intent(in   )              :: ms_fe_f
   integer(kind=I4), intent(in   )              :: code        !! *saved information like P_N*
   logical(kind=I4), intent(in   )              :: nodal       !! *if false : cell value, if true : nodal value*
   real(kind=R8),    intent(inout), allocatable :: mat(:,:)    !! *output matrix containing the information*

      integer(kind=I4) :: ts_nx, ts_ny, ex, ey, e, ne_x, ne_y, nnx, nny, c
      integer(kind=I4), allocatable, dimension(:) :: bs_nx, bs_ny

      ts_nx = ms_fe_f%ts_fe_f%m%nx -1
      ts_ny = ms_fe_f%ts_fe_f%m%ny -1

      allocate(bs_nx(ts_nx*ts_ny), bs_ny(ts_nx*ts_ny))

      bs_nx(:) = ms_fe_f%bs_fe_f(:)%m%nx -1
      bs_ny(:) = ms_fe_f%bs_fe_f(:)%m%ny -1

      nnx = sum( bs_nx(1:ts_nx) )
      nny = sum( bs_ny(1:ts_ny) )
      if (nodal) then
         nnx = nnx +1
         nny = nny +1
      endif

      allocate( mat(1:nnx, 1:nny) ) ; mat = -1.

      c = 0
      if (nodal) c = 1
      do ey = 1, ts_ny
      do ex = 1, ts_nx
         e    = ts_nx*(ey -1) +ex
         ne_x = bs_nx(e)*(ex -1) +1
         ne_y = bs_ny(e)*(ey -1) +1
         mat( ne_x:( ne_x +bs_nx(e) -1 +c), &
              ne_y:( ne_y +bs_ny(e) -1 +c) ) = reshape( ms_fe_f%bs_fe_f(e)%vn(:, code), (/bs_nx(e) +c, bs_ny(e) +c/) )
      enddo
      enddo

      deallocate(bs_nx, bs_ny)
   return
   endsubroutine ms_fe_f_2_mat


   !=========================================================================================
   !< @note Subroutine to transform a matrix into a MS_FE_FILM
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine mat_2_ms_fe_f(ms_fe_f, code, nodal, mat)
   implicit none
   type(MS_FE_FILM), intent(inout)                 :: ms_fe_f
   integer(kind=I4), intent(in   )                 :: code     !! *saved information like P_N*
   logical(kind=I4), intent(in   )                 :: nodal    !! *if false : cell value, if true : nodal value*
   real(kind=R8),    intent(in   ), dimension(:,:) :: mat      !! *input matrix containing the information*

      integer(kind=I4) :: ts_nx, ts_ny, ex, ey, e, ne_x, ne_y, nnx, nny, c
      integer(kind=I4), allocatable, dimension(:) :: bs_nx, bs_ny

      ts_nx = ms_fe_f%ts_fe_f%m%nx -1
      ts_ny = ms_fe_f%ts_fe_f%m%ny -1

      allocate(bs_nx(ts_nx*ts_ny), bs_ny(ts_nx*ts_ny))

      bs_nx(:) = ms_fe_f%bs_fe_f(:)%m%nx -1
      bs_ny(:) = ms_fe_f%bs_fe_f(:)%m%ny -1

      nnx = sum( bs_nx(1:ts_nx) )
      nny = sum( bs_ny(1:ts_ny) )
      if (nodal) then
         nnx = nnx +1
         nny = nny +1
      endif

      c = 0
      if (nodal) c = 1
      do ey = 1, ts_ny
      do ex = 1, ts_nx
         e    = ts_nx*(ey -1) +ex
         ne_x = bs_nx(e)*(ex -1) +1
         ne_y = bs_ny(e)*(ey -1) +1

         ms_fe_f%bs_fe_f(e)%vn(:, code) = reshape( mat( ne_x:( ne_x +bs_nx(e) -1 +c), & !
                                                        ne_y:( ne_y +bs_ny(e) -1 +c) ), (/(bs_nx(e) +c)*(bs_ny(e) +c)/) )
      enddo
      enddo

      deallocate(bs_nx, bs_ny)
   return
   endsubroutine mat_2_ms_fe_f


   !=========================================================================================
   !< @note Subroutine to smooth a MS_FE_FILM field, the pressure for instance.
   !
   !        By default, the smoothing kernel is a 5x5 Gaussian filter
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine smooth_ms_fe_f(ms_fe_f, code, nodal)
   implicit none
   type(MS_FE_FILM), intent(inout) :: ms_fe_f
   integer(kind=I4), intent(in)    :: code      !! *saved information like P_N*
   logical(kind=I4), intent(in)    :: nodal     !! *if false : cell value, if true : nodal value*

      integer(kind=I4) :: nx, ny
      real(kind=R8), allocatable, dimension(:,:) :: mat

      call ms_fe_f_2_mat(ms_fe_f, code, nodal, mat)
      nx = ubound(mat, 1) ; ny = ubound(mat, 2)

      call smooth_mat(mat, nx, ny, s=5)

      call mat_2_ms_fe_f(ms_fe_f, code, nodal, mat)

      deallocate(mat)
   return
   endsubroutine smooth_ms_fe_f


   !=========================================================================================
   !< @note Subroutine to smooth a matrix form field, the pressure for instance.
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine smooth_mat(mat, nx, ny, s)
   implicit none
   integer(kind=I4), intent(in)                       :: nx    !! *matrix x size*
   integer(kind=I4), intent(in)                       :: ny    !! *matrix y size*
   integer(kind=I4), intent(in)                       :: s     !! *kernel size*
   real(kind=R8),    intent(inout), dimension(nx, ny) :: mat   !! *matrix*

      integer(kind=I4), dimension(3,3) :: kernel3
      integer(kind=I4), dimension(5,5) :: kernel5
      real(kind=R8),    dimension(5,5) :: mi_mat
      real(kind=R8),    dimension(:,:), allocatable :: mat_tmp

      integer(kind=I4) :: i, j

      allocate( mat_tmp(nx, ny) )

      if (s==3) then
         kernel3(1:3, 1:3) = reshape((/1,2,1, & !
                                       2,4,2, & !
                                       1,2,1/), (/3,3/))
         do j = 1 +1, nx -1
         do i = 1 +1, ny -1
            mi_mat(1:3, 1:3) = mat(i-1:i+1, j-1:j+1)
            mat_tmp(i, j) = sum( mi_mat(1:3, 1:3)*kernel3(1:3, 1:3) )/16.
         enddo
         enddo

         do j = 1 +1, nx -1
         do i = 1 +1, ny -1
            mat(i, j) = mat_tmp(i, j)
         enddo
         enddo

      endif

      if (s==5) then
         kernel5(1:5, 1:5) = reshape((/ 1, 4, 6, 4, 1, & !
                                        4,16,24,16, 4, & !
                                        6,24,36,24, 6, & !
                                        4,16,24,16, 4, & !
                                        1, 4, 6, 4, 1 /), (/5,5/))
         do j = 1 +2, nx -2
         do i = 1 +2, ny -2
            mi_mat(1:5, 1:5) = mat(i-2:i+2, j-2:j+2)
            mat_tmp(i, j) = sum( mi_mat(1:5, 1:5)*kernel5(1:5, 1:5) )/256.
         enddo
         enddo

         do j = 1 +2, nx -2
         do i = 1 +2, ny -2
            mat(i, j) = mat_tmp(i, j)
         enddo
         enddo

      endif


      deallocate(mat_tmp)

   return
   endsubroutine smooth_mat


   !=========================================================================================
   !< @note Subroutine to assemble the top-scale system, all of the bottom-scale systems being solved.
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine multi_scale_assembly_fe_film_reynolds(ms_fe_f, ms_mat, ass_c)
   implicit none
   type(MS_FE_FILM),  intent(inout) :: ms_fe_f     !!
   type(MS_MAT_SOLV), intent(inout) :: ms_mat      !! *solver type matrices table*
   integer(kind=I4),  intent(in)    :: ass_c       !! *assembly type*

      real(kind=R8),    dimension(:,:,:), allocatable :: ke_ij
      real(kind=R8),    dimension(:,:)  , allocatable ::  be_i
      integer(kind=I4), dimension(:,:)  , allocatable :: ind_e
      integer(kind=I4), dimension(2)                  :: compt

      integer(kind=I4) :: e, i, ii, ne
      real(kind=R8)    :: val_t
      logical(kind=I4) :: first_assembly

      ne = ms_fe_f%ts_fe_f%m%ne

      allocate( ke_ij(ne, MAX_NNC, MAX_NNC) )
      allocate(  be_i(ne, MAX_NNC) )
      allocate( ind_e(ne, MAX_NNC) )

      ! assembly
      compt(:)        = 1
      ms_mat%ts_mat%b = 0._R8

      if (.not.allocated(ms_mat%ass_loc_in_mat)) then
         allocate(ms_mat%ass_loc_in_mat(ne))
         ms_mat%ass_loc_in_mat = -1
      endif

      ! check for first assembly: save the assembly location in ass_loc_in_mat
      if (ms_mat%ass_loc_in_mat(1) == -1) then
         first_assembly = .true.
         ms_mat%ass_loc_in_mat(1) = 1
      else
         first_assembly = .false.
      endif

      do e = 1, ne
         ! copy of the boundary conditions
         do i = 1, 4
            ii = ms_fe_f%ts_fe_f%m%con(e, i)
            be_i(e, i) = ms_fe_f%ts_fe_f%vn(ii, P_N)
         enddo
      enddo

      ! elementary matrices calculation
      !-------------------------------------------
      !open mp instructions (parallele computation)
      !$omp parallel
      !$omp do schedule(runtime)
      do e = 1, ms_fe_f%ts_fe_f%m%ne
         if (VERBOSE >= 3) write(OPU,*) "MS elements ", e, ' thread ', omp_get_thread_num()

         call elementary_full_domain_FE_film_reynolds(fe_f = ms_fe_f%bs_fe_f(e), &
                                                       mat = ms_mat%bs_mat(e),   &
                                                     ke_ij = ke_ij(e, :, :),     &
                                                      be_i = be_i( e, :),        &
                                                     ind_e = ind_e(e, :))
      enddo
      !$OMP end do
      !$OMP end parallel
      !end open mp instructions
      !-------------------------------------------

      do e = 1, ne
         ! consideration of the bc condition
         do i = 1, 4
            ii = ms_fe_f%ts_fe_f%m%con(e, i)
            ind_e(e, i) = ii
            if (ms_fe_f%ts_fe_f%bc(ii, REY) == 0) then
               val_t = ke_ij(e, i, i)
               ke_ij(e, i, :) = 0._R8
               ke_ij(e, i, i) = val_t
               be_i( e, i)    = 0._R8
            endif
         enddo
         ! copy of the rhs member
         do i = 1, 4
            ii = ms_fe_f%ts_fe_f%m%con(e, i)
            ms_mat%ts_mat%b(ii) = ms_mat%ts_mat%b(ii) + be_i(e, i)
         enddo
      enddo

      ! assembly of the elemental matrix in the solver matrix
      if (ass_c == ASS) then

         do e = 1, ne

            call assemble_in_mat_sol(mat = ms_mat%ts_mat,      &
                                     num = e,                  &
                                    nelt = 4,                  &
                                   nline = 4,                  &
                                    tind = ind_e(e, :),        &
                                   m_elt = ke_ij(e, :, :),     &
                                   compt = compt)
         enddo
         ms_mat%ts_mat%eltptr(1) = 1

      endif

      deallocate( ke_ij, be_i, ind_e )

   return
   endsubroutine multi_scale_assembly_fe_film_reynolds


   !=========================================================================================
   !< @note Function to calculate the generated load in a fluid MS film
   !        \( \int_\Omega p d\Omega \)
   !  @endnote
   !-----------------------------------------------------------------------------------------
   real(kind=R8) function ms_fz(ms_fe_f)
   implicit none
   class(MS_FE_FILM), intent(inout) :: ms_fe_f !! *MS FE film*

      integer(kind=I4) :: e
      ms_fz = 0._R8
      do e = 1, ms_fe_f%ts_fe_f%m%ne
          ms_fz = ms_fz + ms_fe_f%bs_fe_f(e)%fz()
      enddo
   return
   endfunction ms_fz


   !=========================================================================================
   !< @note Function to calculate the friction force along \(x\) in a fluid MS film
   !        \( \int_\Omega \tau_{xz} d\Omega \)
   !  @endnote
   !-----------------------------------------------------------------------------------------
   real(kind=R8) function ms_fx(ms_fe_f)
   implicit none
   class(MS_FE_FILM), intent(inout) :: ms_fe_f

      integer(kind=I4) :: e
      ms_fx = 0._R8
      do e = 1, ms_fe_f%ts_fe_f%m%ne
         ms_fx = ms_fx + ms_fe_f%bs_fe_f(e)%fx()
      enddo
   return
   endfunction ms_fx


   !=========================================================================================
   !< @note Function to calculate the friction force along \(y\) in a fluid MS film
   !        \( \int_\Omega \tau_{yz} d\Omega \)
   !  @endnote
   !-----------------------------------------------------------------------------------------
   real(kind=R8) function ms_fy(ms_fe_f)
   implicit none
   class(MS_FE_FILM), intent(inout) :: ms_fe_f

      integer(kind=I4) :: e
      ms_fy = 0._R8
      do e = 1, ms_fe_f%ts_fe_f%m%ne
         ms_fy = ms_fy + ms_fe_f%bs_fe_f(e)%fy()
      enddo
   return
   endfunction ms_fy

endmodule ms_film

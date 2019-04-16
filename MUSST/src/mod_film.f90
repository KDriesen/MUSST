
!! author: Noël Brunetière<br/>&emsp;Arthur Francisco
!! version: 1.0.0
!! date: April,17 2017
!! summary: Definition of a fluid film by FEM

!< <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!        **FE solution of the Reynolds equation**
!  </span>

!<# Description of the film module
!  This is the main module for the finite element solution of the Reynolds equation.
!
!<## Definition of FE_film
! [[FE_FILM]] is a data structure containing a [[FE_MESH]] and some data describing a lubrication problem
!<## Solution procedure
! This module can be used to create a FE_FILM, assemble the FE_FIM and solve it. Some post-treatements like fluxes, forces are available.

module film
use data_arch, only : EPS_R8
use data_film_hd
use fluid_law
use num_param
use mesh
use solver
use elements
use omp_lib
use surfile

implicit none

private

type PRC_TAB
!! <span style="color:green">
!!   *PRC_TAB* stores some precomputed coefficients for the finite element matrices
!! </span>
   integer(kind=I4)                             :: ng          !! *Gauss point number along a direction*
   real(kind=R8), dimension(:),     allocatable :: pg          !! *point coordinates in a direction*
   real(kind=R8), dimension(:),     allocatable :: wg          !! *point weight*
   real(kind=R8), dimension(:,:,:), allocatable :: vni4        !! *for each node, shape function at Gauss points*
   real(kind=R8), dimension(:,:,:), allocatable :: vni4x       !! *for each node, shape function derivative at Gauss points*
   real(kind=R8), dimension(:,:,:), allocatable :: vni4y       !! *for each node, shape function derivative at Gauss points*
   real(kind=R8), dimension(:,:,:), allocatable :: vni4d       !! *for each node, upwind shape function at Gauss points*
   real(kind=R8), dimension(:,:,:), allocatable :: vcal        !! *14 values calculated at the Gauss points*
endtype PRC_TAB

type FE_FILM
!! <span style="color:green">
!!   *FE_FILM* stores the whole stuff related to a film: the nodal variables, the mesh, etc.
!! </span>
   integer(kind=I4)  :: n_vn          !! *number of nodal variables*
   integer(kind=I4)  :: n_vc          !! *number of variables on cells*
   type(FE_MESH)     :: m             !! *mesh of the film*
   type(DATA_FILM)   :: data_f        !! *data of the problem*
   type(PRC_TAB)     :: prc           !! *precomputation*
   type(NUM_PAR)     :: num_p         !! *numerical param for iterative solution*

   real(kind=R8),     dimension(:,:), allocatable :: vn      !! *nodal variables table*
   real(kind=R8),     dimension(:,:), allocatable :: vc      !! *cell variables table*
   integer(kind=I4),  dimension(:,:), allocatable :: bc      !! *boundary nodes code*
   character(len=20), dimension(:),   allocatable :: vn_name !! *nodal variable names*
   character(len=20), dimension(:),   allocatable :: vc_name !! *cell variable names*
   contains
      procedure :: fx                                        !! *force computation along* \(x\)
      procedure :: fy                                        !! *force computation along* \(y\)
      procedure :: fz                                        !! *force computation along* \(z\)
endtype FE_FILM

! codes for variables
integer(kind=I4), parameter :: H1_N     =  1  !! *code for bottom surface height* \(h_1\)
integer(kind=I4), parameter :: H2_N     =  2  !! *code for top surface height* \(h_2\)
integer(kind=I4), parameter :: H_N      =  3  !! *code for film thickness* \(h\)
integer(kind=I4), parameter :: P_N      =  4  !! *code for film absolute pressure* \(p\)
integer(kind=I4), parameter :: RHO_N    =  5  !! *code for fluid density* \(\rho\)
integer(kind=I4), parameter :: T_N      =  6  !! *code for film absolute temperature* \(T\)
integer(kind=I4), parameter :: DRHODP_N =  7  !! *code for film compressibility* \(\frac{\partial \rho}{\partial p}\)
integer(kind=I4), parameter :: MU_N     =  8  !! *code for fluid viscosity* \(\mu\)
integer(kind=I4), parameter :: VX_N     =  9  !! *code for surface velocity along* \(x\) \(V_x\). *Should be modified (surfaces 1 and 2)*
integer(kind=I4), parameter :: VY_N     = 10  !! *code for surface velocity along* \(y\) \(V_y\). *Should be modified (surfaces 1 and 2)*
integer(kind=I4), parameter :: HG_C     =  1  !! *code for groove depth on the stationnary surface* \(h_g\), *cell variable*
integer(kind=I4), parameter :: PEK_C    =  2  !! *code for Peclet number along* \(\xi\)
integer(kind=I4), parameter :: PEE_C    =  3  !! *code for Peclet number along* \(\eta\)

! codes for boundary conditions
integer(kind=I4), parameter :: REY      =  1  !! *code for imposed pressure at the boundary*

! codes for assembly
integer(kind=I4), parameter :: ASS      =  1  !! *code for assembly of the system*
integer(kind=I4), parameter :: NO_ASS   =  0  !! *code for no assemble, computation of the residual only*
integer(kind=I4), parameter :: NO_BC    = -1  !! *code for computation of the residual everywhere, even at the boundaries (fluxes computations)*

! private variable for a '.sur' file creation
type(SCALE_SURF) :: scal_tmp !! *object [[SCALE_SURF]]*

! parameters controlling the iterative process
logical(kind=I4), parameter :: BC_SPLINE = .false. !! *instead of linearly interpolating the boundary pressures, it can be done in a smoother way. NOT TO BE USED YET.*

public :: ASS, DRHODP_N, FE_FILM, H1_N, H2_N, HG_C, H_N, MU_N, NO_ASS, NO_BC, PEE_C, PEK_C, PRC_TAB, &
          P_N, REY, RHO_N, T_N, VX_N, VY_N, apply_bc_FE_film, apply_bc_FE_film_simple, assemble_in_mat_sol, assembly_FE_film_reynolds, compute_corner_fluxes, &
          elementary_full_domain_FE_film_reynolds, init_prc_tab, solve_FE_film, create_rect_FE_film, save_fe_field, BC_SPLINE

contains

   !=========================================================================================
   !< @note Subroutine to create a [[FE_FILM]]
   !<
   !-----------------------------------------------------------------------------------------
   subroutine create_rect_FE_film(data_f, num_p, fe_f)
   implicit none
   type(DATA_FILM), intent(inout) :: data_f !! *data of the film*
   type(NUM_PAR),   intent(in   ) :: num_p  !! *numerical param for iterative solution*
   type(FE_FILM),   intent(inout) :: fe_f   !! *FE film data type*

      ! mesh creation
      call create_rect_x_ymesh(fe_f%m)
      ! copy of data
      fe_f%data_f = data_f
      ! allocation and initialisation of the variables and bc table
      select case (fe_f%data_f%pb_type)
         case(HD)
            ! nodal tables
            fe_f%n_vn = 10
            allocate (fe_f%vn(fe_f%m%n, fe_f%n_vn), fe_f%vn_name(fe_f%n_vn))
            fe_f%vn_name(H1_N)  = 'h1(m)'
            fe_f%vn_name(H2_N)  = 'h2(m)'
            fe_f%vn_name(H_N)   = 'h(m)'
            fe_f%vn_name(P_N)   = 'p(Pa)'
            fe_f%vn_name(RHO_N) = 'rho(kg.m-3)'
            fe_f%vn_name(T_N)   = 'T(K)'
            fe_f%vn_name(DRHODP_N) = 'drhodp(kg.m-3.Pa-1)'
            fe_f%vn_name(MU_N)     = 'mu(Pa.s)'
            fe_f%vn_name(VX_N)     = 'Vx(m.s-1)'
            fe_f%vn_name(VY_N)     = 'Vy(m.s-1)'

            fe_f%vn           = 0._R8
            fe_f%vn(:, H2_N)  = fe_f%data_f%h_0
            fe_f%vn(:, H_N)   = fe_f%vn(:, H2_N) - fe_f%vn(:, H1_N)
            fe_f%vn(:, P_N)   = fe_f%data_f%fl%p_0
            fe_f%vn(:, RHO_N) = fe_f%data_f%fl%rho(fe_f%data_f%fl%p_0, fe_f%data_f%fl%T_0)
            fe_f%vn(:, T_N)   = fe_f%data_f%fl%T_0
            fe_f%vn(:, MU_N)  = fe_f%data_f%fl%mu_0
            fe_f%vn(:, VX_N)  = fe_f%data_f%v_x
            fe_f%vn(:, VY_N)  = fe_f%data_f%v_y
            fe_f%vn(:, DRHODP_N) = fe_f%data_f%fl%drhodp(fe_f%data_f%fl%p_0, fe_f%data_f%fl%T_0)
            ! bondary condition table
            allocate (fe_f%bc(fe_f%m%n, 1))
            ! all the nodes are initialized as unknown
            fe_f%bc = 1
            ! cell variable
            fe_f%n_vc = 3
            allocate (fe_f%vc(fe_f%m%ne, fe_f%n_vc), fe_f%vc_name(fe_f%n_vc))
            fe_f%vc_name(HG_C) = 'h_grv(m)'
            fe_f%vc_name(PEK_C) = 'Pek'
            fe_f%vc_name(PEE_C) = 'Pee'
            fe_f%vc = 0._R8
            ! numerical parameters for iterative problems
            fe_f%num_p = num_p
         case default
            stop 'the problem type is undefined in create_rect_FE_film'
      endselect
   return
   endsubroutine create_rect_FE_film

   !=========================================================================================
   !< @note Subroutine to solve a [[FE_FILM]]
   !<
   !-----------------------------------------------------------------------------------------
   subroutine solve_FE_film(fe_f, mat, bc, flag_ass)
   implicit none
   type(FE_FILM),    intent(inout)                     :: fe_f     !! *FE film data type*
   type(MAT_SOLV),   intent(inout)                     :: mat      !! *matrices for solving*
   real(kind=R8),    intent(in   ), dimension(MAX_NNC) :: bc       !! *boundary conditions at the corners*
   logical(kind=I4), intent(in   ), optional           :: flag_ass !! *optional flag for assembly*

      logical(kind=I4) :: decomp
      integer(kind=I4) :: ass_c, i, it, e
      logical(kind=I4) :: conv
      real(kind=R8)    :: error
      integer(kind=I4), dimension(2) :: compt

      ! solution parameters
      ass_c = ASS
      if (present(flag_ass)) then
         if (.not.flag_ass) ass_c = NO_ASS
      endif
      decomp = (ass_c == ASS)

      ! apply boundary conditions
      if (BC_SPLINE) then
         call apply_bc_FE_film(fe_f = fe_f,    &
                                 bc = bc)
      else
         call apply_bc_FE_film_simple(fe_f = fe_f,    &
                                        bc = bc)
      endif

      if (VERBOSE >= 20) write(OPU,*) 'bc applied'

      ! update fluid properties
      do i = 1, fe_f%m%n
         fe_f%vn(i, RHO_N)    = fe_f%data_f%fl%rho(   fe_f%vn(i, P_N) , fe_f%vn(i, T_N))
         fe_f%vn(i, DRHODP_N) = fe_f%data_f%fl%drhodp(fe_f%vn(i, P_N) , fe_f%vn(i, T_N))
      enddo
      if (VERBOSE >= 20) write(OPU,*) 'fluid properties updated'

      if (mat%first) then
         mat%nn = fe_f%m%n
         mat%ne = fe_f%m%ne
         call solve_syst(mat = mat, &
                        step = 'ini')
         if (VERBOSE >= 10) write(OPU,*) 'Matrix initialized, thread ', omp_get_thread_num()

         ! matrices allocation
         compt(:) = 0
         do e = 1, fe_f%m%ne
            compt(1) = compt(1) +fe_f%m%el_n(e)       ! for each element, the number of lines is added
            do i = 1, fe_f%m%el_t(e)
               compt(2) = compt(2) + fe_f%m%el_t(e)   ! " ", for each node, the number of contributions is added
            enddo
         enddo
         mat%nvar = compt(1)
         mat%nt   = compt(2)
         if (.not.allocated(mat%eltvar)) allocate( mat%eltvar(mat%nvar  ) )
         if (.not.allocated(mat%a_elt )) allocate( mat%a_elt( mat%nt    ) )
         if (.not.allocated(mat%eltptr)) allocate( mat%eltptr(mat%ne +1 ) )
      endif

      ! check of precomputed tables allocation
      if (.not.allocated(fe_f%prc%vcal)) call init_prc_tab(fe_f)

      ! convergence is false at zero iteration
      conv = .false.
      it   = 0
      ! solution loop
      do
         if (VERBOSE >= 20) write(OPU,*) '   Loop FE ', it
         if (conv) exit

         ! assembly of the system
         call assembly_FE_film_reynolds(fe_f, mat, ass_c)
         if (VERBOSE >= 30) write(OPU,*) '   System assembled, thread ', omp_get_thread_num()

         if (ass_c == ASS) then
            ! some stuff can be saved here, provided the reloading of jptr, irow, ... (instead of convert_matrice_format)
            call convert_matrice_format(mat)
            if (VERBOSE >= 30) write(OPU,*) '   Matrix converted, thread ', omp_get_thread_num()
         endif

         if (mat%first) then
            call solve_syst(mat = mat, &
                           step = 'ana')
            mat%first = .false.
            if (VERBOSE >= 10) write(OPU,*) '   Matrix analyzed, thread ', omp_get_thread_num()
         endif

         ! solution of the system
         if (ass_c == ASS) then
            call solve_syst(mat = mat, &
                           step = 'fac')
            if (VERBOSE >= 30) write(OPU,*) '   Matrix factorized, thread ', omp_get_thread_num()
         endif

         call solve_syst(mat = mat, &
                        step = 'sol')
         if (VERBOSE >= 30) write(OPU,*) '   System solved, thread ', omp_get_thread_num()
         if (ass_c == ASS) then
            call solve_syst(mat = mat, &
                           step = 'fre')
            if (VERBOSE >= 30) write(OPU,*) '   Matrix factors freed, thread ', omp_get_thread_num()
         endif

         ! error computation
         error = (sum(mat%x ** 2) / sum(fe_f%vn(:, P_N) ** 2)) ** (0.5_R8)
         it = it + 1
         if (VERBOSE >= 20) write(OPU,*) '   Error ', error

         ! convergence check
         if (error <= fe_f%num_p%eps) conv = .true.

         ! update of variables
         if (fe_f%data_f%fl%fluid_type == MIXT) then
            do i = 1, fe_f%m%n
               if (mat%x(i) < 0.) then
                  fe_f%vn(i, RHO_N) = fe_f%vn(i, RHO_N) + fe_f%vn(i, DRHODP_N) * mat%x(i) * fe_f%num_p%relax
                  if (fe_f%vn(i, RHO_N) < 0.) then
                      fe_f%vn(i,   P_N) = fe_f%data_f%fl%P_0/100
                  else
                      fe_f%vn(i,   P_N) = fe_f%data_f%fl%pres( fe_f%vn(i, RHO_N), & !
                                                               fe_f%vn(i,   T_N) )
                  endif
               else
                  fe_f%vn(i, P_N) = fe_f%vn(i, P_N) + mat%x(i) * fe_f%num_p%relax
               endif
            enddo
         else
            fe_f%vn(:, P_N) = fe_f%vn(:, P_N) + mat%x * fe_f%num_p%relax
         endif

         ! check pressure
         if ( fe_f%data_f%fl%fluid_type == GP ) then
            if (minval(fe_f%vn(:, P_N)) < 0._R8) then
               if (VERBOSE >= 30) write(OPU,*) 'P negative'
            endif
            where (fe_f%vn(:, P_N) < 0._R8) fe_f%vn(:, P_N) = fe_f%data_f%fl%p_0 / 1.e3_R8
         endif

         do i = 1, fe_f%m%n
            fe_f%vn(i, RHO_N)    = fe_f%data_f%fl%rho(   fe_f%vn(i, P_N) , fe_f%vn(i, T_N))
            fe_f%vn(i, DRHODP_N) = fe_f%data_f%fl%drhodp(fe_f%vn(i, P_N) , fe_f%vn(i, T_N))
         enddo

         if (it >= fe_f%num_p%it_max) then
            conv = .true.
            if (VERBOSE >= 20) write(OPU,*) 'maximum number of iteration reached before convergence'
         endif
      enddo

   return
   endsubroutine solve_FE_film

   !=========================================================================================
   !< @note Subroutine to apply boundary conditions on a [[FE_FILM]]
   !<
   !-----------------------------------------------------------------------------------------
   subroutine apply_bc_FE_film(fe_f, bc)
   implicit none
   type(FE_FILM), intent(inout)                     :: fe_f    !! *FE film data type*
   real(kind=R8), intent(in   ), dimension(MAX_NNC) :: bc      !! *boundary conditions at the corners*

      integer(kind=I4) :: e, i, ind, ind_var
      logical(kind=I4) :: mixture
      real(kind=R8)    :: v, v1, v2, x1, x2, y1, y2, l, d

      mixture = (fe_f%data_f%fl%fluid_type == MIXT)

      ! case of hd problem: unknown is p
      if (fe_f%data_f%pb_type == HD) ind_var = P_N

      ! copy of the bc values on the node values
      do i = 1, fe_f%m%nc
         fe_f%vn(fe_f%m%cor(i),ind_var) = fe_f%vn(fe_f%m%cor(i),ind_var) + bc(i)
      enddo

      ! loop on edges
      if (mixture) then
         v1 = bc(1)/fe_f%vn(fe_f%m%cor(1), ind_var)
      else
         v1 = bc(1)
      endif

      x1 = fe_f%m%x(fe_f%m%cor(1))
      y1 = fe_f%m%y(fe_f%m%cor(1))
      do e = fe_f%m%ned, 1, -1
         v2 = v1
         x2 = x1
         y2 = y1

         if (mixture) then
            v1 = bc(e)/fe_f%vn(fe_f%m%cor(e),ind_var)
         else
            v1 = bc(e)
         endif

         x1 = fe_f%m%x(fe_f%m%cor(e))
         y1 = fe_f%m%y(fe_f%m%cor(e))
         ! length of the edge
         l = ((x2 - x1)**2 + (y2 - y1) ** 2) ** 0.5_R8
         do i = 2, fe_f%m%ed(e)%n - 1
            ind = fe_f%m%ed(e)%nm(i)
            ! distance to point 2
            d = ((x2 - fe_f%m%x(ind))**2 + (y2 - fe_f%m%y(ind)) ** 2) ** 0.5_R8
            d = d / l
            v = d * v1 + (1._R8 - d) * v2
            ! linear distribution along the edge
            if (mixture) then
               fe_f%vn(ind, ind_var) = exp( log(fe_f%vn(ind, ind_var)) + v)
            else
               fe_f%vn(ind, ind_var) = fe_f%vn(ind, ind_var) + v
            endif
         enddo

         do i = 1, fe_f%m%ed(e)%n
            ind = fe_f%m%ed(e)%nm(i)
            ! boundary nodes are set as imposed
            fe_f%bc(ind, REY) = 0
         enddo
      enddo

   return
   endsubroutine apply_bc_FE_film

   !=========================================================================================
   !< @note Subroutine to apply boundary conditions on a [[FE_FILM]]
   !<
   !-----------------------------------------------------------------------------------------
   subroutine apply_bc_FE_film_simple(fe_f, bc)
   implicit none
   type(FE_FILM), intent(inout)                     :: fe_f    !! *FE film data type*
   real(kind=R8), intent(in   ), dimension(MAX_NNC) :: bc      !! *boundary conditions at the corners*

      integer(kind=I4) :: e, i, ind, ind_var
      logical(kind=I4) :: mixture
      real(kind=R8)    :: v, v1, v2, x1, x2, y1, y2, l, d

      mixture = (fe_f%data_f%fl%fluid_type == MIXT)

      ! case of hd problem: unknown is p
      if (fe_f%data_f%pb_type == HD) ind_var = P_N

      ! copy of the bc values on the node values
      do i = 1, fe_f%m%nc
         fe_f%vn(fe_f%m%cor(i), ind_var) = bc(i)
      enddo

      ! loop on edges
      if (mixture) then
         v1 = log(fe_f%vn(fe_f%m%cor(1), ind_var))
      else
         v1 = fe_f%vn(fe_f%m%cor(1), ind_var)
      endif

      x1 = fe_f%m%x(fe_f%m%cor(1))
      y1 = fe_f%m%y(fe_f%m%cor(1))
      do e = fe_f%m%ned, 1, -1
         v2 = v1
         x2 = x1
         y2 = y1

         if (mixture) then
            v1 = log(fe_f%vn(fe_f%m%cor(e), ind_var))
         else
            v1 = fe_f%vn(fe_f%m%cor(e), ind_var)
         endif

         x1 = fe_f%m%x(fe_f%m%cor(e))
         y1 = fe_f%m%y(fe_f%m%cor(e))
         ! length of the edge
         l = ((x2 - x1)**2 + (y2 - y1) ** 2) ** 0.5_R8
         do i = 2, fe_f%m%ed(e)%n - 1
            ind = fe_f%m%ed(e)%nm(i)
            ! distance to point 2
            d = ((x2 - fe_f%m%x(ind))**2 + (y2 - fe_f%m%y(ind)) ** 2) ** 0.5_R8
            d = d / l
            v = d * v1 + (1._R8 - d) * v2
            ! linear distribution along the edge
            if (mixture) then
               fe_f%vn(ind, ind_var) = exp(v)
            else
               fe_f%vn(ind, ind_var) = v
            endif
         enddo

         do i = 1, fe_f%m%ed(e)%n
            ind = fe_f%m%ed(e)%nm(i)
            ! boundary nodes are set as imposed
            fe_f%bc(ind, REY) = 0
         enddo
      enddo
   return
   endsubroutine apply_bc_FE_film_simple


   !=========================================================================================
   !< @note Subroutine to to solve a [[FE_FILM]]
   !<
   !-----------------------------------------------------------------------------------------
   subroutine assembly_FE_film_reynolds(fe_f, mat, ass_c)
   implicit none
   type(FE_FILM),    intent(inout) :: fe_f     !! *FE_film data type*
   type(MAT_SOLV),   intent(inout) :: mat      !! *matrices for solving*
   integer(kind=I4), intent(in)    :: ass_c    !! *assembly type*

      integer(kind=I4), dimension(2)                :: compt
      integer(kind=I4), dimension(MAX_NNE)          :: tind4
      real(kind=R8),    dimension(MAX_NNE)          :: b4
      real(kind=R8),    dimension(MAX_NNE, MAX_NNE) :: k4

      integer(kind=I4) :: e, i, ii, el_t, el_n

      ! assembly
      compt(1:2) = 1
      mat%b = 0._R8

      do e = 1, fe_f%m%ne
         call elementary_assembly_FE_film_reynolds(fe_f = fe_f,   & !
                                                  ke_ij = k4,     & ! elementary matrix         : OUT
                                                   be_i = b4,     & ! elementary rhs member     : OUT
                                                  ind_e = tind4,  & ! elementary index member   : OUT
                                                      e = e,      & ! element number
                                                  ass_c = ass_c)    ! assembly type
         do i = 1, 4
            ii = tind4(i)
            mat%b(ii) = mat%b(ii) + b4(i)
         enddo

         el_t = fe_f%m%el_t(e)
         el_n = fe_f%m%el_n(e)

         ! assembly of the elemental matrix in the solver matrix
         if (ass_c == ASS) then
            call assemble_in_mat_sol(mat = mat,    &
                                     num = e,      &
                                    nelt = el_t,   &
                                   nline = el_n,   &
                                    tind = tind4,  &
                                   m_elt = k4,     &
                                   compt = compt)
         endif

      enddo

      if (ass_c == ASS) mat%eltptr(1) = 1
   return
   endsubroutine assembly_FE_film_reynolds

   !=========================================================================================
   !< @note Subroutine to assemble elemental matrices in the solver matrices
   !<
   !-----------------------------------------------------------------------------------------
   subroutine assemble_in_mat_sol(mat, num, nelt, nline, tind, m_elt, compt)
   implicit none
   type(MAT_SOLV),   intent(inout)                        :: mat     !! *mat_solv type*
   integer(kind=I4), intent(in   )                        :: num     !! *element number*
   integer(kind=I4), intent(in   )                        :: nelt    !! *size of elemental matrix*
   integer(kind=I4), intent(in   )                        :: nline   !! *number of lines*
   integer(kind=I4), intent(in   ), dimension(nelt)       :: tind    !! *index table of elemental matrix*
   real(kind=R8),    intent(in   ), dimension(nelt, nelt) :: m_elt   !! *elemental matrix*
   integer(kind=I4), intent(inout), dimension(2)          :: compt   !! *number to index the position in the solver matrix*

      integer(kind=I4) :: i, j

      ! check compt
      if (compt(2) > mat%nt) stop 'compt > matrix size in assemble_in_mat_sol'

      do j = 1, nline
         mat%eltvar(compt(1)) = tind(j)
         compt(1) = compt(1) +1
      enddo
      mat%eltptr(num+1) = compt(1)

      do j = 1, nelt
         do i = 1, nelt
            mat%a_elt(compt(2)) = m_elt(i, j)
            compt(2) = compt(2) + 1
         enddo
      enddo
   return
   endsubroutine assemble_in_mat_sol


   !=========================================================================================
   !< @note Subroutine to solve a [[FE_FILM]]
   !<
   !-----------------------------------------------------------------------------------------
   subroutine elementary_assembly_FE_film_reynolds(fe_f, ke_ij, be_i, ind_e, e, ass_c)
   implicit none
   type(FE_FILM),    intent(inout)                              :: fe_f      !! *FE film*
   real(kind=R8),    intent(out  ), dimension(MAX_NNE, MAX_NNE) :: ke_ij     !! *elementary matrix*
   real(kind=R8),    intent(out  ), dimension(MAX_NNE)          :: be_i      !! *elementary rhs member*
   integer(kind=I4), intent(out  ), dimension(MAX_NNE)          :: ind_e     !! *elementary index member*
   integer(kind=I4), intent(in   )                              :: e         !! *element number*
   integer(kind=I4), intent(in   )                              :: ass_c     !! *assembly type*

      real(kind=R8),    dimension(MAX_NNG, MAX_NNG) :: vni4x, vni4y
      integer(kind=I4), dimension(MAX_NNE)          :: con

      integer(kind=I4) :: i, j, ii, jj, ng, ne

      ng = fe_f%prc%ng
      ne = fe_f%m%el_t(e)
      ! check of precomputed tables allocation
      !~ if (.not.allocated(fe_f%prc%vcal)) call init_prc_tab(fe_f)
      ! assembly
      call compute_prc_tables_reynolds_supg(fe_f, e)

      con(1:ne) = fe_f%m%con(e, 1:ne)

      ke_ij = 0._R8
      ind_e = 0
      be_i  = 0._R8
      do i = 1, ne
         ii = con(i)
         ind_e(i) = ii
         ! case of nodes on the boundary
         if  ((fe_f%bc(ii, REY) == 0) .and. (ass_c /= NO_BC)) then
            ke_ij(i, i)=1.0e10
          ! case of nodes where the pressure is unknown
         else

            vni4x(1:ng, 1:ng) = fe_f%prc%vni4x(i, 1:ng, 1:ng)
            vni4y(1:ng, 1:ng) = fe_f%prc%vni4y(i, 1:ng, 1:ng)

            be_i(i) = be_i(i) - sum( vni4x(1:ng, 1:ng) * (fe_f%prc%vcal(6, 1:ng, 1:ng) - fe_f%prc%vcal(8, 1:ng, 1:ng)) ) & !
                              - sum( vni4y(1:ng, 1:ng) * (fe_f%prc%vcal(7, 1:ng, 1:ng) - fe_f%prc%vcal(9, 1:ng, 1:ng)) )
            do j = 1, fe_f%m%el_t(e)
               jj = con(j)
               ke_ij(i, j) = sum( (vni4x(1:ng, 1:ng) * fe_f%prc%vni4x(j, 1:ng, 1:ng)   &
                                  +vni4y(1:ng, 1:ng) * fe_f%prc%vni4y(j, 1:ng, 1:ng))* fe_f%prc%vcal(2, 1:ng, 1:ng) ) &
                           + sum( (vni4x(1:ng, 1:ng) * fe_f%prc%vcal(10, 1:ng, 1:ng)    &
                                 + vni4y(1:ng, 1:ng) * fe_f%prc%vcal(11, 1:ng, 1:ng))* fe_f%prc%vni4(j, 1:ng, 1:ng))* fe_f%vn(jj, DRHODP_N) &
                           - sum( (vni4x(1:ng, 1:ng) * fe_f%prc%vcal( 3, 1:ng, 1:ng)                             &
                                 + vni4y(1:ng, 1:ng) * fe_f%prc%vcal( 4, 1:ng, 1:ng))                            &
                                 * fe_f%prc%vni4d(j, 1:ng, 1:ng) * fe_f%vn(jj, DRHODP_N) ) !better convergence
              enddo
          endif
      enddo
   return
   endsubroutine elementary_assembly_FE_film_reynolds


   !=========================================================================================
   !< @note Calculate the elementary matrices on a full domain
   !<
   !-----------------------------------------------------------------------------------------
   subroutine elementary_full_domain_FE_film_reynolds(fe_f, mat, ke_ij, be_i, ind_e)
   implicit none
   type(FE_FILM),    intent(inout)                              :: fe_f   !! *FE film*
   type(MAT_SOLV),   intent(inout)                              :: mat    !! *solver type matrices*
   real(kind=R8),    intent(out  ), dimension(MAX_NNC, MAX_NNC) :: ke_ij  !! *elementary matrix*
   real(kind=R8),    intent(inout), dimension(MAX_NNC         ) :: be_i   !! *elementary rhs member*
   integer(kind=I4), intent(out  ), dimension(MAX_NNC         ) :: ind_e  !! *elementary index member*

      integer(kind=I4) :: e, ee, i, j, ii, nbc, nc
      logical(kind=I4) :: flag

      integer(kind=I4), dimension(MAX_NBS) :: ind_n
      real(kind=R8),    dimension(MAX_NBS) :: sav_vn
      real(kind=R8),    dimension(MAX_NNC) :: bc, dp

      nc = fe_f%m%nc

      ! check of precomputed tables allocation
      if (.not.allocated(fe_f%prc%vcal)) call init_prc_tab(fe_f)

      ! boundary condition initialization
      if (BC_SPLINE) then
         bc(1:nc)  = 0
      else
         bc(1:nc) = be_i(1:nc)
      endif

      j = 1
      do e = 1, fe_f%m%ned
      do ee = 1, fe_f%m%ed(e)%n
         ind_n(j)  = fe_f%m%ed(e)%nm(ee)
         sav_vn(j) = fe_f%vn(ind_n(j), P_N)
         j = j +1
      enddo
      enddo
      nbc = j -1

      ! matrices
      ke_ij = 0._R8
      ind_e = 0
      be_i  = 0._R8

      ! calculation of the elementary matrix
      do i = 1, nc ! number of corners
         ii = fe_f%m%cor(i)
         ind_e(i) = ii

         dp(i) = fe_f%vn(ii, P_N)/1000
         if (BC_SPLINE) then
            bc(i) = dp(i)
         else
            bc(i) = bc(i) + dp(i)
         endif

         if (i==1) flag = .true.

         call solve_FE_film(fe_f = fe_f,     &  !
                             mat = mat,      &  !
                              bc = bc,       &  ! in
                        flag_ass = flag)        !

         flag = .false.

         call compute_corner_fluxes(fe_f = fe_f,      &  !
                                     mat = mat,       &  !
                                      bf = be_i)         ! out

         ke_ij(1:nc, i) = - be_i(1:nc) / dp(i)
         if (BC_SPLINE) then
            bc(i) = 0
         else
            bc(i) = bc(i) -dp(i)
         endif

         do j = 1, nbc
            fe_f%vn(ind_n(j), P_N) = sav_vn(j)
         enddo

      enddo

      ! calculation of the RHS table
      call solve_FE_film(fe_f = fe_f,     &  !
                          mat = mat,      &  !
                           bc = bc,       &  !
                     flag_ass = .false.)     !

      call compute_corner_fluxes(fe_f = fe_f,      &  !
                                  mat = mat,       &  !
                                   bf = be_i)         ! out

      ! update of the matrix: the derivative is the difference
      ! between the residuals divided by the delta p
      do i = 1, nc
         ke_ij(1:nc, i) = ke_ij(1:nc, i) + be_i(1:nc)/dp(i)
      enddo
   return
   endsubroutine elementary_full_domain_FE_film_reynolds

   !=========================================================================================
   !< @note Subroutine to calculate the precomputed tables on a FE_film
   !<
   !-----------------------------------------------------------------------------------------
   subroutine init_prc_tab(fe_f)
   implicit none
   type(FE_FILM), intent(inout) :: fe_f         ! *FE film*

      integer(kind=I4) ::  i, j, k, ng
      real(kind=R8), dimension(:), allocatable :: pg, wg

      ! number of Gauss points
      ng = 2

      ! table allocation
      allocate(pg(ng), wg(ng))
      fe_f%prc%ng = ng
      allocate( fe_f%prc%pg(ng), &
                fe_f%prc%wg(ng), &
                fe_f%prc%vni4 (4, ng, ng), &
                fe_f%prc%vni4x(4, ng, ng), &
                fe_f%prc%vni4y(4, ng, ng), &
                fe_f%prc%vni4d(4, ng, ng), &
                fe_f%prc%vcal(20, ng, ng) )
      ! Gauss points definition
      pg(1) = - 1._R8 / sqrt(3._R8)
      pg(2) = - pg(1)
      wg    = 1._R8

      fe_f%prc%pg = pg
      fe_f%prc%wg = wg

      ! precomputation of shape functions and derivatives at the Gauss points
      do i = 1, ng
         do j = 1, ng
            do k = 1, 4
               fe_f%prc%vni4(k, i, j) = ni4(k, pg(i), pg(j))
            enddo
         enddo
      enddo
      deallocate(pg, wg)
   return
   endsubroutine init_prc_tab


   !=========================================================================================
   !> @note Subroutine to calculate the precomputed tables on a [[FE_FILM]]
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine compute_prc_tables_reynolds_supg(fe_f, e)
   implicit none
   type(FE_film),    intent(inout) :: fe_f      !! *FE film*
   integer(kind=I4), intent(in   ) :: e         !! *element number*

      integer(kind=I4) :: i, j, k, ng, nne
      logical(kind=I4) :: gaz
      real(kind=R8), dimension(MAX_NNE) :: tvx, tvy, th, tp, trho, tmu, tdrdp, tx, ty, tnix, tniy
      real(kind=R8), dimension(MAX_NNG) :: wg, pg
      real(kind=R8) :: c5, vc1
      real(kind=R8) :: Pe_ux, Pe_uy, Pe_u, Pe_vx, Pe_vy, Pe_v, coeff_x, coeff_y, kk

      real(kind=R8) :: sx, sy, alpha_u, alpha_ux, alpha_uy, alpha_vx, alpha_vy, alpha_v
      real(kind=R8) :: dNidx_p, dNidy_p, Nid_rho_h, Ni_p, Ni_h, Ni_h3, Ni_vx, Ni_vy, Ni_mu, Ni_rho, Ni_inv_h, Ni_inv_mu

      real(kind=R8) :: lu, wu, lv, wv

      real(kind=R8) :: Ni_drhodp, gradh, gradhx, gradhy, gradp, gradpx, gradpy
      real(kind=R8) :: drdp, h, mu, rho, u, ux, uy, v, vx, vy

      real(kind=R8), dimension(MAX_NNE) :: vni4, vni4x, vni4y, vni4d

      real(kind=R8), dimension(14) :: vcal

      !============================================
      !> {!MUSST/src/inc_doc/Reynolds_discretization.md!}
      !============================================


      !============================================
      !> {!MUSST/src/inc_doc/upwinding_coefficients.md!}
      !============================================

      gaz = (fe_f%data_f%fl%fluid_type==GP)
      nne = fe_f%m%el_t(e)

      ! values on the nodes
      do i = 1, nne
         j = fe_f%m%con(e, i)
         tx(i)    = fe_f%m%x(j)           ! coordinates
         ty(i)    = fe_f%m%y(j)
         tvx(i)   = fe_f%vn (j, VX_N)     ! velocities
         tvy(i)   = fe_f%vn (j, VY_N)
         th(i)    = fe_f%vn (j, H_N)      ! heigth
         tp(i)    = fe_f%vn (j, P_N)      ! pressure
         trho(i)  = fe_f%vn (j, RHO_N)    ! density
         tmu(i)   = fe_f%vn (j, MU_N)     ! viscosity
         tdrdp(i) = fe_f%vn (j, DRHODP_N) ! viscosity derivative regarding P
      enddo

      c5 = dj4(ksi=0._R8, eta=0._R8, x=tx(1:nne), y=ty(1:nne))
      call calc_ni4_xy_derivatives(ni4x = tnix(1:nne), &
                                   ni4y = tniy(1:nne), &
                                    ksi = 0._R8,       &
                                    eta = 0._R8,       &
                                      x = tx(1:nne),   &
                                      y = ty(1:nne),   &
                                     dj = c5)

      ! addition of the groove depth
      th(1:nne) = th(1:nne) + fe_f%vc(e, HG_C)

      ux = sum(tvx(1:nne))/nne
      uy = sum(tvy(1:nne))/nne
      u  = sqrt( ux**2 + uy**2 )

      h    = sum(th   (1:nne))/nne
      drdp = sum(tdrdp(1:nne))/nne
      mu   = sum(tmu  (1:nne))/nne
      rho  = sum(trho (1:nne))/nne

      gradpx = sum( tp(1:nne)*tnix(1:nne) )
      gradpy = sum( tp(1:nne)*tniy(1:nne) )

      gradp = sqrt( gradpx**2 +gradpy**2 )

      gradhx = sum( th(1:nne)*tnix(1:nne) )
      gradhy = sum( th(1:nne)*tniy(1:nne) )

      gradh = sqrt( gradhx**2 +gradhy**2 )

      kk = 6*(drdp/rho)*mu*(1./h**2)

      ! Gauss points
      ng = fe_f%prc%ng

      pg(1:ng) = fe_f%prc%pg(1:ng)
      wg(1:ng) = fe_f%prc%wg(1:ng)

      if (gaz) then

         v  = -(h**2)/(6*mu)
         vx = v*gradpx
         vy = v*gradpy
         v  = sqrt( vx**2 + vy**2 )

         call length_width_elem(spdx = ux, & !  in, speed u along x
                                spdy = uy, & !  in, speed u along y
                                   x = tx, & !  in, element x coordinates
                                   y = ty, & !  in, element y coordinates
                              length = lu, & ! out, length along u
                               width = wu)   ! out,  width

         call length_width_elem(spdx = vx, & !  in, speed v along x
                                spdy = vy, & !  in, speed v along y
                                   x = tx, & !  in, element x coordinates
                                   y = ty, & !  in, element y coordinates
                              length = lv, & ! out, length along v
                               width = wv)   ! out,  width

         ! Peclet number related to u
         Pe_ux = ux * (lu/2) * kk
         Pe_uy = uy * (lu/2) * kk
         Pe_u  = sqrt(Pe_ux**2 +Pe_uy**2)

         ! Peclet number related to v
         Pe_vx = vx * (wu/2) * kk
         Pe_vy = vy * (wu/2) * kk
         Pe_v  = sqrt(Pe_vx**2 +Pe_vy**2)

         alpha_ux = Pe_ux
         alpha_uy = Pe_uy
         alpha_u  = Pe_u

         alpha_v  = (h/1.e6)*(u/1.e2)*(mu/1.e-5) * (wu/lu) * ( Pe_vx * (-uy) + Pe_vy * (+ux) )/u
         alpha_vx = alpha_v * (-uy)/u
         alpha_vy = alpha_v * (+ux)/u
         alpha_v  = sqrt(alpha_vx**2 +alpha_vy**2)

         coeff_x = alpha_ux +alpha_vx ; sx = sign(1._R8, coeff_x) ; coeff_x = abs(coeff_x)
         coeff_y = alpha_uy +alpha_vy ; sy = sign(1._R8, coeff_y) ; coeff_y = abs(coeff_y)

         fe_f%vc(e, PEK_C) = sx*coeff_x
         fe_f%vc(e, PEE_C) = sy*coeff_y

      else

         lu = maxval(tx(1:nne)) -minval(tx(1:nne))
         wu = maxval(ty(1:nne)) -minval(ty(1:nne))
         lv = lu * ux + wu * uy
         Pe_ux = kk * (lv/2)
         Pe_uy = 0._R8
         if (abs(Pe_ux) < 1.e-2_R8) then
            Pe_ux = Pe_ux / 3
         else
            Pe_ux = 1._R8 / (tanh(Pe_ux)) - 1._R8 / Pe_ux
         endif
         alpha_vx = 1._r8
         alpha_vy = 1._r8
         if (u > 0._r8) then
            alpha_vx = lv * ux / (2 * (u ** 2))
            alpha_vy = lv * uy / (2 * (u ** 2))
         endif

         fe_f%vc(e, Pek_c) = alpha_vx
         fe_f%vc(e, Pee_c) = alpha_vy

      endif

      !=============================================
      !> {!MUSST/src/inc_doc/precomputed_integrations.md!}
      !=============================================

      do i = 1, ng
         do j = 1, ng

            vni4(1:nne) = fe_f%prc%vni4(1:nne, i, j)

            ! computation of the shape function derivatives
            c5 = dj4(ksi=pg(i), eta=pg(j), x=tx(1:nne), y=ty(1:nne))
            call calc_ni4_xy_derivatives(ni4x = vni4x(1:nne), &
                                         ni4y = vni4y(1:nne), &
                                          ksi = pg(i),        &
                                          eta = pg(j),        &
                                            x = tx(1:nne),    &
                                            y = ty(1:nne),    &
                                           dj = c5)

            if (c5 < 0) stop 'compute_prc_tables_reynolds_supg: jacobian negative for elt'

            if (gaz) then
               do k = 1, nne
                  vni4d(k) = ni4_up_2d(k, pg(i), pg(j), (/coeff_x, coeff_y/), (/sx, sy/))
               enddo
            else
               vni4d(1:nne) = vni4(1:nne) -alpha_vx*Pe_ux*vni4x(1:nne) -alpha_vy*Pe_ux*vni4y(1:nne)
            endif

            fe_f%prc%vni4x(1:nne, i, j) = vni4x(1:nne)
            fe_f%prc%vni4y(1:nne, i, j) = vni4y(1:nne)

            fe_f%prc%vni4d(1:nne, i, j) = vni4d(1:nne)

            ! computation of the coefficients for the Reynolds equation
            vc1 = wg(i) * wg (j) * c5

            dNidx_p     = sum(vni4x(1:nne) * tp(1:nne))
            dNidy_p     = sum(vni4y(1:nne) * tp(1:nne))

            Nid_rho_h   = sum(vni4d(1:nne) * trho(1:nne) * th(1:nne))

            Ni_p        = sum( vni4(1:nne) * tp(1:nne))
            Ni_h        = sum( vni4(1:nne) * th(1:nne))
            Ni_h3       = sum( vni4(1:nne) * th(1:nne)**3)

            Ni_vx       = sum( vni4(1:nne) * tvx(1:nne))
            Ni_vy       = sum( vni4(1:nne) * tvy(1:nne))
            Ni_mu       = sum( vni4(1:nne) * tmu(1:nne))
            Ni_rho      = sum( vni4(1:nne) * trho(1:nne))

            Ni_inv_h    = sum( vni4(1:nne) * (1._R8 / th(1:nne)))
            Ni_inv_mu   = sum( vni4(1:nne) * (1._R8 / tmu(1:nne)))

            Ni_drhodp   = sum( vni4(1:nne) * tdrdp(1:nne))


            vcal( 1) = vc1
            vcal( 2) = vc1 * Ni_h3 * Ni_rho * Ni_inv_mu

            vcal( 3) = vc1 * 6*Ni_vx * Ni_h
            vcal( 4) = vc1 * 6*Ni_vy * Ni_h

            vcal( 5) = 0

            vcal( 6) = vc1 * Ni_h3 * Ni_inv_mu * Ni_rho * dNidx_p
            vcal( 7) = vc1 * Ni_h3 * Ni_inv_mu * Ni_rho * dNidy_p

            vcal( 8) = vc1 * 6*Ni_vx * Nid_rho_h
            vcal( 9) = vc1 * 6*Ni_vy * Nid_rho_h

            vcal(10) = vc1 * Ni_h3 * Ni_inv_mu * dNidx_p
            vcal(11) = vc1 * Ni_h3 * Ni_inv_mu * dNidy_p

            vcal(12) = vc1 * Ni_p

            vcal(13) = vc1 * (-dNidx_p * Ni_h/2 - Ni_mu * Ni_vx * Ni_inv_h)
            vcal(14) = vc1 * (-dNidy_p * Ni_h/2 - Ni_mu * Ni_vy * Ni_inv_h)

            fe_f%prc%vcal(1:14, i ,j) = vcal(1:14)
         enddo
      enddo

      !=============================================
      !> {!MUSST/css/button.html!}
      !=============================================
   return
   endsubroutine compute_prc_tables_reynolds_supg


   !=========================================================================================
   !> @note Subroutine to calculate the characteristic length of an element, giving
   !        \(x\) and \(y\) velocities
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine length_width_elem(spdx, spdy, x, y, length, width)
   implicit none
   real(kind=R8), intent(in )               :: spdx   !! *fluid velocity along \(x\) axis*
   real(kind=R8), intent(in )               :: spdy   !! *fluid velocity along \(y\) axis*
   real(kind=R8), intent(in ), dimension(4) :: x      !! *corner abscissae*
   real(kind=R8), intent(in ), dimension(4) :: y      !! *corner ordinates*
   real(kind=R8), intent(out)               :: length !! *fluid element length*
   real(kind=R8), intent(out)               :: width  !! *fluid element width*

      real(kind=R8) ::  abx, adx, cbx, cdx,               &
                        aby, ady, cby, cdy,               &
                        acx, bdx, acy, bdy,               &
                        el_S, px, py, pm, spd,            &
                        ac_dot_p, bd_dot_p, ab_dot_p,     &
                        ad_dot_p, cb_dot_p, cd_dot_p,     &
                        x1, x2, x3, x4, y1, y2, y3, y4

      length = 0.
      width  = 0.

      spd = sqrt(spdx**2 +spdy**2)
      if (spd > EPS_R8) then

         px = -spdy  ! direction perpendicular to vector v
         py = +spdx
         pm =  spd

         x1 = x(1) ; y1 = y(1)
         x2 = x(2) ; y2 = y(2)
         x3 = x(3) ; y3 = y(3)
         x4 = x(4) ; y4 = y(4)

         acx = x3 -x1 ; acy = y3 -y1
         bdx = x4 -x2 ; bdy = y4 -y2
         abx = x2 -x1 ; aby = y2 -y1
         adx = x4 -x1 ; ady = y4 -y1
         cbx = x2 -x3 ; cby = y2 -y3
         cdx = x4 -x3 ; cdy = y4 -y3

         ac_dot_p = abs( px*acx +py*acy )
         bd_dot_p = abs( px*bdx +py*bdy )
         ab_dot_p = abs( px*abx +py*aby )
         ad_dot_p = abs( px*adx +py*ady )
         cb_dot_p = abs( px*cbx +py*cby )
         cd_dot_p = abs( px*cdx +py*cdy )

         ! quadrangle projection on the perpendicular to v
         width = max(ac_dot_p, bd_dot_p, ab_dot_p, &
                     ad_dot_p, cb_dot_p, cd_dot_p)/pm

         ! quadrangle area
         el_S = 0.5*( abs(cbx*cdy -cdx*cby) +abs(abx*ady -adx*aby) )

         ! element length
         length = el_S/width

      endif
   return
   endsubroutine length_width_elem


   !=========================================================================================
   !< @note Subroutine to calculate the fluxes at the corner of the domain
   !<
   !-----------------------------------------------------------------------------------------
   subroutine compute_corner_fluxes(fe_f, mat, bf)
   implicit none
   type(FE_FILM),  intent(inout)                     :: fe_f      !! *FE film*
   type(MAT_SOLV), intent(inout)                     :: mat       !! *solver type matrices*
   real(kind=R8),  intent(out  ), dimension(MAX_NNC) :: bf        !! *table of fluxes at the corner*

      integer(kind=I4) :: ind2, ind1, ind, e, i
      real(kind=R8)    :: x1, x2, y1, y2, d, l

      ! calculation of the fluxes on the boundaries by assembly ass_c = no_bc
      call assembly_FE_film_reynolds(fe_f = fe_f,     & !
                                      mat = mat,      & !
                                    ass_c = NO_BC)      !

      ! projection of the fluxes on the corners
      ! loop on edges
      bf   = 0.0_R8
      ind1 = 1
      x1 = fe_f%m%x(fe_f%m%cor(1))
      y1 = fe_f%m%y(fe_f%m%cor(1))
      do e = fe_f%m%ned, 1, -1
         ind2 = ind1
         x2 = x1
         y2 = y1
         ind1 = e
         x1 = fe_f%m%x(fe_f%m%cor(e))
         y1 = fe_f%m%y(fe_f%m%cor(e))
         ! length of the edge
         l = ((x2 - x1)**2 + (y2 - y1) ** 2) ** 0.5_R8
         ! loop on all the nodes - 1 on the edge
         ! without -1, the corner contribution is counted two times
         do i = 1, fe_f%m%ed(e)%n - 1
            ind = fe_f%m%ed(e)%nm(i)
            ! distance to point 2
            d = ((x2 - fe_f%m%x(ind))**2 + (y2 - fe_f%m%y(ind)) ** 2) ** 0.5_R8
            d = d / l
            ! linear projection along the edge of the local flux
            bf(ind1) = bf(ind1) + d * mat%b(ind)
            bf(ind2) = bf(ind2) + (1._R8 - d) * mat%b(ind)
         enddo
      enddo
   return
   endsubroutine compute_corner_fluxes


   !=========================================================================================
   !< @note function to calculate the generated load in a fluid film
   !< \( \int_\Omega p d\Omega \)
   !<
   !-----------------------------------------------------------------------------------------
   real(kind=R8) function fz(fe_f)
   implicit none
   class(FE_film), intent(inout) :: fe_f !! *fluid type*

      integer(kind=I4) :: e
      fz = 0.0_R8
      do e = 1, fe_f%m%ne
          call compute_prc_tables_reynolds_supg(fe_f, e)
          fz = fz + sum(fe_f%prc%vcal(12, :, :))
      enddo
   return
   endfunction fz


   !=========================================================================================
   !< @note function to calculate the friction force along \(x\) in a fluid film
   !< \( \int_\Omega \tau_{xz} d\Omega \)
   !<
   !-----------------------------------------------------------------------------------------
   real(kind=R8) function fx(fe_f)
   implicit none
   class(FE_FILM), intent(inout) :: fe_f !! *FE film*

      integer(kind=I4) :: e
      fx = 0._R8
      do e = 1, fe_f%m%ne
         call compute_prc_tables_reynolds_supg(fe_f, e)
         fx = fx + sum(fe_f%prc%vcal(13, :, :))
      enddo
   return
   endfunction fx


   !=========================================================================================
   !< @note function to calculate the friction force along \(y\) in a fluid film
   !< \( \int_\Omega \tau_{yz} d\Omega \)
   !<
   !-----------------------------------------------------------------------------------------
   real(kind=R8) function fy(fe_f)
   implicit none
   class(FE_film), intent(inout) :: fe_f !! *FE film*

      integer(kind=I4) :: e
      fy = 0._R8
      do e = 1, fe_f%m%ne
         call compute_prc_tables_reynolds_supg(fe_f, e)
         fy = fy + sum(fe_f%prc%vcal(14, :, :))
      enddo
   return
   endfunction fy


   !=========================================================================================
   !> @note Subroutine to save a [[FE_FILM]] in a ```.sur``` file
   !  @endnote
   !-----------------------------------------------------------------------------------------
   subroutine save_fe_field(fe_f, file_name, code, nodal)
   implicit none
   type(FE_FILM),    intent(in) :: fe_f
   character(len=*), intent(in) :: file_name
   integer(kind=I4), intent(in) :: code
   logical(kind=I4), intent(in) :: nodal !! if false : cell value, if true : nodal value

      integer(kind=I4) :: nx, ny
      integer(kind=I4), dimension(1) :: i1, i2
      real(kind=R8), allocatable, dimension(:,:) :: tab_s
      real(kind=R8) :: lx, ly
      character(len=8) :: unit_z

      nx = fe_f%m%nx
      ny = fe_f%m%ny
      lx = fe_f%m%lx
      ly = fe_f%m%ly

      call empty(unit_z)

      if (.not.nodal) then
         nx = nx -1
         ny = ny -1
      endif

      allocate( tab_s(1:nx, 1:ny) ) ; tab_s = -1.

      if (nodal) then
         tab_s = reshape( fe_f%vn(:, code), (/nx, ny/) )
         i1 = index(fe_f%vn_name(code), '(') +1
         i2 = index(fe_f%vn_name(code), ')') -1
         unit_z = fe_f%vn_name(code)(i1(1):i2(1))
      else
         tab_s = reshape( fe_f%vc(:, code), (/nx, ny/) )
         i1 = index(fe_f%vc_name(code), '(') +1
         i2 = index(fe_f%vc_name(code), ')') -1
         unit_z = fe_f%vc_name(code)(i1(1):i2(1))
      endif

      call init_scal(scal = scal_tmp, & !
                       nx = nx,       & !
                       ny = ny,       & !
                       lx = lx,       & ! default unit : m
                       ly = ly,       & !
                   unit_z = unit_z    ) !

      call write_surf(nom_fic = file_name,  & !
                        tab_s = tab_s,      & !
                         scal = scal_tmp    )

      deallocate(tab_s)
   return
   endsubroutine save_fe_field

endmodule film

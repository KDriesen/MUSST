
!! author: Noël Brunetière<br/>&emsp;Arthur Francisco
!! version: 1.0.0
!! date: March,22 2017
!! summary: Definition of a finite element mesh

!< <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<       **Definition of a finite element mesh**
!< </span>

!<### FE_edge type
   !< The edge type define a line which is the boundary of the domain. It contains a number of node, of elements and a connectivity table.
   !< A table is used to give a link with the 2D domain nodes number
!<### Fe_mesh type
   !< The mesh is a 2D mesh (only 4 nodes quadrangles in the current version)
!<### Mesh generation
   !< It is possible to create a structured rectanguler mesh with \ref module_fe_mesh::create_rect_x_ymesh

module mesh
use data_arch, only : I4, R8
implicit none

private

! FE_edge type definition
type FE_EDGE
   integer(kind=I4) :: n                                 !! *number of nodes*
   integer(kind=I4) :: ne                                !! *number of elements*
   integer(kind=I4), dimension(:),   allocatable :: nm   !! *mesh node (numbers in the 2d mesh)*
   integer(kind=I4), dimension(:,:), allocatable :: con  !! *connectivity table*
endtype FE_EDGE

! FE_mesh type definition
type FE_MESH
   real(kind=R8)    :: lx, ly                               !! *size of rectangle*
   real(kind=R8)    :: zx, zy                               !! *coordinates of first point*
   integer(kind=I4) :: nx, ny                               !! *number of nodes in \(x\), \(y\) directions*

   integer(kind=I4) :: n                                    !! *number of nodes*
   integer(kind=I4) :: ne                                   !! *number of elements*
   integer(kind=I4) :: ned                                  !! *number of edges*
   integer(kind=I4) :: nc                                   !! *number of corners*
   real(kind=R8),    dimension(:),   allocatable :: x, y, z !! *nodes coordinates*
   integer(kind=I4), dimension(:,:), allocatable :: con     !! *connectivity table*
   integer(kind=I4), dimension(:),   allocatable :: el_t    !! *element type*
   integer(kind=I4), dimension(:),   allocatable :: el_n    !! *element number of lines*
   type(FE_EDGE),    dimension(:),   allocatable :: ed      !! *edges of the mesh*
   integer(kind=I4), dimension(:),   allocatable :: cor     !! *number of the corner node*
endtype FE_MESH

integer(kind=I4), parameter :: MAX_NNE = 4      !! *maximum number of nodes per element*
integer(kind=I4), parameter :: MAX_NNC = 4      !! *maximum number of corners per element*
integer(kind=I4), parameter :: MAX_NNG = 2      !! *maximum number of Gauss points in a direction*
integer(kind=I4), parameter :: MAX_NBS = 512    !! *maximum number of nodes per BS element in a direction*

public :: FE_EDGE, FE_MESH, create_rect_x_ymesh, MAX_NNE, MAX_NNC, MAX_NNG, MAX_NBS

contains
   !=========================================================================================
   !< @note Subroutine to create a rectangular mesh in the \(x\), \(y\) directions
   !<
   !-----------------------------------------------------------------------------------------
   subroutine create_rect_x_ymesh(m)
   implicit none
   type(FE_MESH), intent(inout)   :: m  !! *FE mesh*
      integer(kind=I4) :: i, j, ind, inde
      real(kind=R8)    :: lx, ly, zx, zy
      integer(kind=I4) :: nx, ny

      lx = m%lx
      ly = m%ly
      zx = m%zx
      zy = m%zy
      nx = m%nx
      ny = m%ny

      ! mesh size definition
      ! number of nodes
      m%n = nx * ny
      ! number of elements
      m%ne = (nx - 1) * (ny - 1)
      ! nodes and elements table allocation
      allocate( m%x(m%n), m%y(m%n), m%z(m%n), m%con(m%ne, 4), m%el_t(m%ne), m%el_n(m%ne) )
      ! all the elements are qua_4
      m%el_t = 4 ! 4 nodes
      m%el_n = 4 ! number of different lines
      ! tables initialisation
      m%con = 0
      m%x   = 0._R8
      m%y   = 0._R8
      m%z   = 0._R8
      ! nodes coordinates
      do j = 1, ny
         do i = 1, nx
            ind =  (j - 1) * nx + i
            m%x(ind) = zx + lx * (i - 1) / (nx - 1)
            m%y(ind) = zy + ly * (j - 1) / (ny - 1)
         enddo
      enddo
      ! connectivity table
      do j = 1, ny - 1
         do i = 1, nx -1
            inde =  (j - 1) * (nx - 1) + i
            ind  =  (j - 1) * nx + i
            m%con(inde, 1) = ind
            m%con(inde, 2) = ind + 1
            m%con(inde, 3) = ind + 1 + nx
            m%con(inde, 4) = ind + nx
         enddo
      enddo
      !----
      ! edges definition
      ! number of edges
      m%ned = 4
      ! allocattion of the edges table
      allocate( m%ed(m%ned) )
      ! number of nodes and elements of each edge
      m%ed(1)%n  = nx
      m%ed(1)%ne = nx - 1
      m%ed(2)%n  = ny
      m%ed(2)%ne = ny - 1
      m%ed(3)%n  = nx
      m%ed(3)%ne = nx - 1
      m%ed(4)%n  = ny
      m%ed(4)%ne = ny - 1
      ! allocation of the edges nodes table
      do j = 1, m%ned
         allocate( m%ed(j)%nm(m%ed(j)%n) )
         m%ed(j)%nm = 0
      enddo
      ! nodes number of the edges from the 2D mesh
      do i = 1, nx
         m%ed(1)%nm(i) = i
         m%ed(3)%nm(i) = (ny -1 ) *nx + (nx - i + 1)
      enddo
      do j = 1, ny
         m%ed(2)%nm(j) = nx * j
         m%ed(4)%nm(j) = nx * (ny - j) + 1
      enddo
      ! edges connectivity table allocation and creation
      do j = 1, m%ned
         allocate (m%ed(j)%con(m%ed(j)%ne, 2))
         m%ed(j)%con=0
         do i = 1, m%ed(j)%ne
            m%ed(j)%con(i,1) = m%ed(j)%nm(i)
            m%ed(j)%con(i,2) = m%ed(j)%nm(i) + 1
         enddo
      enddo
      !----
      ! corners of the mesh
      ! number of corners ( = 4)
      m%nc = m%ned
      ! allocation of corner nodes table
      allocate(m%cor(m%nc))
      ! value of the corner nodes table
      do j = 1, m%nc
         m%cor(j) = m%ed(j)%nm(1)
      enddo

   return
   endsubroutine create_rect_x_ymesh

endmodule mesh


!< author: Arthur Francisco
!  version: 1.0.0
!  date: july, 23 2018
!
!  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!     **Routines to handle Digital Surf binary format (.sur)**
!  </span>

module surfile
use, intrinsic :: ISO_C_BINDING, only : C_CHAR, C_NULL_CHAR, C_FLOAT, C_INT, C_SHORT
use data_arch
use sort_arrays, only : sort_real, sort_real_1real, sort_real_2real
implicit none

private

!< <span style="color:green">Surface object: header and heights</span>
! @note
!    Adapted from 'surffile.c', 'gwyddion' software, Copyright (C) 2005 David Necas, Petr Klapetek.
! @endnote
! @warning
!    Must be 512 bytes long
! @endwarning


type SCALE_SURF
   ! bytes below: 8+10+2*12+9*16+2*30+34+128 = 408
   character(len =  12) :: signature
   character(len =  16) :: xlength_unit
   character(len =  16) :: ylength_unit
   character(len =  16) :: zlength_unit
   character(len =  16) :: xaxis
   character(len =  16) :: yaxis
   character(len =  16) :: zaxis
   character(len =  16) :: dx_unit
   character(len =  16) :: dy_unit
   character(len =  16) :: dz_unit
   character(len =  30) :: object_name
   character(len =  30) :: operator_name
   character(len = 128) :: client_zone
   character(len =   8) :: reserved
   character(len =  34) :: reservedzone
   character(len =  12) :: obsolete
   character(len =  10) :: obsolete2

   ! bytes below: 10*4 = 40
   real(kind=R4) :: dx
   real(kind=R4) :: dy
   real(kind=R4) :: dz
   real(kind=R4) :: xunit_ratio
   real(kind=R4) :: yunit_ratio
   real(kind=R4) :: zunit_ratio
   real(kind=R4) :: XOffset
   real(kind=R4) :: YOffset
   real(kind=R4) :: ZOffset
   real(kind=R4) :: measurement_duration

   ! bytes below: 5*4 = 20
   integer(kind=I4) :: zmin
   integer(kind=I4) :: zmax
   integer(kind=I4) :: xres
   integer(kind=I4) :: yres
   integer(kind=I4) :: nofpoints

   ! bytes below: 22*2 = 44
   integer(kind=2) :: format
   integer(kind=2) :: version
   integer(kind=2) :: material_code
   integer(kind=2) :: type
   integer(kind=2) :: range
   integer(kind=2) :: special_points
   integer(kind=2) :: absolute
   integer(kind=2) :: pointsize
   integer(kind=2) :: imprint
   integer(kind=2) :: inversion
   integer(kind=2) :: leveling
   integer(kind=2) :: seconds
   integer(kind=2) :: minutes
   integer(kind=2) :: hours
   integer(kind=2) :: day
   integer(kind=2) :: month
   integer(kind=2) :: year
   integer(kind=2) :: dayof
   integer(kind=2) :: comment_size
   integer(kind=2) :: private_size
   integer(kind=2) :: nobjects
   integer(kind=2) :: acquisition

   real(kind=R8)   :: lx       !! *surface length*
   real(kind=R8)   :: ly       !! *surface width*
   real(kind=R8)   :: lz       !! *surface height (max -min)*

   real(kind=R8)   :: mu       !! *surface mean height*
   real(kind=R8)   :: si       !! *surface mean height*
endtype SCALE_SURF

type OBJ_SURF
   ! bytes below: 8+10+2*12+9*16+2*30+34+128 = 408
   character(kind=C_CHAR), dimension( 12) :: signature
   character(kind=C_CHAR), dimension( 16) :: xlength_unit
   character(kind=C_CHAR), dimension( 16) :: ylength_unit
   character(kind=C_CHAR), dimension( 16) :: zlength_unit
   character(kind=C_CHAR), dimension( 16) :: xaxis
   character(kind=C_CHAR), dimension( 16) :: yaxis
   character(kind=C_CHAR), dimension( 16) :: zaxis
   character(kind=C_CHAR), dimension( 16) :: dx_unit
   character(kind=C_CHAR), dimension( 16) :: dy_unit
   character(kind=C_CHAR), dimension( 16) :: dz_unit
   character(kind=C_CHAR), dimension( 30) :: object_name
   character(kind=C_CHAR), dimension( 30) :: operator_name
   character(kind=C_CHAR), dimension(128) :: client_zone
   character(kind=C_CHAR), dimension(  8) :: reserved
   character(kind=C_CHAR), dimension( 34) :: reservedzone
   character(kind=C_CHAR), dimension( 12) :: obsolete
   character(kind=C_CHAR), dimension( 10) :: obsolete2

   ! bytes below: 10*4 = 40
   real(kind=C_FLOAT) :: dx
   real(kind=C_FLOAT) :: dy
   real(kind=C_FLOAT) :: dz
   real(kind=C_FLOAT) :: xunit_ratio
   real(kind=C_FLOAT) :: yunit_ratio
   real(kind=C_FLOAT) :: zunit_ratio
   real(kind=C_FLOAT) :: XOffset
   real(kind=C_FLOAT) :: YOffset
   real(kind=C_FLOAT) :: ZOffset
   real(kind=C_FLOAT) :: measurement_duration

   ! bytes below: 5*4 = 20
   integer(kind=C_INT) :: zmin
   integer(kind=C_INT) :: zmax
   integer(kind=C_INT) :: xres
   integer(kind=C_INT) :: yres
   integer(kind=C_INT) :: nofpoints

   ! bytes below: 22*2 = 44
   integer(kind=C_SHORT) :: format
   integer(kind=C_SHORT) :: version
   integer(kind=C_SHORT) :: material_code
   integer(kind=C_SHORT) :: type
   integer(kind=C_SHORT) :: range
   integer(kind=C_SHORT) :: special_points
   integer(kind=C_SHORT) :: absolute
   integer(kind=C_SHORT) :: pointsize
   integer(kind=C_SHORT) :: imprint
   integer(kind=C_SHORT) :: inversion
   integer(kind=C_SHORT) :: leveling
   integer(kind=C_SHORT) :: seconds
   integer(kind=C_SHORT) :: minutes
   integer(kind=C_SHORT) :: hours
   integer(kind=C_SHORT) :: day
   integer(kind=C_SHORT) :: month
   integer(kind=C_SHORT) :: year
   integer(kind=C_SHORT) :: dayof
   integer(kind=C_SHORT) :: comment_size
   integer(kind=C_SHORT) :: private_size
   integer(kind=C_SHORT) :: nobjects
   integer(kind=C_SHORT) :: acquisition

   integer(kind=C_INT), allocatable :: val(:) !! *heights*
endtype OBJ_SURF

integer(kind=4), parameter :: SURF_DAT = 1 !! *'.dat' format, txt*
integer(kind=4), parameter :: SURF_SUR = 2 !! *'.sur' format, binary*

public :: read_surf, write_surf, init_scal, unit2IUc, unit2IUf, &
          trans_surf_txt, scal2surf, surf2scal, empty, SCALE_SURF, OBJ_SURF

contains

   subroutine scal2surf(scal, surf)
   implicit none
   type(SCALE_SURF), intent(in ) :: scal !! *object [[SCALE_SURF]]*
   type(OBJ_SURF),   intent(out) :: surf !! *object [[OBJ_SURF]]*

      call f_c_string(fs=trim(scal%signature),     &
                      cs=     surf%signature)

      call f_c_string(fs=trim(scal%xlength_unit),  &
                      cs=     surf%xlength_unit)

      call f_c_string(fs=trim(scal%ylength_unit),  &
                      cs=     surf%ylength_unit)

      call f_c_string(fs=trim(scal%zlength_unit),  &
                      cs=     surf%zlength_unit)

      call f_c_string(fs=trim(scal%xaxis),         &
                      cs=     surf%xaxis)

      call f_c_string(fs=trim(scal%yaxis),         &
                      cs=     surf%yaxis)

      call f_c_string(fs=trim(scal%zaxis),         &
                      cs=     surf%zaxis)

      call f_c_string(fs=trim(scal%dx_unit),       &
                      cs=     surf%dx_unit)

      call f_c_string(fs=trim(scal%dy_unit),       &
                      cs=     surf%dy_unit)

      call f_c_string(fs=trim(scal%dz_unit),       &
                      cs=     surf%dz_unit)

      call f_c_string(fs=trim(scal%object_name),   &
                      cs=     surf%object_name)

      call f_c_string(fs=trim(scal%operator_name), &
                      cs=     surf%operator_name)

      call f_c_string(fs=trim(scal%client_zone),   &
                      cs=     surf%client_zone)

      call f_c_string(fs=trim(scal%reserved),      &
                      cs=     surf%reserved)

      call f_c_string(fs=trim(scal%reservedzone),  &
                      cs=     surf%reservedzone)

      call f_c_string(fs=trim(scal%obsolete),      &
                      cs=     surf%obsolete)

      call f_c_string(fs=trim(scal%obsolete2),     &
                      cs=     surf%obsolete2)

      surf%dx                    = scal%dx
      surf%dy                    = scal%dy
      surf%dz                    = scal%dz
      surf%xunit_ratio           = scal%xunit_ratio
      surf%yunit_ratio           = scal%yunit_ratio
      surf%zunit_ratio           = scal%zunit_ratio
      surf%XOffset               = scal%XOffset
      surf%YOffset               = scal%YOffset
      surf%ZOffset               = scal%ZOffset
      surf%measurement_duration  = scal%measurement_duration

      surf%zmin      = scal%zmin
      surf%zmax      = scal%zmax
      surf%xres      = scal%xres
      surf%yres      = scal%yres
      surf%nofpoints = scal%nofpoints

      surf%format          = scal%format
      surf%version         = scal%version
      surf%material_code   = scal%material_code
      surf%type            = scal%type
      surf%range           = scal%range
      surf%special_points  = scal%special_points
      surf%absolute        = scal%absolute
      surf%pointsize       = scal%pointsize
      surf%imprint         = scal%imprint
      surf%inversion       = scal%inversion
      surf%leveling        = scal%leveling
      surf%seconds         = scal%seconds
      surf%minutes         = scal%minutes
      surf%hours           = scal%hours
      surf%day             = scal%day
      surf%month           = scal%month
      surf%year            = scal%year
      surf%dayof           = scal%dayof
      surf%comment_size    = scal%comment_size
      surf%private_size    = scal%private_size
      surf%nobjects        = scal%nobjects
      surf%acquisition     = scal%acquisition
   return
   endsubroutine scal2surf


   subroutine surf2scal(surf, scal)
   implicit none
   type(OBJ_SURF),   intent(in ) :: surf !! *object [[OBJ_SURF]]*
   type(SCALE_SURF), intent(out) :: scal !! *object [[SCALE_SURF]]*
      integer(kind=I4) :: i

      call c_f_string(cs = surf%signature,     &
                      fs = scal%signature,     &
                 borne_s = i)

      call c_f_string(cs = surf%xlength_unit,  &
                      fs = scal%xlength_unit,  &
                 borne_s = i)

      call c_f_string(cs = surf%ylength_unit,  &
                      fs = scal%ylength_unit,  &
                 borne_s = i)

      call c_f_string(cs = surf%zlength_unit,  &
                      fs = scal%zlength_unit,  &
                 borne_s = i)

      call c_f_string(cs = surf%xaxis,         &
                      fs = scal%xaxis,         &
                 borne_s = i)

      call c_f_string(cs = surf%yaxis,         &
                      fs = scal%yaxis,         &
                 borne_s = i)

      call c_f_string(cs = surf%zaxis,         &
                      fs = scal%zaxis,         &
                 borne_s = i)

      call c_f_string(cs = surf%dx_unit,       &
                      fs = scal%dx_unit,       &
                 borne_s = i)

      call c_f_string(cs = surf%dy_unit,       &
                      fs = scal%dy_unit,       &
                 borne_s = i)

      call c_f_string(cs = surf%dz_unit,       &
                      fs = scal%dz_unit,       &
                 borne_s = i)

      call c_f_string(cs = surf%object_name,   &
                      fs = scal%object_name,   &
                 borne_s = i)

      call c_f_string(cs = surf%operator_name, &
                      fs = scal%operator_name, &
                 borne_s = i)

      call c_f_string(cs = surf%client_zone,   &
                      fs = scal%client_zone,   &
                 borne_s = i)

      call c_f_string(cs = surf%reserved,      &
                      fs = scal%reserved,      &
                 borne_s = i)

      call c_f_string(cs = surf%reservedzone,  &
                      fs = scal%reservedzone,  &
                 borne_s = i)

      call c_f_string(cs = surf%obsolete,      &
                      fs = scal%obsolete,      &
                 borne_s = i)

      call c_f_string(cs = surf%obsolete2,     &
                      fs = scal%obsolete2,     &
                 borne_s =i)

      scal%dx                    = surf%dx
      scal%dy                    = surf%dy
      scal%dz                    = surf%dz
      scal%xunit_ratio           = surf%xunit_ratio
      scal%yunit_ratio           = surf%yunit_ratio
      scal%zunit_ratio           = surf%zunit_ratio
      scal%XOffset               = surf%XOffset
      scal%YOffset               = surf%YOffset
      scal%ZOffset               = surf%ZOffset
      scal%measurement_duration  = surf%measurement_duration

      scal%zmin      = surf%zmin
      scal%zmax      = surf%zmax
      scal%xres      = surf%xres
      scal%yres      = surf%yres
      scal%nofpoints = surf%nofpoints

      scal%format          = surf%format
      scal%version         = surf%version
      scal%material_code   = surf%material_code
      scal%type            = surf%type
      scal%range           = surf%range
      scal%special_points  = surf%special_points
      scal%absolute        = surf%absolute
      scal%pointsize       = surf%pointsize
      scal%imprint         = surf%imprint
      scal%inversion       = surf%inversion
      scal%leveling        = surf%leveling
      scal%seconds         = surf%seconds
      scal%minutes         = surf%minutes
      scal%hours           = surf%hours
      scal%day             = surf%day
      scal%month           = surf%month
      scal%year            = surf%year
      scal%dayof           = surf%dayof
      scal%comment_size    = surf%comment_size
      scal%private_size    = surf%private_size
      scal%nobjects        = surf%nobjects
      scal%acquisition     = surf%acquisition
   return
   endsubroutine surf2scal

   !=========================================================================================
   !>@note
   !   Just empties a string
   ! @endnote
   !-----------------------------------------------------------------------------------------
   subroutine empty(charinout)
   implicit none
   character(len=*), intent(inout) :: charinout
      charinout = repeat(' ', len(charinout))
   return
   endsubroutine empty


   !=========================================================================================
   !< @note
   !   Converts a C string to a Fortran string
   !
   !   A From a memory viewpoint, a C string is like a character vector ending with a ```C_NULL_CHAR```, so as long as
   !   it is not found, the characters are copied one by one in a fortran string
   ! @endnote
   !-----------------------------------------------------------------------------------------
   subroutine c_f_string(cs, fs, borne_s)
   implicit none
   character(kind=C_CHAR), dimension(:), intent(in) :: cs   !! *C string*
   character(len=*), intent(out) :: fs                      !! *Fortran string*
   integer(kind=I4), intent(out) :: borne_s                 !! *resulting Fortran string length*
      integer(kind=I4) :: i, ucs
      ucs = size(cs) ! vector length
      borne_s = ucs  ! resulting string default length
      i = 1
      do
         if (i>ucs) exit
         if (cs(i)==C_NULL_CHAR) then  ! fin de chaîne c rencontrée ; s'il n'y a pas de null_char
            borne_s = i-1              ! c'est qu'on utilise tout le vecteur
            exit
         endif
         i = i + 1
      enddo
      call empty(fs)
      do i = 1, borne_s ! the C string is translated into fortran
         fs(i:i) = cs(i)
      enddo
   return
   endsubroutine c_f_string


   !=========================================================================================
   !< @note
   !   Converts a Fortran string to a C string
   !
   !   A From a memory viewpoint, a C string is like a character vector ending with a ```C_NULL_CHAR```,
   !   so the characters are copied one by one in a ```C_CHAR``` vector that ends with a ```C_NULL_CHAR```
   ! @endnote
   !-----------------------------------------------------------------------------------------
   subroutine f_c_string(fs, cs)
   implicit none
   character(len=*), intent(in) :: fs                       !! *fortran string*
   character(kind=C_CHAR), dimension(:), intent(out) :: cs  !! *resulting C string*
      integer(kind=I4) :: i, ufs
      ufs = len_trim(fs)   ! longueur de la chaîne fortran sans les null, les blancs, ...
      if (ufs==0) then     ! si la chaîne est vide
         cs(1) = C_NULL_CHAR
      else
         do i = 1, ufs
            cs(i) = fs(i:i)
         enddo
         if (ufs<size(cs)) cs(ufs+1) = C_NULL_CHAR ! si la fin du vecteur n'est pas atteinte
      endif
   return
   endsubroutine f_c_string


   !=========================================================================================
   !> @note
   !   [[OBJ_SURF]] initialization, every unit is m
   ! @endnote
   !-----------------------------------------------------------------------------------------
   subroutine init_scal(scal, nx, ny, lx, ly, unit_z)
   implicit none
   type(SCALE_SURF), intent(out) :: scal !! *object [[SCALE_SURF]]*
   integer(kind=I4), optional, intent(in) :: nx
   integer(kind=I4), optional, intent(in) :: ny
   real(kind=R8),    optional, intent(in) :: lx
   real(kind=R8),    optional, intent(in) :: ly
   character(*),     optional, intent(in) :: unit_z
      integer(kind=I4), dimension(1:8) :: time_val
      character(len=256) :: string

      call date_and_time(values=time_val)

      scal%format                = 0
      scal%nobjects              = 1
      scal%version               = 1
      scal%type                  = 2
      scal%material_code         = 1
      scal%acquisition           = 0
      scal%range                 = 0
      scal%special_points        = 0
      scal%absolute              = 1
      scal%pointsize             = 32
      scal%zmin                  = 0
      scal%zmax                  = 0
      scal%xres                  = 0
      scal%yres                  = 0
      scal%nofpoints             = 0
      scal%xunit_ratio           = 1.
      scal%yunit_ratio           = 1.
      scal%zunit_ratio           = 1.
      scal%imprint               = 0
      scal%inversion             = 0
      scal%leveling              = 0
      scal%seconds               = time_val(7)
      scal%minutes               = time_val(6)
      scal%hours                 = time_val(5)
      scal%day                   = time_val(3)
      scal%month                 = time_val(2)
      scal%year                  = time_val(1)
      scal%dayof                 = 0
      scal%measurement_duration  = 0.0
      scal%comment_size          = 0
      scal%private_size          = 0
      scal%XOffset               = 0.
      scal%YOffset               = 0.
      scal%ZOffset               = 0.

      call empty(scal%reserved)
      call empty(scal%obsolete2)
      call empty(scal%obsolete)
      call empty(scal%reservedzone)
      call empty(scal%client_zone)

      call empty(scal%object_name)
      call empty(scal%signature)
      call empty(scal%operator_name)
      scal%object_name     = 'HOME MADE'
      scal%signature       = 'DIGITAL SURF'
      scal%operator_name   = 'MOD_SURFILE'

      call empty(scal%xaxis)
      call empty(scal%yaxis)
      call empty(scal%zaxis)
      scal%xaxis  = "X"
      scal%yaxis  = "Y"
      scal%zaxis  = "Z"

      call empty(scal%xlength_unit)
      call empty(scal%ylength_unit)
      call empty(scal%zlength_unit)
      call empty(scal%dx_unit)
      call empty(scal%dy_unit)
      call empty(scal%dz_unit)
      scal%xlength_unit = "m" ; scal%dx_unit = trim(scal%xlength_unit) ; scal%dx = 1.0
      scal%ylength_unit = "m" ; scal%dy_unit = trim(scal%ylength_unit) ; scal%dy = 1.0
      scal%zlength_unit = "m" ; scal%dz_unit = trim(scal%zlength_unit) ; scal%dz = 1.0

      scal%mu = 0
      scal%si = 0

      if (present(nx)) scal%xres = nx
      if (present(ny)) scal%yres = ny

      if (present(nx).and.present(lx)) scal%dx = lx/nx
      if (present(ny).and.present(ly)) scal%dy = ly/ny

      if (present(unit_z)) then ; scal%zlength_unit = trim(unit_z)
                                  scal%dz_unit      = trim(unit_z) ; endif

   return
   endsubroutine init_scal


   !=========================================================================================
   !< @note
   !   Writes an [[OBJ_SURF]] object in a text file
   !
   !   The object components are first written in a fortran string, then it is written into
   !   the file with a comment
   ! @endnote
   !-----------------------------------------------------------------------------------------
   subroutine trans_surf_txt(surf, fichier, xyz)
   implicit none
   type(OBJ_SURF), intent(in)   :: surf      !! *object [[OBJ_SURF]]*
   character(len=*), intent(in) :: fichier   !! *text file to write*
   logical(kind=I4), intent(in) :: xyz       !! *whether to also write the heights (maybe huge)*
      integer(kind=I4) :: i, k, s
      character(len=512) :: string, cc

      call get_unit(k)

      open(k, file=trim(fichier))

         call c_f_string(cs=surf%signature, fs=string, borne_s=s)
         write(cc,*) '"',trim(string(1:s)),'"'  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "signature              " ; call empty(cc)
         write(cc,*) surf%format                ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "format                 " ; call empty(cc)
         write(cc,*) surf%nobjects              ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "nobjects               " ; call empty(cc)
         write(cc,*) surf%version               ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "version                " ; call empty(cc)
         write(cc,*) surf%type                  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "type                   " ; call empty(cc)
         call c_f_string(cs=surf%object_name, fs=string, borne_s=s)
         write(cc,*) '"',trim(string(1:s)),'"'  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "object_name            " ; call empty(cc)
         call c_f_string(cs=surf%operator_name, fs=string, borne_s=s)
         write(cc,*) '"',trim(string(1:s)),'"'  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "operator_name          " ; call empty(cc)
         write(cc,*) surf%material_code         ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "material_code          " ; call empty(cc)
         write(cc,*) surf%acquisition           ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "acquisition            " ; call empty(cc)
         write(cc,*) surf%range                 ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "range                  " ; call empty(cc)
         write(cc,*) surf%special_points        ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "special_points         " ; call empty(cc)
         write(cc,*) surf%absolute              ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "absolute               " ; call empty(cc)
         call c_f_string(cs=surf%reserved, fs=string, borne_s=s)
         write(cc,*) '"',trim(string(1:s)),'"'  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "reserved               " ; call empty(cc)
         write(cc,*) surf%pointsize             ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "pointsize              " ; call empty(cc)
         write(cc,*) surf%zmin                  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "zmin                   " ; call empty(cc)
         write(cc,*) surf%zmax                  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "zmax                   " ; call empty(cc)
         write(cc,*) surf%xres                  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "xres                   " ; call empty(cc)
         write(cc,*) surf%yres                  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "yres                   " ; call empty(cc)
         write(cc,*) surf%nofpoints             ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "nofpoints              " ; call empty(cc)
         write(cc,*) surf%dx                    ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "dx                     " ; call empty(cc)
         write(cc,*) surf%dy                    ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "dy                     " ; call empty(cc)
         write(cc,*) surf%dz                    ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "dz                     " ; call empty(cc)
         call c_f_string(cs=surf%xaxis, fs=string, borne_s=s)
         write(cc,*) '"',trim(string(1:s)),'"'  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "xaxis                  " ; call empty(cc)
         call c_f_string(cs=surf%yaxis, fs=string, borne_s=s)
         write(cc,*) '"',trim(string(1:s)),'"'  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "yaxis                  " ; call empty(cc)
         call c_f_string(cs=surf%zaxis, fs=string, borne_s=s)
         write(cc,*) '"',trim(string(1:s)),'"'  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "zaxis                  " ; call empty(cc)
         call c_f_string(cs=surf%dx_unit, fs=string, borne_s=s)
         write(cc,*) '"',trim(string(1:s)),'"'  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "dx_unit                " ; call empty(cc)
         call c_f_string(cs=surf%dy_unit, fs=string, borne_s=s)
         write(cc,*) '"',trim(string(1:s)),'"'  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "dy_unit                " ; call empty(cc)
         call c_f_string(cs=surf%dz_unit, fs=string, borne_s=s)
         write(cc,*) '"',trim(string(1:s)),'"'  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "dz_unit                " ; call empty(cc)
         call c_f_string(cs=surf%xlength_unit, fs=string, borne_s=s)
         write(cc,*) '"',trim(string(1:s)),'"'  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "xlength_unit           " ; call empty(cc)
         call c_f_string(cs=surf%ylength_unit, fs=string, borne_s=s)
         write(cc,*) '"',trim(string(1:s)),'"'  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "ylength_unit           " ; call empty(cc)
         call c_f_string(cs=surf%zlength_unit, fs=string, borne_s=s)
         write(cc,*) '"',trim(string(1:s)),'"'  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "zlength_unit           " ; call empty(cc)
         write(cc,*) surf%xunit_ratio           ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "xunit_ratio            " ; call empty(cc)
         write(cc,*) surf%yunit_ratio           ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "yunit_ratio            " ; call empty(cc)
         write(cc,*) surf%zunit_ratio           ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "zunit_ratio            " ; call empty(cc)
         write(cc,*) surf%imprint               ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "imprint                " ; call empty(cc)
         write(cc,*) surf%inversion             ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "inversion              " ; call empty(cc)
         write(cc,*) surf%leveling              ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "leveling               " ; call empty(cc)
         call c_f_string(cs=surf%obsolete, fs=string, borne_s=s)
         write(cc,*) '"',trim(string(1:s)),'"'  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "obsolete               " ; call empty(cc)
         write(cc,*) surf%seconds               ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "seconds                " ; call empty(cc)
         write(cc,*) surf%minutes               ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "minutes                " ; call empty(cc)
         write(cc,*) surf%hours                 ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "hours                  " ; call empty(cc)
         write(cc,*) surf%day                   ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "day                    " ; call empty(cc)
         write(cc,*) surf%month                 ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "month                  " ; call empty(cc)
         write(cc,*) surf%year                  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "year                   " ; call empty(cc)
         write(cc,*) surf%dayof                 ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "dayof                  " ; call empty(cc)
         write(cc,*) surf%measurement_duration  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "measurement_duration   " ; call empty(cc)
         call c_f_string(cs=surf%obsolete2, fs=string, borne_s=s)
         write(cc,*) '"',trim(string(1:s)),'"'  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "obsolete2              " ; call empty(cc)
         write(cc,*) surf%comment_size          ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "comment_size           " ; call empty(cc)
         write(cc,*) surf%private_size          ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "private_size           " ; call empty(cc)
         call c_f_string(cs=surf%client_zone, fs=string, borne_s=s)
         write(cc,*) '"',trim(string(1:s)),'"'  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "client_zone            " ; call empty(cc)
         write(cc,*) surf%XOffset               ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "XOffset                " ; call empty(cc)
         write(cc,*) surf%YOffset               ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "YOffset                " ; call empty(cc)
         write(cc,*) surf%ZOffset               ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "ZOffset                " ; call empty(cc)
         call c_f_string(cs=surf%reservedzone, fs=string, borne_s=s)
         write(cc,*) '"',trim(string(1:s)),'"'  ; write(k,'(1x,a,T130,a)') adjustl(trim(cc)), "reservedzone           " ; call empty(cc)

         if (xyz) then
            do i = 0, surf%nofpoints -1
               write(k,*) mod(i, surf%xres)*surf%dx, (i/surf%xres)*surf%dy, surf%val(i+1)*surf%dz
            enddo
         endif

      close(k)

   return
   endsubroutine trans_surf_txt


   !=========================================================================================
   !< @note
   !   Subroutine that opens a ```.sur``` file and transfers it contents into an object [[OBJ_SURF]]
   ! @endnote
   ! @warning
   !    By default here, the heights are not written with ```dump=.true.```
   ! @endwarning
   !-----------------------------------------------------------------------------------------
   subroutine open_surffile(fichier, surf, scal, dump)
   implicit none
   character(len=*), intent(in)  :: fichier         !! *file to be read*
   type(OBJ_SURF), intent(out)   :: surf            !! *object that will contain the file infos and heights*
   type(SCALE_SURF), intent(out) :: scal            !! *object [[SCALE_SURF]]*
   logical(kind=I4), optional, intent(in) :: dump   !! *whether to transform the data in a text file*
      integer(kind=I4) :: i, k
      real(kind=R8)    :: scal_x, scal_y, scal_z
      character(kind=C_CHAR) :: charact

      call get_unit(k)
      open(k , file=trim(fichier),     & !
               form='unformatted',     & !
               access="stream",        & ! beware the "frecord-marker" in other modes
               action="read",          & !
               position="rewind",      & !
               convert='little_endian',& !
               status='old')

         read(k)  surf%signature, surf%format, surf%nobjects, surf%version, surf%type, surf%object_name,                &
                  surf%operator_name, surf%material_code, surf%acquisition, surf%range, surf%special_points,            &
                  surf%absolute, surf%reserved, surf%pointsize, surf%zmin, surf%zmax, surf%xres, surf%yres,             &
                  surf%nofpoints, surf%dx, surf%dy, surf%dz, surf%xaxis, surf%yaxis, surf%zaxis, surf%dx_unit,          &
                  surf%dy_unit, surf%dz_unit, surf%xlength_unit, surf%ylength_unit, surf%zlength_unit,                  &
                  surf%xunit_ratio, surf%yunit_ratio, surf%zunit_ratio, surf%imprint, surf%inversion, surf%leveling,    &
                  surf%obsolete, surf%seconds, surf%minutes, surf%hours, surf%day, surf%month, surf%year, surf%dayof,   &
                  surf%measurement_duration, surf%obsolete2, surf%comment_size, surf%private_size, surf%client_zone,    &
                  surf%XOffset, surf%YOffset, surf%ZOffset, surf%reservedzone

         do i = 1, surf%comment_size
            read(k) charact
         enddo

         allocate( surf%val(1:surf%nofpoints) )
         do i = 1, surf%nofpoints
            read(k) surf%val(i)
         enddo

      close(k)

      call surf2scal(surf, scal)

      if (present(dump).and.dump) call trans_surf_txt(surf, trim(fichier)//'.txt', xyz=.false.)

   return
   endsubroutine open_surffile


   !=========================================================================================
   !< @note
   !   Converts uppercase to lowercase, adapted from [here](http://fortranwiki.org/fortran/show/String_Functions)
   ! @endnote
   !-----------------------------------------------------------------------------------------
   function lower(s1) result (s2)
   character(*), intent(in) :: s1   !! *string to transform to lower case*
   character(len(s1))  :: s2        !! *result: same string but each character is lower case*
      character(len=1) :: ch
      integer(kind=I4), parameter :: duc = ichar('A') - ichar('a')
      integer(kind=I4) :: i
      do i = 1, len(s1)
         ch = s1(i:i)
         if (ch >= 'A'.and.ch <= 'Z') ch = char(ichar(ch)-duc)
         s2(i:i) = ch
      enddo
   return
   endfunction lower

   function unit2IUc(string) result (met)
   implicit none
   real(kind=R8) :: met
   character(kind=C_CHAR), dimension(:), intent(in) :: string
      character(len=2) :: chaine
      chaine = string(1)//string(2)
      met    = unit2IUf(chaine)
   return
   endfunction unit2IUc

   function unit2IUf(string) result (met)
   implicit none
   real(kind=R8) :: met
   character(*), intent(in) :: string
      select case(string)
         case('m ')
            met = 1.e+00_R8
         case('cm')
            met = 1.e-02_R8
         case('mm')
            met = 1.e-03_R8
         case('µm')
            met = 1.e-06_R8
         case('nm')
            met = 1.e-09_R8
         case('pm')
            met = 1.e-12_R8
         case('Pa')
            met = 1.e+00_R8
         case('MP')
            met = 1.e+06_R8
         case('GP')
            met = 1.e+09_R8
         case default
            met = 1.e+00_R8
      endselect
   return
   endfunction unit2IUf

   !=========================================================================================
   !< @note
   !   Subroutine that writes the heights of an [[OBJ_SURF]] object into a 2D array
   ! @endnote
   !-----------------------------------------------------------------------------------------
   subroutine trans_surf_tab(surf, tab)
   implicit none
   type(OBJ_SURF), intent(inout) :: surf                            !! *object ```OBJ_SURF``` that contains the heights*
   real(kind=R8), dimension(:, :), allocatable, intent(out) :: tab  !! *height array*
      integer(kind=I4) :: long, larg, i, j, k
      real(kind=R8)    :: unit_z

      long = surf%xres
      larg = surf%yres
      allocate( tab(1:long, 1:larg) )

      unit_z = unit2IUc(surf%dz_unit)

      do j = 1, larg
      do i = 1, long
         k = (j-1)*long +i
         tab(i, j) = (surf%val(k)*surf%dz+ surf%Zoffset)*unit_z
      enddo
      enddo
      deallocate(surf%val)
   return
   endsubroutine trans_surf_tab


   !=========================================================================================
   !< @note
   !   Subroutine that opens a surface file ```.sur``` or ```.dat```
   !
   !   The heights are centred, scaled then put into a vector.
   ! @endnote
   ! @warning
   !    If the scale factor ```sq``` is negative, the heights are not scaled when reading ```.sur```
   ! @endwarning
   ! @warning
   !    By default, the ```.sur``` header is dumped
   ! @endwarning
   !-----------------------------------------------------------------------------------------
   subroutine read_surf(nom_fic, mu, sq, tab_s, scal)
   implicit none
   character(len=*), intent(in)  :: nom_fic                     !! *file name*
   real(kind=R8),    intent(in)  :: mu                          !! *if > 0, values are centered*
   real(kind=R8),    intent(in)  :: sq                          !! *desired height standard deviation*
   type(SCALE_SURF), intent(out) :: scal                        !! *object [[SCALE_SURF]]*
   real(kind=R8), dimension(:,:), allocatable, intent(out) :: tab_s  !! *height array*
      integer(kind=I4) :: style, i, ii, j, jj, k, nb, eof, tmp, nx, ny
      real(kind=R8)    :: sqs, mean, dz, lx, ly, lz
      type(OBJ_SURF)   :: surf
      character(len=3) :: ext
      real(kind=R8), dimension(:), allocatable :: x, y, z

      i   = len_trim(nom_fic)
      ext = lower( nom_fic(i-2:i) )

      if (ext == 'dat') style = SURF_DAT
      if (ext == 'sur') style = SURF_SUR

      select case (style)
         case (SURF_SUR)
            call open_surffile(fichier=trim(nom_fic), surf=surf, scal=scal, dump=.false.)
            call trans_surf_tab(surf=surf, tab=tab_s)

         case (SURF_DAT)
            call get_unit(k)
            open(unit=k, file=trim(nom_fic), status='old')
               nb = 0
               do
                  read(k, *, iostat=eof)
                  if ( eof/=0 ) exit
                  nb = nb +1
               enddo
            rewind(k)
               allocate( x(1:nb) )
               allocate( y(1:nb) )
               allocate( z(1:nb) )
               do i = 1, nb
                  read(k, *) x(i), y(i), z(i)
               enddo
            close(k)

            ! the triplet x, y, z is sorted according x
            !-----------------------------------------------------------
            call sort_real_2real(g=1, d=nb, rtabref=x(1:nb), rtab1=y(1:nb), rtab2=z(1:nb))
            i = 1
            do
               if ( abs(x(i) -x(1)) > 100*EPS_R8 ) exit
               i = i +1
            enddo
            scal%dx = abs(x(i) -x(1))

            j = 1
            do
               if (abs(x(j +1) -x(j))>1.0e-10) exit
               j = j +1
            enddo

            ny = j ! number of same abscissae for a given column
            if (mod(nb, ny) /=0 ) STOP 'READ_SURF, non rectangular mesh'

            nx = nb/ny

            scal%xres = nx
            scal%yres = ny

            do i = 1, nx
               ii = (i-1)*ny +1
               jj = ii +ny -1
               call sort_real_1real(g=1, d=ny, rtabref=y(ii:jj), rtab1=z(ii:jj))
            enddo

            j = 1
            do
               if ( abs(y(j) -y(1)) > 100*EPS_R8 ) exit
               j = j +1
            enddo
            scal%dy = abs(y(j) -y(1))

            allocate( tab_s(1:nx, 1:ny) )
            k = 0
            do i = 1, nx
               do j = 1, ny
                  k = k +1
                  tab_s(i, j) = z(k)
               enddo
            enddo

            lx = maxval(x)     -minval(x)
            ly = maxval(y)     -minval(y)
            lz = maxval(tab_s) -minval(tab_s)

            scal%dz = lz/(nx*ny)

            scal%xlength_unit = 'm ' ; scal%dx_unit = 'm '
            scal%ylength_unit = 'm ' ; scal%dy_unit = 'm '
            scal%zlength_unit = 'm ' ; scal%dz_unit = 'm '
            !call sort_real(g=1, d=nx*ny, rtabref=z(1:nx*ny))
            !scal%dz = 1.e+10
            !do i = 2, nx*ny
            !   dz = z(i) -z(i-1)
            !   if (abs(dz) < 100*EPS_R8) cycle
            !   scal%dz = min(scal%dz, dz)
            !enddo
            !scal%dz = dz

            deallocate(x, y, z)
      endselect

      nx = scal%xres
      ny = scal%yres

      ! centering ?
      if (mu > 0.) then
         mean    = sum(tab_s(1:nx, 1:ny)) / (nx * ny)
         scal%mu = mean
         tab_s(1:nx, 1:ny) = tab_s(1:nx, 1:ny) -mean
      endif

      ! scaling ?
      if (sq > 0.) then
         sqs     = (sum(tab_s(1:nx, 1:ny) ** 2) / (nx * ny)) ** (0.5_R8)
         scal%si = sqs
         tab_s(1:nx, 1:ny) = tab_s(1:nx, 1:ny) * sq / sqs
      endif

   return
   endsubroutine read_surf


   !=========================================================================================
   !< @note
   !   Subroutine that creates an object [[OBJ_SURF]]
   ! @endnote
   !-----------------------------------------------------------------------------------------
   subroutine build_surf(surf, tab)
   implicit none
   type(OBJ_SURF), intent(inout) :: surf  !! *resulting object ```OBJ_SURF```*
   real(kind=R8), dimension(1:surf%xres, 1:surf%yres), intent(in) :: tab
      integer(kind=I4)   :: i, j, k, nx, ny
      real(kind=R8)      :: max_n, min_t, max_t, mil_t, amp_t, unit_x, unit_y, unit_z

      nx = surf%xres
      ny = surf%yres

      surf%nofpoints = nx*ny

      unit_x = unit2IUc(surf%dx_unit)
      unit_y = unit2IUc(surf%dy_unit)
      unit_z = unit2IUc(surf%dz_unit)

      if (allocated(surf%val)) deallocate(surf%val)
      allocate(surf%val(1:surf%nofpoints))

      min_t = minval( tab(1:nx, 1:ny) )/unit_z
      max_t = maxval( tab(1:nx, 1:ny) )/unit_z

      mil_t = 0.5_R8*(min_t +max_t) ! middle of the range
      amp_t = max_t -min_t          ! range amplitude

      surf%ZOffset = mil_t

      max_n = 0.5*huge(1)     ! the heights are integers, allowed to span
                              ! half the positive integer range.
      surf%dz = amp_t/max_n   ! subsequent dz

      k = 0
      do j = 1, ny
      do i = 1, nx
         k = k +1
         surf%val(k) = nint( (tab(i, j)/unit_z -surf%ZOffset)/surf%dz ) !
      enddo
      enddo
      surf%zmin = minval(surf%val)
      surf%zmax = maxval(surf%val)

   return
   endsubroutine build_surf


   !=========================================================================================
   !< @note
   !   Subroutine that writes an object [[OBJ_SURF]] in a file
   ! @endnote
   !-----------------------------------------------------------------------------------------
   subroutine write_surffile(fichier, surf)
   implicit none
   character(len=*), intent(in)  :: fichier  !! *file to be written*
   type(OBJ_SURF), intent(inout) :: surf     !! *object ```OBJ_SURF``` to write*
      integer(kind=I4) :: i, k

      call get_unit(k)
      open(k,  file=trim(fichier),     & !
               form='unformatted',     & !
               access="stream",        & ! beware the "frecord-marker" in other modes
               action="write",         & !
               position="rewind",      & !
               status="replace",       & !
               convert='little_endian')

         write(k) surf%signature, surf%format, surf%nobjects, surf%version, surf%type, surf%object_name,                &
                  surf%operator_name, surf%material_code, surf%acquisition, surf%range, surf%special_points,            &
                  surf%absolute, surf%reserved, surf%pointsize, surf%zmin, surf%zmax, surf%xres, surf%yres,             &
                  surf%nofpoints, surf%dx, surf%dy, surf%dz, surf%xaxis, surf%yaxis, surf%zaxis, surf%dx_unit,          &
                  surf%dy_unit, surf%dz_unit, surf%xlength_unit, surf%ylength_unit, surf%zlength_unit,                  &
                  surf%xunit_ratio, surf%yunit_ratio, surf%zunit_ratio, surf%imprint, surf%inversion, surf%leveling,    &
                  surf%obsolete, surf%seconds, surf%minutes, surf%hours, surf%day, surf%month, surf%year, surf%dayof,   &
                  surf%measurement_duration, surf%obsolete2, surf%comment_size, surf%private_size, surf%client_zone,    &
                  surf%XOffset, surf%YOffset, surf%ZOffset, surf%reservedzone,                                          &
                  (surf%val(i), i=1, surf%nofpoints)
      close(k)

   return
   endsubroutine write_surffile


   !=========================================================================================
   !< @note
   !   Subroutine that writes a height array into a surface file ```.sur``` or ```.dat```
   ! @endnote
   !-----------------------------------------------------------------------------------------
   subroutine write_surf(nom_fic, tab_s, scal)
   implicit none
   character(len=*), intent(in)    :: nom_fic   !! *file name*
   type(SCALE_SURF), intent(inout) :: scal      !! *object [[SCALE_SURF]]*
   real(kind=R8), dimension(1:scal%xres, 1:scal%yres), intent(in) :: tab_s
      character(len=3) :: ext
      integer(kind=I4) :: style, i, j, k
      type(OBJ_SURF)   :: surf_s
      real(kind=R8)    :: dx, dy

      i   = len_trim(nom_fic)
      ext = lower( nom_fic(i-2:i) )

      if (ext == 'dat') style = SURF_DAT
      if (ext == 'sur') style = SURF_SUR

      select case (style)
         case (SURF_SUR)
            call scal2surf(scal, surf_s)
            call build_surf(surf=surf_s, tab=tab_s(1:scal%xres, 1:scal%yres))
            surf_s%comment_size = 0 ! to increase compatibility with mountains
            call write_surffile(fichier=trim(nom_fic), surf=surf_s)
            call surf2scal(surf_s, scal)

         case (SURF_DAT)
            dx = scal%dx
            dy = scal%dy
            call get_unit(k)
            open(k, file=trim(nom_fic))
               do i = 1, scal%xres
                  do j = 1, scal%yres
                     write(k,*) (i-1)*dx, (j-1)*dy, tab_s(i, j)
                  enddo
               enddo
            close(k)
      endselect
   return

   endsubroutine write_surf


endmodule surfile

!============ EN TETE TYPIQUE
!~  "DIGITAL SURF"
!~  0
!~  1
!~  1
!~  2
!~  "* 16-774-lm1-pnm bouches *"
!~  ""
!~  0
!~  0
!~  0
!~  0
!~  1
!~  ""
!~  32
!~  -19122091
!~  4341882
!~  512
!~  512
!~  262144
!~  3.55476077E-04
!~  3.54291551E-04
!~  3.25564076E-10
!~  "X"
!~  "Y"
!~  "Z"
!~  "mm"
!~  "mm"
!~  "mm"
!~  "mm"
!~  "mm"
!~  "mm"
!~  1.00000000
!~  1.00000000
!~  1.00000000
!~  0
!~  0
!~  0
!~  ""
!~  0
!~  0
!~  0
!~  0
!~  0
!~  0
!~  0
!~  0.00000000
!~  ""
!~  0
!~  0
!~  ""
!~   0.00000000
!~   0.00000000
!~  -0.00000000
!~  ""
!~

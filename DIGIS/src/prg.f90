
!< author: Arthur Francisco
!  version: 1.0.0
!  date: April, 16 2019
!
!  <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!     **Routines to handle Digital Surf binary format (.sur). Example of use**
!  </span>

program test_surfile
use data_arch
use surfile
implicit none

real(kind=R8), dimension(:,:), allocatable :: tab  !! *height array*
type(SCALE_SURF) :: scal                           !! *object [[SCALE_SURF]], fortran type version*
type(OBJ_SURF)   :: surf                           !! *object [[OBJ_SURF]],         c type version*

integer(kind=I4) :: i, j

   ! --- reads a xyz ascii file, writes the corresponding "sur" file and dumps the header
   call init_scal(scal)                          ! creates an empty surface type (fortran form)

   call read_surf(nom_fic = "sur/600x300.dat", & !  in; three columns in ascii format : x y f(x,y); no header; tab separation
                       mu = -1._R8,            & !  in; no centering required
                       sq = -1._R8,            & !  in; no normalization required
                    tab_s = tab,               & ! out; array containing the surface
                     scal = scal)                ! out; surface type containing some informations like length, width, etc.

   call write_surf(nom_fic = "out/600x300_dat_to_not-scaled.sur",  & !    in; filename of the ".sur" to be created
                     tab_s = tab,                                  & !    in; surface array
                      scal = scal)                                   ! inout; surface type
   call scal2surf(scal, surf)                                        !        surface type transformation: fortran form to c form
   call trans_surf_txt(surf = surf,                                & !    in; dumps the surface header ...
                    fichier = "out/600x300_dat_to_not-scaled.txt", & !    in; ... in a specified file ...
                        xyz = .false.)                               !    in; ... and no f(x,y) dump.

   deallocate(tab)

   ! --- reads a "sur" file, writes its scaled form in a "sur" file and a xyz file
   call read_surf(nom_fic = "sur/600x300.sur",  & !  in; Digital Surf format
                       mu = 1._R8,              & !  in; data will be centered
                       sq = 1._R8,              & !  in; data will be normalized (with the standard deviation)
                    tab_s = tab,                & ! out; array containing the surface
                     scal = scal)                 ! out; surface type containing some informations like length, width, etc.

   call write_surf(nom_fic = "out/600x300_resu_scaled.sur",       & !    in; filename of the ".sur" to be created
                     tab_s = tab,                                 & !    in; surface array
                      scal = scal)                                  ! inout; surface type

   call write_surf(nom_fic = "out/600x300_resu_scaled.dat",       & !    in; filename of the ascii ".dat" to be created
                     tab_s = tab,                                 & !    in; surface array
                      scal = scal)                                  ! inout; surface type

   deallocate(tab)

   ! --- creates a surface and writes a ".sur" file
   allocate(tab(600,300))
   do j = 1, 300
      do i = 1, 600
         tab(i, j) = 1.e+9*cos(6 *2*PI_R8*i/600)*cos(3 *2*PI_R8*j/300) +1.e+8
      enddo
   enddo

   call init_scal(scal = scal,      &           ! out; creates a surface type, containing ...
                    nx = 600,       &           !  in; ... the number of points along x ...
                    ny = 300,       &           !  in; ... the number of points along y ...
                    lx = 1.0e-3_R8, &           !  in; ... the length (default unit : m) ...
                    ly = 0.5e-3_R8, &           !  in; ... the width ...
                unit_z = 'Pa'       )           !  in; ... and the unit along z.

   call write_surf(nom_fic = "out/cos.sur",  &  !    in; filename of the ".sur" to be created
                     tab_s = tab,            &  !    in; surface array
                      scal = scal)              ! inout; surface type

   call read_surf(nom_fic = "out/cos.sur", &    !  in; Digital Surf format
                       mu = -1._R8,        &    !  in; data will not be centered
                       sq = -1._R8,        &    !  in; data will not be normalized (with the standard deviation)
                    tab_s = tab,           &    ! out; array containing the surface
                     scal = scal)               ! out; surface type containing some informations like length, width, etc.
   call scal2surf(scal, surf)                   !      surface type transformation: fortran form to c form
   call trans_surf_txt(surf = surf,          &  ! in; dumps the surface header ...
                    fichier = "out/cos.txt", &  ! in; ... in a specified file ...
                        xyz = .false.)          ! in; ... and no f(x,y) dump.

   deallocate( tab )

stop
endprogram test_surfile


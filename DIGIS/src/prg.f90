program test_surfile
use data_arch
use surfile
implicit none

real(kind=R8), dimension(:,:), allocatable :: tab  !! *height array*
type(SCALE_SURF) :: scal                           !! *object [[SCALE_SURF]]*
type(OBJ_SURF)   :: surf
integer(kind=I4) :: i, j

   ! --- reads a xyz ascii file, writes the corresponding "sur" file and dumps the header
   call init_scal(scal)

   call read_surf(nom_fic = "sur/600x300.dat", & !
                       mu = -1._R8,            & !
                       sq = -1._R8,            & !
                    tab_s = tab,               & !
                     scal = scal)

   call write_surf(nom_fic = "out/600x300_dat_to_not-scaled.sur", & !
                     tab_s = tab,                                 & !
                      scal = scal)
   call scal2surf(scal, surf)
   call trans_surf_txt(surf = surf,                                & !
                    fichier = "out/600x300_dat_to_not-scaled.txt", & !
                        xyz = .false.)

   deallocate(tab)

   ! --- reads a "sur" file, writes its scaled form in a "sur" file and a xyz file
   call read_surf(nom_fic = "sur/600x300.sur", & !
                       mu = 1._R8,              & !
                       sq = 1._R8,              & !
                    tab_s = tab,                & !
                     scal = scal)

   call write_surf(nom_fic = "out/600x300_resu_scaled.sur",       & !
                     tab_s = tab,                                 & !
                      scal = scal)

   call write_surf(nom_fic = "out/600x300_resu_scaled.dat",       & !
                     tab_s = tab,                                 & !
                      scal = scal)

   deallocate(tab)

   ! --- creates a surface and writes a "sur" file
   allocate(tab(600,300))
   do j = 1, 300
      do i = 1, 600
         tab(i, j) = 1.e+9*cos(6 *2*PI_R8*i/600)*cos(3 *2*PI_R8*j/300) +1.e+8
      enddo
   enddo

   call init_scal(scal = scal,      & !
                    nx = 600,       & !
                    ny = 300,       & !
                    lx = 1.0e-3_R8, & ! default unit : m
                    ly = 0.5e-3_R8, & !
                unit_z = 'Pa'       ) !

   call write_surf(nom_fic = "out/cos.sur",  & !
                     tab_s = tab,            & !
                      scal = scal)

   call read_surf(nom_fic = "out/cos.sur", & !
                       mu = -1._R8,        & !
                       sq = -1._R8,        & !
                    tab_s = tab,           & !
                     scal = scal)
   call scal2surf(scal, surf)
   call trans_surf_txt(surf = surf,          & !
                    fichier = "out/cos.txt", & !
                        xyz = .false.)

   deallocate( tab )

stop
endprogram test_surfile


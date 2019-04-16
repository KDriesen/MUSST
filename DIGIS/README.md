DIGIS
=====

DIGIS is a modern fortran program that handles the [Digital Surf](https://www.digitalsurf.com/) format ```.sur```. The ```.sur``` format has been implemented thanks [Gwyddion](http://gwyddion.net/) source file ```surffile.c```.

The ```.sur``` is a little endian binary format which consists of a 512 bytes header, followed by the heights (32 bit signed) 

The two principal subroutines are ```read_surf``` and ```write_surf``` which 
* reads a ```.sur``` file and stores the data in an array
* writes an array in a ```.sur``` file.

Table of Contents
-----------------
- [Installation](#Installation)
- [Execution](#Execution)
	- [Example 1](#Example-1)
	- [Example 2](#Example-2)
	- [Example 3](#Example-3)
- [License](#License)

Installation
============

* simple compilation : ```$ make```
* debug compilation : ```$ make debug```

Execution
=========

To run the main program : ```$./prg```

The main program ```prg.f90``` provides 3 simple examples of use.

Example 1
---------

```fortran
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
                      
   call scal2surf(scal, surf)                                        ! surface type transformation: fortran form to c form
   
   call trans_surf_txt(surf = surf,                                & ! in; dumps the surface header ...
                    fichier = "out/600x300_dat_to_not-scaled.txt", & ! in; ... in a specified file ...
                        xyz = .false.)                               ! in; ... and no f(x,y) dump.

   deallocate(tab)
```
[top](#table-of-contents)

Example 2
---------

```fortran
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
```
[top](#table-of-contents)

Example 3
---------

```fortran
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
                     
   call scal2surf(scal, surf)                   ! surface type transformation: fortran form to c form
   
   call trans_surf_txt(surf = surf,          &  ! in; dumps the surface header ...
                    fichier = "out/cos.txt", &  ! in; ... in a specified file ...
                        xyz = .false.)          ! in; ... and no f(x,y) dump.

   deallocate( tab )
```
[top](#table-of-contents)

License
=======

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but **without any warrenty**; without even the implied warranty of **merchantability** or **fitness for a particular purpose**. 
See the GNU General Public License for more details.

[top](#table-of-contents)


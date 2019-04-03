
!! author: Noël Brunetière<br/>&emsp;Arthur Francisco
!! version: 1.0.0
!! date: April,17 2017
!! summary: main program to run MUSST on a configuration file

!< <span style="color: #337ab7; font-family: cabin; font-size: 1.5em;">
!<       **Subroutine to run a test with MUSST**
!< </span>

program main
use omp_lib
use test_musst, only : run_test
use data_arch,  only : I4, NB_THREADS_MAX, OUT_U=>OPU, get_unit
use num_param,  only : VERBOSE, OUTPUT_FILE, OPU
use solver,     only : SOLVER_BS, SOLVER_TS, SOLV_MESS
implicit none

character(len=128) :: prg_arg
character(len=128) :: job_file

character(len= 8) :: date
character(len=10) :: time
character(len=15) :: repos

integer(kind=I4) :: ARCHIVE

! a result directory is created under "/out" with the date as name
! ----------------------------------------------------------------
call date_and_time(date, time)
repos = date//'_'//time(1:6)
call execute_command_line("mkdir out/"//repos, wait=.true.)

! the job file is the program argument
! ------------------------------------
prg_arg  = repeat(' ', len(prg_arg ))
job_file = repeat(' ', len(job_file))
call get_command_argument(1, prg_arg)     ! argument one: job file
if (len_trim(prg_arg) == 0) then          ! if there is no job file, stop
   stop 'no job file, stop'
else
   job_file = trim(adjustl(prg_arg))
endif

! the job file is copied in the result directory
! ----------------------------------------------
call execute_command_line("cp "//trim(job_file)//" out/"//repos, wait=.true.)

! the job file is processed
! -------------------------
call read_config

stop

contains

   subroutine read_config()
   implicit none
      integer(kind=I4)   :: jf, err_read, test_num
      character(len=032) :: word
      character(len=128) :: job_copy

      job_copy = "out/"//repos//job_file(4:len_trim(job_file))

      call get_unit(jf) ; open(jf, file = trim(job_file), status = 'old')
         do
            word = repeat(' ', len(word))
            read(jf, *, iostat = err_read) word

            if ( index(word, 'END_FILE'      ) /= 0 ) then
               if (OPU/=OUT_U) close(OPU)
               exit
            endif

            if ( index(word, 'VERBOSE'       ) /= 0 ) then
               read(jf, *) VERBOSE
            endif

            if ( index(word, 'SOLV_MESS'     ) /= 0 ) then
               read(jf, *) SOLV_MESS
            endif

            if ( index(word, 'OUTPUT_UNIT'   ) /= 0 ) then
               read(jf, *) OPU
               if (OPU==0) then
                  OPU = OUT_U
               else
                  read(jf, *) OUTPUT_FILE
                  call get_unit(OPU)
                  open(unit = OPU,                          &
                       file = trim(adjustl(OUTPUT_FILE)),   &
                     status = 'unknown')
               endif
            endif

            if ( index(word, 'SOLVER_BS'     ) /= 0 ) then
               read(jf, *) SOLVER_BS
               if (SOLVER_BS==1) stop "some remaining problems with SuLU, choose '3' instead"
               if (SOLVER_BS==2) stop "MUMPS is multithreaded, it is designed for TS, choose '3' instead"
            endif

            if ( index(word, 'SOLVER_TS'     ) /= 0 ) then
               read(jf, *) SOLVER_TS
            endif

            if ( index(word, 'EXEC_MUSST'    ) /= 0 ) then
               read(jf, *) test_num
               call run_test(test_num, jf, repos)
            endif

            if ( index(word, 'NB_THREADS_MAX') /= 0 ) then
               read(jf, *) NB_THREADS_MAX
               if (NB_THREADS_MAX<0) then
                  !$ NB_THREADS_MAX = omp_get_num_procs()
               endif
            endif

            if ( index(word, 'ARCHIVE') /= 0 ) then
               read(jf, *) ARCHIVE
               if (ARCHIVE==1) then
                  call execute_command_line("sh bin/save/save.sh", wait=.true.)
                  call execute_command_line("mv *.7z out/"//repos//"/")
               endif
            endif

         enddo
      close(jf)
   return
   endsubroutine read_config

endprogram main

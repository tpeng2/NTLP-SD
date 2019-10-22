      subroutine close_his
c
c ---- close history files
c
      use pars
      logical there
c
c ---- root closes and checks the files
c
      close(nhis1)
      close(nhisp)
      close(nrad)
      inquire(file=path_sav_h,exist=there)
      if(.not.there) then
         write(6,8000) path_sav_h
         call mpi_finalize(ierr)
         stop
      endif
      inquire(file=path_sav_hp,exist=there)
      if(.not.there) then
         write(6,8100) path_sav_hp
         call mpi_finalize(ierr)
         stop
      endif
      inquire(file=path_sav_hist,exist=there)
      if(.not.there) then
         write(6,8200) path_sav_hist
         call mpi_finalize(ierr)
         stop
      endif
      write(6,7000) path_sav_h
      write(6,7100) path_sav_hp
      write(6,7200) path_sav_hist
c
      return
c -------------------- process write errors
 7000 format(' HISTORY DATA IS WRITTEN IN FILE  ',a80)
 7100 format(' PROFILE HISTORY DATA IS WRITTEN IN FILE  ',a80)
 7200 format(' RADHIST DATA IS WRITTEN IN FILE  ',a80)
 8000 format(' SR. SAVE_HIS: Trouble history file not in path =',a80)
 8100 format(' SR. SAVE_HIS: Trouble profile history file',
     +       ' not in path =',a80)
 8200 format(' SR. SAVE_RADPDF: Trouble rad histogram file =',a80)
      end

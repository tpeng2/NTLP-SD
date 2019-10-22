      subroutine read_input_file
      use pars
      use particles
      use con_data
      implicit none

      character(48) :: label
      character(180) :: params_dir
      namelist /step_params/ iti,itmax,imean,ihst,itape,
     +                       itstr,it_his,it_viz,i_viz,itn

      namelist /grid_params/ ncpu_s, Uo, Ttop, Tbot,
     +         qstar, tsfcc, ugcont, vgcont,
     +         zi, zl, xl, yl, zw1, dpdx,
     +         cfl,dt_new,surf_RH,psurf

      namelist /path_names/ path_seed,path_part,path_res,
     +         path_sav,path_his,path_viz_xy,path_viz_xz,
     +         path_viz_yz,path_stuf,path_ran,path_histog

      namelist /flags/ ismlt,ifree,isfc,iradup,
     +         iupwnd,ibuoy,ifilt,itcut,isubs,ibrcl,iocean,
     +         method,idebug,iz_space,ivis0,ifix_dt,new_vis,iDNS,
     +         icouple,iTcouple,iHcouple,ievap,ifields,ilin,
     +         ineighbor,icoalesce,ipart_method,
     +         ipartdiff,isfs,iexner


      namelist /constants/ rhoa, nuf, Cpa, Pra, Sc,
     +         tnumpart,mult_init,rhow,part_grav,
     +         Cpp,Mw,Ms,Ru,Gam,Ion,Os,Sal,Rd,
     +         radius_init,
     +         grav, t00,fcor,zo,
     +         vp_init,Tp_init,qf_init


      !params.in contains namelists to read
      !open(12, file="./params.in", status="old")
      call get_command_argument(1,params_dir)
      open(12,file=params_dir,status="old")

      read(12,nml=step_params)
      if (myid==0) print step_params

      read(12,nml=flags)
      if (myid==0) print flags

      read(12,nml=grid_params)
      if (myid==0) print grid_params

      read(12,nml=path_names)
      if (myid==0) print path_names

      read(12,nml=constants)
      if (myid==0) print constants
      CpaCpp = Cpa/Cpp



      end subroutine read_input_file

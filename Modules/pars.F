c     ----------------------------------------------------------------------
      module pars
c     ----------------------------------------------------------------------
      ! Define constant
      real :: pi
      ! Time step parameters
      integer :: iti,itmax,imean,ihst,itape
      integer :: itstr,it_his,it_viz,i_viz
      integer :: itcplbegin !time when 2-way coupling begins
      integer :: itclsbegin !time when particle coalesce begins

      integer :: ismlt,ifree,iradup,
     +     iupwnd,ibuoy,ifilt,itcut,isubs,ibrcl,iocean,method,
     +     idebug,iz_space,ivis0,ifix_dt,new_vis,iDNS,
     +     icouple,iTcouple,iHcouple,ievap,ifields,isurfUlog,
     +     itempfld,ilin,ineighbor,
     +     icoalesce,ipart_method,ipartdiff,isfs,iexner
c     
      integer :: iparsdefault ! for default setting
      integer :: ifixnumpart ! for fixed particle numbers
      !For particle generation
      integer :: iSSGF,ivpran,irpsmp,iFnumsmp,iszcls,impscl,impmnl
c     -------- j_recl = record length in "bytes" for history file
c     k8     = kind parameter for integers in mpi_io routines
c     
      integer, parameter :: j_recl=4, k8=8

      integer, parameter :: nscl = 2, nvar = (4+nscl)
      integer, parameter :: maxnx=128,maxny=128,maxnz=128
      integer, parameter :: maxnz1=maxnz+1,
     + maxnz2=maxnz+2,maxnx2=maxnx+2,maxny2=maxny+2

      !Boundary conditions and other parameters
      real :: Uo, Ttop(nscl), Tbot(nscl),qstar(nscl),tsfcc(nscl)
      real :: grav,t00,fcor,zi,zl,xl,yl,zw1,dpdx
      real :: ugcont,vgcont,zo,surf_RH,psurf
      integer :: isfc(nscl)

      !Material parameters, both air and particles
      real :: Lv = (25.0 - 0.02274*26.0)*10.0**5
      real :: Sal,Ms,Gam,Ion,Os
      real :: rhow,nuf 
      real :: CpaCpp,Cpa,Cpp,Pra,Sc,Mw,Ru,rhoa,Rd
      real :: m_s, m_p
      
      !Number-distribution of particle sizes
      integer :: Nrp
      ! SSGF for subroutine create_random_radius
      real,allocatable :: Sn(:)
      real,allocatable :: Snscl(:),rpcrop(:),CDF(:)
      

c     
c     ------------ leave following definitions as is
c     
c     ----------------------------------------------------------------------
      integer ::    nnx, nny, nnz, nxy, ncx, nnxp1, nnyp1, ncy,
     +     nnxp2, nnyp2, nnzp1, ivis, nnzm1, isize, krec,
     +     izs, ize, ixs, ixe, jxs, jxe, kxs, kxe,
     +     mxs, mxe, iss, ise, iys, iye, jys, jye
c     ----------------------------------------------------------------------
      character case*3
c     ----------------------------------------------------------------------
      integer  ::   nvel, npre, nhis1, nprt,nrad,
     +     nhisp, nvelc, 
     +     nviz_z, nviz_y, 
     +     nviz_x, nviz_s,
     +     kfile, jfile, ibcl, ibcu,
     +     igrdr, ichannel,  imach, itn, it_his_nxt, it_viz_nxt
      logical ::    mnout, micut, mtape, mhis, msave,mtrans,
     +     l_root, l_debug, msave_v, mviz
c     ----------------------------------------------------------------------
      real    ::    windm,u1xy,v1xy,t1xy(nscl),
     +     t10xy(nscl),au13m,au23m,aut3m(nscl),tsfcm(nscl),
     +     thstar(nscl), eavg(maxnz), tr_tau(0:maxnz),
     +     pxym(0:maxnz1), zi_min
      integer ::    izi, iz_min
      real, allocatable ::
     +     wind(:,:), tau13m(:,:), tau23m(:,:), 
     +     taut3m(:,:,:), t_grnd(:,:,:)
c     ----------------------------------------------------------------------
      real ::       u_mn(0:maxnz1), v_mn(0:maxnz1),
     +     w_mn(0:maxnz1), t_mn(0:maxnz1,nscl), alphaC(0:maxnz1,nscl)
c     ----------------------------------------------------------------------
      real ::       dzw(0:maxnz2), dzu(0:maxnz2),
     +     dzw_i(0:maxnz2), dzu_i(0:maxnz2)
c     ----------------------------------------------------------------------
      real ::       t_factor, t_ref, c_rate, t_surf_i
c     ----------------------------------------------------------------------
      real ::       dfac(maxnz), dsl_z(0:maxnz1),
     +     xksurf, viscon, vise, almin_c,stabmin,
     +     ck,ceps,csmag,stab_c,vis_mean(0:maxnz)
      integer ::    nmatch
c     ----------------------------------------------------------------------
      real ::       zetas(3), gama(3), etas(4), dt_new,
     +     umax,vmax,wmax, wabs, vismax,
     +     cfl, tzero,
     +     ucfl, vcfl, wcfl
c     ----------------------------------------------------------------------
      character*80  path_res, path_sav, path_his, path_prt,
     +     path_hp, path_sav_hp, path_part, path_ran,
     +     path_v, path_c, path_p, path_h,
     +     path_sav_v, path_sav_c,
     +     path_sav_p, path_sav_h, path_sav_part,
     +     path_histog, path_sav_hist,
     +     bad_news
      character case_inp*3
      character*80 path_viz_xy, path_viz_xz, path_viz_yz, path_stuf,
     +             path_seed
c     ----------------------------------------------------------------------
      integer ::    myid, numprocs, i_root,
     +     ziloc, myid_newvis, ncpu_s, ncpu_z, maxp
      integer, allocatable, dimension(:) :: 
     +     ix_s, ix_e, jx_s, jx_e,
     +     kx_s, kx_e, mx_s, mx_e,
     +     iy_s, iy_e, jy_s, jy_e,
     +     is_s, is_e, iz_s, iz_e
      end module pars
c     ----------------------------------------------------------------------
      module particles
      integer :: rproc,trproc,tproc,tlproc,lproc,blproc,bproc,brproc
      integer :: pr_r,pl_r,pt_r,pb_r,ptr_r,ptl_r,pbl_r,pbr_r
      integer :: pr_s,pl_s,pt_s,pb_s,ptr_s,ptl_s,pbl_s,pbr_s
      real :: ymin,ymax,zmin,zmax,xmax,xmin
      real, allocatable :: uext(:,:,:), vext(:,:,:), wext(:,:,:)
      real, allocatable :: u_t(:,:,:), v_t(:,:,:), w_t(:,:,:)
      real, allocatable :: Text(:,:,:),T_t(:,:,:)
      real, allocatable :: T2ext(:,:,:),T2_t(:,:,:)
      real, allocatable :: partTsrc(:,:,:),partTsrc_t(:,:,:)
      real, allocatable :: partHsrc(:,:,:),partHsrc_t(:,:,:) 
      real, allocatable :: partTEsrc(:,:,:),partTEsrc_t(:,:,:) 
      real, allocatable :: partcount_t(:,:,:),partsrc_t(:,:,:,:)
      real, allocatable :: vpsum_t(:,:,:,:),vpsqrsum_t(:,:,:,:)
      real, allocatable :: upwp_t(:,:,:),upwp(:,:,:)
      real, allocatable :: partcount(:,:,:),partsrc(:,:,:,:)
      real, allocatable :: vpsum(:,:,:,:),vpsqrsum(:,:,:,:)
      real, allocatable :: Tpsum(:,:,:),Tpsum_t(:,:,:)
      real, allocatable :: Tpsqrsum(:,:,:),Tpsqrsum_t(:,:,:)
      real, allocatable :: Tfsum(:,:,:),Tfsum_t(:,:,:)
      real, allocatable :: qfsum(:,:,:),qfsum_t(:,:,:)
      real, allocatable :: wpTpsum(:,:,:),wpTpsum_t(:,:,:)
      real, allocatable :: radsum(:,:,:),radsum_t(:,:,:) 
      real, allocatable :: rad2sum(:,:,:),rad2sum_t(:,:,:) 
      real, allocatable :: mpsum(:,:,:),mpsum_t(:,:,:) 
      real, allocatable :: mwsum(:,:,:),mwsum_t(:,:,:)
      real, allocatable :: qstarsum(:,:,:),qstarsum_t(:,:,:) 

c     --- SFS velocity calculation ---------
      real, allocatable :: sigm_s(:,:,:),sigm_sdx(:,:,:),sigm_sdy(:,:,:)
      real, allocatable :: sigm_sdz(:,:,:),sigm_sext(:,:,:)
      real, allocatable :: sigm_sdxext(:,:,:),sigm_sdyext(:,:,:)
      real, allocatable :: sigm_sdzext(:,:,:)
      real, allocatable :: vis_ss(:,:,:),vis_sext(:,:,:) 
!     --- Count reintroduced particle
      integer :: num_reintro,tot_reintro,his_reintro !Tianze 9/11/15

      integer :: particletype,pad_diff
      integer :: numpart,tnumpart,ngidx,massfrac
!     --- 10-m wind speed for SSGF only
      integer :: U10
!     --- manual size input, defined in parameter files
!     --- for 2-3 sizes
      integer :: phmt1,phmt2,phmt3
!     --- for 6 sizes
      integer :: phm1, phm2,phm3,phm4,phm5,phm6
      integer :: iseed
      integer :: num100=0, num1000=0, numimpos=0
      integer :: tnum100, tnum1000, tnumimpos
      integer :: denum, actnum, tdenum, tactnum
      integer :: num_destroy=0,tnum_destroy=0

      real :: Rep_avg,part_grav(3)
      real :: radavg,radmin,radmax,tempmin,tempmax,qmin,qmax
      real :: vp_init(3),Tp_init,qf_init,radius_init
      integer*8 :: mult_init

      real :: avgres=0,tavgres=0

      integer, parameter :: histbins = 512
      real :: hist_rad(histbins+2)
      real :: bins_rad(histbins+2)      


      !for smallest rp incremental
      real :: rp_inc 
      !Following three lines are all set in the input file 9/18/18
      real,allocatable :: rpsmp(:)  !candidate size array
      integer,allocatable::irpind(:) !size index
      real :: rpmax,rpmed,rpmin,Vpmax
      
      
      !REMEMBER: IF ADDING ANYTHING, MUST UPDATE MPI DATATYPE!
      type :: particle
      integer :: pidx,procidx,nbr_pidx,nbr_procidx
      real :: vp(3),xp(3),uf(3),xrhs(3),vrhs(3),Tp,Tprhs_s
      real :: Tprhs_L,Tf,radius,radrhs,qinf,qstar,dist
      real :: res,m_s
      real :: u_sub(3),sigm_s
      integer*8 :: mult
      type(particle), pointer :: prev,next
      end type particle

      type(particle), pointer :: part,first_particle
      end module particles
c     --------------------------------------------------------------------- 
      module fields
      real, allocatable :: 
     +     u(:,:,:), v(:,:,:), w(:,:,:), t(:,:,:,:), e(:,:,:), 
     +     r1(:,:,:), r2(:,:,:), r3(:,:,:), r4(:,:,:,:), r5(:,:,:)
      real, allocatable :: 
     +     ux(:,:,:), uy(:,:,:), vx(:,:,:), vy(:,:,:), 
     +     wx(:,:,:), wy(:,:,:),
     +     p(:,:,:), ptop(:,:,:), vis_m(:,:,:), vis_s(:,:,:,:) 
      real, allocatable :: 
     +     ubc(:,:,:), vbc(:,:,:), wbc(:,:,:), tbc(:,:,:,:), 
     +     ebc(:,:,:), pbc(:,:,:), pbc2(:,:,:)
      end module fields
c     ----------------------------------------------------------------------
      module fftwk
      real, allocatable :: trigx(:,:), trigc(:)
      end module fftwk
c     ----------------------------------------------------------------------
      module con_data
c     ----------------------------------------------------------------------
      use pars, only : nscl
      type con_d
      sequence
      real ::  vk, vkin, vk74, vk74in, 
     +     gcp, fcor_h, pi2,
     +     batagk, batag,
     +     cdbtm, dtdzf(nscl), dtjump, ugal, divgls,
     +     z1, utausv, dx, dy, dz, dt, 
     +     fnxy, dzdz, dsl, c23, dtgama, dtzeta, xkmax,
     +     time, zody, zody74,
     +     utau, wtsfc(nscl),
     +     uwsfc, vwsfc, amonin,
     +     zol, hol, smal_e, sml_eg
      end type con_d
      type(con_d), target :: c_c
      real, pointer ::
     +     vk, vkin, vk74, vk74in, 
     +     gcp, fcor_h, pi2,
     +     batagk, batag,
     +     cdbtm, dtdzf(:), dtjump, ugal, divgls,
     +     z1, utausv, dx, dy, dz, dt,
     +     fnxy, dzdz, dsl, c23, dtgama, dtzeta, xkmax, 
     +     time, zody, zody74,
     +     utau, wtsfc(:), 
     +     uwsfc, vwsfc, amonin,
     +     zol, hol, smal_e, sml_eg
      contains
      subroutine fill_cc
c     
c     --------------- pointer associations for constant variables
c     
      vk     => c_c%vk 
      vkin   => c_c%vkin
      vk74   => c_c%vk74
      vk74in => c_c%vk74in
      gcp    => c_c%gcp
      fcor_h => c_c%fcor_h
      pi2    => c_c%pi2
      batagk => c_c%batagk
      batag  => c_c%batag
      cdbtm  => c_c%cdbtm
      dtdzf  => c_c%dtdzf
      dtjump => c_c%dtjump
      ugal   => c_c%ugal
      divgls => c_c%divgls
      z1     => c_c%z1
      utausv => c_c%utausv
      dx     => c_c%dx
      dy     => c_c%dy
      dz     => c_c%dz
      dt     => c_c%dt
      fnxy   => c_c%fnxy
      dzdz   => c_c%dzdz
      dsl    => c_c%dsl
      c23    => c_c%c23
      dtgama => c_c%dtgama
      dtzeta => c_c%dtzeta
      xkmax  => c_c%xkmax
      time   => c_c%time
      zody   => c_c%zody
      zody74 => c_c%zody74
      utau   => c_c%utau
      wtsfc  => c_c%wtsfc
      uwsfc  => c_c%uwsfc
      vwsfc  => c_c%vwsfc
      amonin => c_c%amonin
      zol    => c_c%zol
      hol    => c_c%hol
      smal_e => c_c%smal_e
      sml_eg => c_c%sml_eg
      return
      end subroutine fill_cc
      end module con_data
c ----------------------------------------------------------------------
      module con_stats
        use pars
        type con_s
        sequence
        real ::  wwsb(maxnz),engz(0:maxnz1),
     +           engsbz(0:maxnz1),
     +           englez(maxnz),uxym(0:maxnz1),
     +           vxym(0:maxnz1),wxym(0:maxnz1),
     +           txym(0:maxnz1,nscl),divz(0:maxnz1),
     +           utle(maxnz,nscl), utsb(maxnz,nscl),
     +           vtle(maxnz,nscl), vtsb(maxnz,nscl),
     +           wtle(maxnz,nscl), wtsb(0:maxnz,nscl),
     +           wt_tot(maxnz,nscl),
     +           z(0:maxnz1),zz(0:maxnz1),
     +           shrz(maxnz),buyz(maxnz),
     +           triz(maxnz),
     +           uwsb(0:maxnz),vwsb(0:maxnz),
     +           uwle(maxnz),vwle(maxnz),
     +           uw_tot(maxnz),vw_tot(maxnz),
     +           wcube(maxnz), wfour(maxnz),
     +           tcube(maxnz,nscl),
     +           ups(maxnz), vps(maxnz),
     +           wps(maxnz), tps(maxnz,nscl),
     +           t_rprod(maxnz), t_wq(maxnz),
     +           t_wp(maxnz), t_tau(maxnz),
     +           t_tran(maxnz), t_buoy(maxnz),
     +           t_diss(maxnz), t_sprod(maxnz),
     +           zconc(maxnz),
     +           vp1mean(maxnz),vp2mean(maxnz),vp3mean(maxnz),
     +           vp1msqr(maxnz),vp2msqr(maxnz),vp3msqr(maxnz),
     +           m1src(maxnz),m2src(maxnz),m3src(maxnz),
     +           upwpm(maxnz),
     +           Tpmean(maxnz),Tpmsqr(maxnz),
     +           Tfmean(maxnz),qfmean(maxnz),
     +           wpTpm(maxnz),Tpsrc(maxnz),
     +           Tv_tran(maxnz,nscl),Tv_prod(maxnz,nscl),
     +           Tv_diss(maxnz,nscl),
     +           Tv_part1(maxnz,nscl),Tv_part2(maxnz),
     +           radmean(maxnz),rad2mean(maxnz),
     +           mpmean(maxnz),mwmean(maxnz), 
     +           phimp(maxnz),phimw(maxnz),  
     +           Hpsrc(maxnz),TEpsrc(maxnz),
     +           qstarm(maxnz),
     +           trip(maxnz),TpSq(maxnz),Tpdqp(maxnz)
        real ::  xkn(maxnx),ykn(maxny),
     +           xk(maxnx),yk(maxny),
     +           xks(maxnx2,maxny),wavexy(maxnx2,maxny)
        real ::  ug(maxnz),vg(maxnz),
     +           wls(maxnz),uls(maxnx)
        real ::  udrift,vdrift,
     +           stokesw,stokesa,
     +           stokess,stokes(maxnz1)
        real ::  dtg, dslg, dzg
        real ::  lwc,phiw,phiv 

        end type con_s
        type(con_s), target :: c_s
        real, pointer ::
     +           wwsb(:), engz(:), engsbz(:),
     +           englez(:), uxym(:), vxym(:), wxym(:),
     +           txym(:,:), divz(:), utle(:,:), utsb(:,:),
     +           vtle(:,:), vtsb(:,:), wtle(:,:), wtsb(:,:),
     +           wt_tot(:,:), z(:), zz(:), shrz(:), buyz(:),
     +           triz(:), uwsb(:), vwsb(:), uwle(:), vwle(:),
     +           uw_tot(:), vw_tot(:), wcube(:), wfour(:),
     +           tcube(:,:), ups(:), vps(:),
     +           wps(:), tps(:,:), t_rprod(:), t_wq(:),
     +           t_wp(:), t_tau(:), t_tran(:), t_buoy(:),
     +           t_diss(:), t_sprod(:),
     +           zconc(:),
     +           vp1mean(:),vp2mean(:),vp3mean(:),
     +           vp1msqr(:),vp2msqr(:),vp3msqr(:),
     +           m1src(:),m2src(:),m3src(:),
     +           upwpm(:),
     +           Tpmean(:),Tpmsqr(:),
     +           Tfmean(:),qfmean(:),
     +           wpTpm(:),Tpsrc(:),
     +           Tv_tran(:,:),Tv_prod(:,:),Tv_diss(:,:),
     +           Tv_part1(:,:),Tv_part2(:),
     +           radmean(:),rad2mean(:),
     +           mpmean(:),mwmean(:),
     +           phimp(:),phimw(:), 
     +           Hpsrc(:), TEpsrc(:),
     +           qstarm(:),
     +           trip(:),TpSq(:),Tpdqp(:)
        real, pointer ::  
     +           xkn(:), ykn(:), xk(:), yk(:), xks(:,:), wavexy(:,:)
        real, pointer ::  
     +           ug(:), vg(:), wls(:), uls(:)
        real, pointer ::  
     +           udrift, vdrift, stokesw, stokesa,
     +           stokess, stokes(:)
        real, pointer ::  
     +           dtg, dslg, dzg
        real, pointer ::  
     +           lwc,phiw,phiv 
      contains
         subroutine fill_cs
c
c -------------- pointer association for stat arrays and get size
c                of stat arrays isize for history files
c
             isize = 0 
             wwsb    => c_s%wwsb     ; isize = isize + size(wwsb)
             engz    => c_s%engz     ; isize = isize + size(engz)
             engsbz  => c_s%engsbz   ; isize = isize + size(engsbz)
             englez  => c_s%englez   ; isize = isize + size(englez)
             uxym    => c_s%uxym     ; isize = isize + size(uxym)
             vxym    => c_s%vxym     ; isize = isize + size(vxym)
             wxym    => c_s%wxym     ; isize = isize + size(wxym)
             txym    => c_s%txym     ; isize = isize + size(txym)
             divz    => c_s%divz     ; isize = isize + size(divz)
             utle    => c_s%utle     ; isize = isize + size(utle)
             utsb    => c_s%utsb     ; isize = isize + size(utsb)
             vtle    => c_s%vtle     ; isize = isize + size(vtle)
             vtsb    => c_s%vtsb     ; isize = isize + size(vtsb)
             wtle    => c_s%wtle     ; isize = isize + size(wtle)
             wtsb    => c_s%wtsb     ; isize = isize + size(wtsb)
             wt_tot  => c_s%wt_tot   ; isize = isize + size(wt_tot)
             z       => c_s%z        ; isize = isize + size(z)
             zz      => c_s%zz       ; isize = isize + size(zz)
             shrz    => c_s%shrz     ; isize = isize + size(shrz)
             buyz    => c_s%buyz     ; isize = isize + size(buyz)
             triz    => c_s%triz     ; isize = isize + size(triz)
             uwsb    => c_s%uwsb     ; isize = isize + size(uwsb)
             vwsb    => c_s%vwsb     ; isize = isize + size(vwsb)
             uwle    => c_s%uwle     ; isize = isize + size(uwle)
             vwle    => c_s%vwle     ; isize = isize + size(vwle)
             uw_tot  => c_s%uw_tot   ; isize = isize + size(uw_tot)
             vw_tot  => c_s%vw_tot   ; isize = isize + size(vw_tot)
             wcube   => c_s%wcube    ; isize = isize + size(wcube)
             wfour   => c_s%wfour    ; isize = isize + size(wfour)
             tcube   => c_s%tcube    ; isize = isize + size(tcube)
             ups     => c_s%ups      ; isize = isize + size(ups)
             vps     => c_s%vps      ; isize = isize + size(vps)
             wps     => c_s%wps      ; isize = isize + size(wps)
             tps     => c_s%tps      ; isize = isize + size(tps)
             t_rprod => c_s%t_rprod  ; isize = isize + size(t_rprod)
             t_wq    => c_s%t_wq     ; isize = isize + size(t_wq)
             t_wp    => c_s%t_wp     ; isize = isize + size(t_wp)
             t_tau   => c_s%t_tau    ; isize = isize + size(t_tau)
             t_tran  => c_s%t_tran   ; isize = isize + size(t_tran)
             t_buoy  => c_s%t_buoy   ; isize = isize + size(t_buoy)
             t_diss  => c_s%t_diss   ; isize = isize + size(t_diss)
             t_sprod => c_s%t_sprod  ; isize = isize + size(t_sprod)
             zconc   => c_s%zconc    ; isize = isize + size(zconc)
             vp1mean  => c_s%vp1mean ; isize = isize + size(vp1mean)
             vp2mean  => c_s%vp2mean ; isize = isize + size(vp2mean)
             vp3mean  => c_s%vp3mean ; isize = isize + size(vp3mean)
             vp1msqr  => c_s%vp1msqr ; isize = isize + size(vp1msqr)
             vp2msqr  => c_s%vp2msqr ; isize = isize + size(vp2msqr)
             vp3msqr  => c_s%vp3msqr ; isize = isize + size(vp3msqr)
             m1src   => c_s%m1src    ; isize = isize + size(m1src)
             m2src   => c_s%m2src    ; isize = isize + size(m2src)
             m3src   => c_s%m3src    ; isize = isize + size(m3src)
             upwpm   => c_s%upwpm    ; isize = isize + size(upwpm)
             Tpmean  => c_s%Tpmean   ; isize = isize + size(Tpmean)
             Tpmsqr  => c_s%Tpmsqr   ; isize = isize + size(Tpmsqr)
             Tfmean  => c_s%Tfmean   ; isize = isize + size(Tfmean)
             qfmean  => c_s%qfmean   ; isize = isize + size(qfmean)
             wpTpm   => c_s%wpTpm    ; isize = isize + size(wpTpm)
             Tpsrc   => c_s%Tpsrc    ; isize = isize + size(Tpsrc)
             Tv_tran => c_s%Tv_tran  ; isize = isize + size(Tv_tran)
             Tv_prod => c_s%Tv_prod  ; isize = isize + size(Tv_prod)
             Tv_diss => c_s%Tv_diss  ; isize = isize + size(Tv_diss)
             Tv_part1 => c_s%Tv_part1  ; isize = isize + size(Tv_part1)
             Tv_part2 => c_s%Tv_part2  ; isize = isize + size(Tv_part2)
             radmean => c_s%radmean  ; isize = isize + size(radmean)
             rad2mean => c_s%rad2mean  ; isize = isize + size(rad2mean)
             mpmean => c_s%mpmean    ; isize = isize + size(mpmean)
             mwmean => c_s%mwmean    ; isize = isize + size(mwmean) 
             phimp => c_s%phimp    ; isize = isize + size(phimp) 
             phimw => c_s%phimw    ; isize = isize + size(phimw)
             Hpsrc   => c_s%Hpsrc    ; isize = isize + size(Hpsrc)
             TEpsrc  => c_s%TEpsrc  ; isize = isize + size(TEpsrc)
             qstarm => c_s%qstarm  ; isize = isize + size(qstarm)
             trip => c_s%trip; isize = isize + size(trip)
             TpSq => c_s%TpSq; isize = isize + size(TpSq)
             Tpdqp => c_s%Tpdqp; isize = isize + size(Tpdqp)


             xkn     => c_s%xkn
             ykn     => c_s%ykn
             xk      => c_s%xk
             yk      => c_s%yk
             xks     => c_s%xks 
             wavexy  => c_s%wavexy
             ug      => c_s%ug
             vg      => c_s%vg
             wls     => c_s%wls
             uls     => c_s%uls
             udrift  => c_s%udrift
             vdrift  => c_s%vdrift
             stokesw => c_s%stokesw
             stokesa => c_s%stokesa
             stokess => c_s%stokess
             stokes  => c_s%stokes
             dtg     => c_s%dtg
             dslg    => c_s%dslg 
             dzg     => c_s%dzg
             lwc    => c_s%lwc
             phiw    => c_s%phiw 
             phiv    => c_s%phiv
         return
         end subroutine fill_cs
      end module con_stats
c ----------------------------------------------------------------------
      program les_mpi_2d
c
      use pars
      use fields
      use particles
      use con_data
      use con_stats
      include 'mpif.h'
c
c ------------- definition of internal flags
c
c
c       iDNS    = 0; call the subgrid computation of vis_m and vis_s
c               = 1; call the molecular viscosity and diffusivity
c
c       igrdr   =  3; data comes from restart file
c               =  2; data comes from initialization (random)
c               =  1; data comes from coarser grid (or otherwise)
c
c       ibcu    = 1 ; upper boundary condition set by radiation bc
c               = 0 ; fixed value = 0.
c               = -1; value defined by coarser mesh for all variables
c
c       ibcl    = 0 ; lower boundary condition set by similarity theory (sr. setup)
c               = -1; value defined by coarser mesh for all variables
c
c       ifix_dt = 0 ; variable time step with fixed cfl number in setcon
c               = 1 ; fixed time step set in sr. get_dt
c
c       ifree   = 0 ; use spatially averaged surface conditions for MO (call lower)
c               = 1 ; use point-by-point conditions for MO free convection (call lower_free)
c
c       ihst    = nn; frequency at which global variables are output in history file
c               < 0 ; no history files
c
c       it_his  = time step where history files start, incremented by itape
c
c       it_viz  = time step where viz files start, incremented by itape
c
c       ismlt   = 1 ; use businger formulas in MO 
c                 0 ; use large and everyone elses formulas in MO 
c
c       iupwnd  = 0;  use skew symmetric formulas for all derivatives
c                     in scalar equations
c               = 1;  use hybrid upwind scheme for all derivatives
c                     in scalar equations
c
c       ivis0   = 0; old eddy viscosity model 
c               = 1; new eddy viscosity model 
c
c       new_vis = step; the iteration step for which the new model
c                       is turned on when ivis0=1
c               < 0; new model is on at all steps for ivis0=1
c
c       nscl  .ge. 1   number of scalars to be followed set in parameter statements
c                      change entries in sr. init, and sr. suft for surface bc's
c
c -------------------------------------------------------------------------------
c    ---------- Define PI()--------------

      pi = 4.0*atan(1.0)

c
c ---------- initialize MPI, get myid, numprocs, 
c            test if on root process
c
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,myid,ierr)
      call mpi_comm_size(mpi_comm_world,numprocs,ierr)
c
      i_root = 0
      l_root = .false.
      if(myid .eq. i_root) l_root = .true.
c
      l_debug = .false.
      if(idebug .eq. 1) l_debug = .true.
c
      ts_mpi = mpi_wtime()

!----- Read the input file for all necessary parameters
      call read_input_file

c
c -------- set number of x-y slab cpus
c
!      ncpu_s = 8
c
      case_inp = 'cou'
c
      call get_units
      call gridd
      call setcon
c
c -------------- scratch run
c
      if (iti.eq.0)  then
         igrdr = 2
         case = case_inp
         call init
         call setup(it)

         call particle_setup
         call particle_init 
         !call read_part_res !comment out usually 

c
c ---------- choose routine for getting initial guess
c
         if(iocean .eq. 1) then
            call randoc
         else
            if (ifields .eq. 1) then
            call get_fields
            else
            call random   
            end if
         endif
         !Call dns_vis even when doing LES since vis_m needs to be initialized with something
         call dns_vis
         call get_max
      else

         igrdr = 3
         call restart
         call setup(it)

         !Call dns_vis even when doing LES since vis_m needs to be initialized with something
         call dns_vis
         call get_max

 
         call particle_setup
         call read_part_res
      endif


c
c --------------- time loop ------------
c
      tzero = time
      call get_dt
 9000 continue
      call set_sav(it,iti)
      if (myid==0 .and. it==0) then
      write(*,*) 'Starting time loop'
      end if
      if (myid==0) then
      write(*,*) 'it,time = ',it,time
      end if
      if (it == 1) numpart = 0
      part => first_particle
      do while (associated(part))
      !particle sample 1
      if (part%pidx == 1 .and. part%procidx == 0) then
      write(*,'(a6,4e15.6)') 'xp1:  ',time,part%xp(1:3)
      write(*,'(a6,4e15.6)') 'vp1:  ',time,part%vp(1:3)
      write(*,'(a6,4e15.6)') 'uf1:  ',time,part%uf(1:3)
      write(*,'(a6,2e15.6)') 'Tp1:  ',time,part%Tp
      write(*,'(a6,2e15.6)') 'Tf1:  ',time,part%Tf
      write(*,'(a6,2e15.6)') 'rad1: ',time,part%radius 
      write(*,'(a6,2e15.6)') 'qinf1:',time,part%qinf 
      write(*,'(a6,2e15.6)') 'qstr1:',time,part%qstar 
      write(*,'(a6,e15.6,i)') 'mult1:',time,part%mult
      end if

      !particle sample 2
      if (part%pidx == 2 .and. part%procidx == 0) then
      write(*,'(a6,4e15.6)') 'xp2:  ',time,part%xp(1:3)
      write(*,'(a6,4e15.6)') 'vp2:  ',time,part%vp(1:3)
      write(*,'(a6,4e15.6)') 'uf2:  ',time,part%uf(1:3)
      write(*,'(a6,2e15.6)') 'Tp2:  ',time,part%Tp
      write(*,'(a6,2e15.6)') 'Tf2:  ',time,part%Tf
      write(*,'(a6,2e15.6)') 'rad2: ',time,part%radius 
      write(*,'(a6,2e15.6)') 'qinf2:',time,part%qinf 
      write(*,'(a6,2e15.6)') 'qstr2:',time,part%qstar 
      write(*,'(a6,e15.6,i)') 'mult2:',time,part%mult
      end if

      !particle sample 3
      if (part%pidx == 1 .and. part%procidx == 1) then
      write(*,'(a6,4e15.6)') 'xp3:  ',time,part%xp(1:3)
      write(*,'(a6,4e15.6)') 'vp3:  ',time,part%vp(1:3)
      write(*,'(a6,4e15.6)') 'uf3:  ',time,part%uf(1:3)
      write(*,'(a6,2e15.6)') 'Tp3:  ',time,part%Tp
      write(*,'(a6,2e15.6)') 'Tf3:  ',time,part%Tf
      write(*,'(a6,2e15.6)') 'rad3: ',time,part%radius 
      write(*,'(a6,2e15.6)') 'qinf3:',time,part%qinf 
      write(*,'(a6,2e15.6)') 'qstr3:',time,part%qstar 
      write(*,'(a6,e15.6,i)') 'mult3:',time,part%mult
      end if
      part => part%next

      if (it == 1) numpart = numpart + 1
      end do
      if(myid==0) then

      write(*,'(a15,3e15.6)')'avg U(zl/2),U(zl):',
     +time,uxym(maxnz/2),uxym(maxnz)
      write(*,'(a15,3e15.6)')'avg T(zl/2),T(zl):',
     +time,txym(maxnz/2,1),txym(maxnz,1)
      write(*,'(a15,e15.6,i12)')'tnumpart:      ',time,tnumpart
      write(*,'(a15,e15.6,i12)')'tnum_des:      ',time,tnum_destroy
      write(*,'(a15,3e15.6)')   'radmin,radmax: ',time,radmin,radmax
      write(*,'(a15,3e15.6)')   'tempmin,tempma:',time,tempmin,tempmax
      write(*,'(a15,3e15.6)')   'qmin,qmax:     ',time,qmin,qmax
      write(*,'(a15,2e15.6)')   'radavg:        ',time,radavg
      write(*,'(a15,e15.6,2i12)')'t,100,1000:    ',time,tnum100,tnum1000
      write(*,'(a15,e15.6,i12)')  'impos:         ',time,tnumimpos
      write(*,'(a15,e15.6,i12)')  'actnum:        ',time,tactnum
      write(*,'(a15,e15.6,i12)')  'denum:         ',time,tdenum
      end if

      if(it .ge. new_vis .and. ivis0 .eq. 1) then
          ivis = 1
      else
          ivis = 0
      endif


c
c ---------------- 3 stage runge-kutta time stepping
c
      t_stage_s = mpi_wtime()
      do  8999 istage=1,3
c
      dtzeta = dt*zetas(istage)
      dtgama = dt*gama(istage)
c
c ---------- compute derivatives of (u,v,w)
c

      call exchange
c

      call get_derv

c
c --------- new eddy viscosity, and bcs
c
      if(iss .eq. 0 .and. ifree .eq. 0) then
         if (iDNS .eq. 1) then
            call lower_dns(it)
         else
            call lower(it)
         end if
      elseif(ifree .eq. 1) then
         if (iDNS .eq. 1) then
            call lower_dns(it)
         else
            call lower_free(it)
         end if
      endif
      if(ise .eq. numprocs-1) then
         if (iDNS .eq. 1) then
            call upper_dns
         else
            call upper
         end if
      endif
      call bcast_pbc
      call get_means(istage)
      if(ivis .eq. 1) then
         call iso(it)
         call surfvis(it)
      endif
      if(istage .eq. 1)then
        call xy_stats
        call tke_budget
        do iscl=1,nscl
        call Tvar_budget(iscl)
        enddo
        call extra_flux_terms
        call pbltop(itop)
      endif


c
c ------------ save velocity field
c
      if(msave .and. istage .eq. 1) then
         call save_v(it)
         call save_particles
      endif
      if(msave_v .and. istage .eq. 1) then
         call save_viz(it)
      endif
c
c ------------ save pressure field
c
      if(msave .and. istage .eq. 1) then
         call save_p
      endif

      if(mhis  .and. istage .eq. 1) then
          call write_histograms
      end if

c
c --------- get rhs for all equations
c
      !t_s = mpi_wtime()
      call comp1(istage,it)
      !t_f = mpi_wtime()
      !write(*,*) 'proc, comp1: ',myid,t_f-t_s
      !call mpi_barrier(mpi_comm_world,ierr)
      !if (myid==5) write(*,*) 'time comp1: ',t_f-t_s
      if(istage .eq. 1) then
         if(msave .and. l_root) call save_c(it)
      endif
c
c --------- solve for pressure
c
      !t_s = mpi_wtime()
      call comp_p
      !t_f = mpi_wtime()
      !write(*,*) 'proc, comp_p: ',myid,t_f-t_s
      !call mpi_barrier(mpi_comm_world,ierr)
      !if (myid==5) write(*,*) 'time comp_p: ',t_f-t_s
c
c --------- add pressure gradient and dealias
c
      !t_s = mpi_wtime()
      call comp2
      !t_f = mpi_wtime()
      !!write(*,*) 'proc, comp2: ',myid,t_f-t_s
      !call mpi_barrier(mpi_comm_world,ierr)
      !if (myid==5) write(*,*) 'time comp2: ',t_f-t_s
      if(micut) then
         call dealias
      endif



c      
c -------- update particles
c
         if (ipart_method .eq. 1) then
         t_s = mpi_wtime()
         !call particle_update_cn(it,istage) 
         call particle_update_rk3(it,istage) 

         t_f = mpi_wtime()
c        call mpi_barrier(mpi_comm_world,ierr)
c        if (myid==5) write(*,*) 'time part: ',t_f-t_s
         end if


      if(mnout .and. istage .eq. 1)  then
          if(l_debug) then
             call print(nprt,it,izs,ize)
          endif
          if(l_root) call print(6,it,1,nnz)
      endif
      if(l_root) then
         if(mhis  .and. istage .eq. 1)  call write_his(itop)
         if(mhis  .and. istage .eq. 1 .and. mtape) call close_his
      endif
c
 8999 continue
      t_stage_f = mpi_wtime()
      !if (myid==5) write(*,*) 'time stage: ',t_stage_f - t_stage_s

      if (ipart_method .eq. 2) then
         t_s = mpi_wtime()

         call particle_update_BE(it)

         t_f = mpi_wtime()
         call mpi_barrier(mpi_comm_world,ierr)
         if (myid==5) write(*,*) 'time part: ',t_f-t_s
      end if

      !Call coalescence outside of RK loop since it's not appropriate as
      !part of RK scheme
      if (icoalesce) then
      t_s = mpi_wtime()

      call particle_coalesce

      call mpi_barrier(mpi_comm_world,ierr)
      t_f = mpi_wtime()
      if (myid==5) write(*,*) 'time coalescence:', t_f - t_s
      end if

      if (ipartdiff) then
      t_s = mpi_wtime()

        if (isfs == 2) then
           !Call Weil et al. (2004) Lagrangian stochastic model
           call SFS_velocity
        elseif(isfs == 1) then
           !Call stochastic model for the position based on vis_s
           call SFS_position
        end if

      call mpi_barrier(mpi_comm_world,ierr)
      t_f = mpi_wtime()
      if (myid==5) write(*,*) 'time partdiff:', t_f - t_s
      end if

      !call radius_histogram
     
      call get_max
      call get_dt
      if (it.ge.itmax) go to 99000
      go to 9000
c
99000 continue
      te_mpi = mpi_wtime()
      write(6,9997) (te_mpi - ts_mpi)
 9997 format(' Job Execution Time = ',e15.6)
c
 9998 continue
      call mpi_finalize(ierr)
c
      stop
      end

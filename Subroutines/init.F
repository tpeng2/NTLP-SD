      subroutine init
c
      use pars
      use fields
      use con_data
      use con_stats
      use particles
c
      real :: RHT,RHB,mod_Magnus

      pi2  = 2.0*pi
      bfac = 1.0
      if(ibuoy.eq.0) bfac = 0.
c
c -------------------- case specific data
c
      if(iocean .eq. 1) then
         t00     = 283.
         t00b    = 5000.0
         cp      = 4.20e03
         gcp     = grav/cp
         batag   = bfac*grav/t00b
c        fcor    = 0.0
         fcor    = 1.39e-04
         fcor_h  = 0.0
c        wtsfc(1)=0.00
c        wtsfc(1)=4.96e-07
         wtsfc(1)=1.190476e-06
         qstar(1)=wtsfc(1)
c        dtdzf(1)=0.000
         dtdzf(1)=0.2548
         dtjump  = 0.
         divgls  = 0.
         zo      = 0.0001
         zi      = -5.
c        izi     = (55*nnz)/100
c        izi     = nnz
         izi     = 55
         xl      = 50.
         yl      = 50.
         zl      = -20.
c
c ---------- if stretched grid specify location of first point
c
         zw1 = -0.5
      else
        !--- if (icoean.ne.1)
         gcp     = grav/Cpa
         batag   = bfac*grav/t00
         fcor_h  = 0.0

         !!!If doing DNS, BCs are set somewhat differently than when doing LES
         if (iDNS .eq. 1) then
           !Get BCs of temperature in Celsius for use in specific expressions
           topC=(Ttop(1)-273.15)
           botC=(Tbot(1)-273.15)   
           !Normalize the RH given in input file
           RHT = Ttop(2)/100.0 
           RHB = Tbot(2)/100.0
           !Convert RH given in input file into specific humidity
           Ttop(2) =  RHT*Mw/Ru/Ttop(1)*610.94*
     +                 EXP(17.6257*topC/(topC+243.04))/rhoa
           Tbot(2) =  RHB*Mw/Ru/Tbot(1)*610.94*
     +                 EXP(17.6257*botC/(botC+243.04))/rhoa
           wtsfc(1) = qstar(1)
           wtsfc(2) = qstar(2)

         else !(doing LES)
           !qstar specified in input file for isfc = 0
           wtsfc = qstar

           !tsfcc specified in input file for isfc = 1

           !Assuming tsfcc(2) is giving RH:
           surf_RH=tsfcc(2)
           RHB = surf_RH/100.0
           !Convert to specific huminity
           tsfcc(2) = RHB*Mw/Ru/tsfcc(1)*mod_Magnus(tsfcc(1))/rhoa

           if (myid==0) then
              write(*,*) 'SURFACE HUMIDITY: %RH, q = ',surf_RH,tsfcc(2)
           end if
      end if
c
c
         !Gradients of temperature and scalars above inversion (upper BC)
         dtdzf(1)=0.003
         dtdzf(2)=0.000
c
         dtjump  = 0.0
         divgls  = 0.0

  
c
c ---------- if stretched grid specify location of first point
c
         !Particle channel:
         izi     = (60*nnz)/100
      endif
c
      time  = 0.0
c 
c ---------- outermost coarse grid  indicies are bounds of grid
c
      izlow = 1
      izup  = nnz
      dz    = zl/nnz
      dzg   = abs(dz)
      if(l_root) write(6,4040) zl,nnz,dzg
c
c --------------- generate z grids for particular mesh from
c                 iz = 0,1,...,nnz+1; this allows indexing
c                 to array elements z(0), etc.
c
      zwstrt = 0.0
c
c ------------ if uniform vertical spacing then
c
      if(iz_space .eq. 0) then
c
c ------------ build z grid for w points
c
         do iz=0,nnz+1
            z(iz) = dz*float(iz) + zwstrt
         enddo
      else
       if(ichannel.eq.1) then
       call vgrid_channel(zw1,zi,zl,nnz,z(0),l_root,l_debug)
       else if (ichannel .eq.0) then
       call vgrid(zw1,zi,zl,nnz,z(0),l_root,l_debug)
       else if (ichannel .eq.2) then
       call vgrid_channel_fstrm(zw1,zi,zl,nnz,z(0),l_root,l_debug) !Tianze 9/1/15
       end if
      endif
c
      call get_dz
c
      if(l_root) then
         write(6,8002) zwstrt
         write(6,8003) (iz,z(iz),zz(iz),iz=0,nnz+1)
      endif
c
      nnzm1 = nnz-1
      dx    = xl/nnx
      dy    = yl/nny
      fnxy  = 1./float(nxy)
      dzdz  = dzw(1)*dzw(1)
      z1    = zz(1)
c
      c23  = 2.0/3.0
      dsl  = (dx*1.5*dy*1.5*abs(dzw(1)))**(1./3.)
      dslg = dsl
      cs   = 0.2
c
      vk     = 0.4
      batagk = batag*vk
      vkin   = 1./vk
      ttmean = 0.
      zody   = alog(abs(z1/zo))
      write(nprt, 9901) z1,zo,zody
 9901 format(' 9901 z1 = ',e15.6,' zo = ',e15.6,/,
     +       ' zody = ',e15.6)
      zodyin = 1./zody
      wstar  = abs(batag*zi*wtsfc(1))**(1./3.)
      if(ismlt .eq. 1) then
c
c ---- set constants for businger similarity functions
c
         vk74   = vk*0.74
         vk74in = 0.74/vk
         zody74 = zody*0.74
      else 
c
c ---- set constants for large similarity functions
c
        vk74    = vk
        vk74in  = 1.0/vk
        zody74  = zody
      endif
      ugal   = 0.0
c      ugal   = ugcont*0.5
c     ugcont = ugcont - ugal
      cdbtm  = vk*vk/zody/zody
      if(iocean .eq. 1) then
c ----------- set surface friction velocity here and in sr. sufto
c        utau = 4.29e-03
         utau = 7.00e-03
      else
         ufree = 0.07*(abs(batag*qstar(1)*dzw(1)))**(1./3.)
c
c ---- note : new estimate for utau !!!
c
         utau  = vk*(ufree+ugcont)/zody
c        utau  = vk*(ufree)/zody
      endif
      utau2    = utau*utau
      if(ibuoy .eq. 0 .or. qstar(1) .eq. 0.) then
        amonin = 1000.0
      else
        amonin = -utau2*utau/(batagk*qstar(1))
      endif
      hol   = abs(zi)/amonin
      zol   = abs(z1)/amonin
      uwsfc = -utau*utau
      vwsfc = -utau*utau
c
      if(l_root) then
         write(6,80)
         write(6,2)wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +         ,cdbtm,ugcont
      endif
c
      if(l_debug) then
         write(nprt,80)
         write(nprt,2)wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +         ,cdbtm,ugcont
      endif
c
      return
c ------------------------
   2  format(10x,' WT =',e12.4,',  U* =',e12.4,',  L =',e12.4,/,
     +       10x,' DTDZ FREE =',e12.4,',  ZODY=',e12.4,/,10x,
     +       ' ZO(BTM) =',e12.4,',  CDBTM=',e12.4,
     +       ',  UG = ',e12.4)
  80  format(///,' ***** SCRATCH RUN ***** ',//)
 4040 format(' zl = ',e15.6,' nnz = ',i5,' dzg = ',e15.6)
 4043 format(' znest = ',e15.6,' nnz = ',i5,' dzg = ',e15.6)
 8002 format(' zwstrt = ',e12.4)
 8003 format(' iz ',5x,' zw',5x,' zu ',5x,/,(i3,2e12.4))
      end
      subroutine randoc
c
c -------- random initial conditions for an
c          ocean simulation
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real psi(nnx,iys:iye), psix(nnx,iys:iye),
     +     psiy(nnx,iys:iye,izs:izs), uxx(nnx,iys:iye),
     +     vyy(nnx,iys:iye,izs:izs)
c
      izi=(5*nnz)/20
      zi=z(izi)
      tmixed = 283.0
      do iz=izs,ize
         if (iz.le.izi) then
            do iy=iys,iye
            do ix=1,nnx
               u(ix,iy,iz)   = ugcont-ugal
               v(ix,iy,iz)   = vgcont
               w(ix,iy,iz)   = 0.0
               t(ix,iy,1,iz) = tmixed
               e(ix,iy,iz)   = 0.0
            enddo
            enddo
         endif
         if (iz.gt.izi) then
            do iy=iys,iye
            do ix=1,nnx
               u(ix,iy,iz)   = ugcont-ugal
               v(ix,iy,iz)   = vgcont
               w(ix,iy,iz)   = 0.0
               t(ix,iy,1,iz) = tmixed + dtdzf(1)*(zz(iz)-zi)
               e(ix,iy,iz)   = 0.0
            enddo
            enddo
         endif
         do iy=iys,iye
         do ix=1,nnx
            w(ix,iy,iz)    = 0.0
            r1(ix,iy,iz)   = 0.0
            r2(ix,iy,iz)   = 0.0
            r3(ix,iy,iz)   = 0.0
            r4(ix,iy,1,iz) = 0.0
            r5(ix,iy,iz)   = 0.0
         enddo
         enddo
      enddo
c
c ------------- set initial random field to be
c               divergence free
c
      idum = -1
      do iz=izs,ize
      if (iz.le.4) then
c
c ----------- ampv and ampt are max amplitudes of random 
c             velocity and temperature fields
c
         ampv = 0.01
c        ampt = 0.00
         ampt = 0.0001
c  
c ------- simple random field scaled between 0 and 1
c
         do iy=iys,iye
         do ix=1,nnx
            psi(ix,iy) = ran1(idum)
         enddo
         enddo
c
         do iy=iys,iye
         do ix=1,nnx
            psix(ix,iy) = psi(ix,iy)
            psiy(ix,iy,izs) = psi(ix,iy)
         enddo
         enddo
         call xderivp(psix(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(psiy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
c
         vmaxx = 0.0
         do iy=iys,iye
         do ix=1,nnx
            vmag = sqrt(psix(ix,iy)**2 + psiy(ix,iy,izs)**2)
            if(vmag .gt. vmaxx) vmaxx = vmag
         enddo
         enddo
         facv = ampv/vmaxx
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz) = u(ix,iy,iz) - psiy(ix,iy,izs)*facv
            v(ix,iy,iz) = v(ix,iy,iz) + psix(ix,iy)*facv
            t(ix,iy,1,iz) = t(ix,iy,1,iz) + psi(ix,iy)*ampt
            e(ix,iy,iz) = 0.0001
         enddo
         enddo
      endif
c
c ---------- check divergence of initial field
c
         do iy=iys,iye
         do ix=1,nnx
            uxx(ix,iy) = u(ix,iy,iz)
            vyy(ix,iy,izs) = v(ix,iy,iz)
         enddo
         enddo
         call xderivp(uxx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(vyy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
         do iy=iys,iye
         do ix=1,nnx
            divz(iz) = divz(iz) + (uxx(ix,iy) + vyy(ix,iy,izs))**2
         enddo
         enddo
         divz(iz) = divz(iz)*fnxy
c
c -------- end z loop
c
      enddo
c
      call mpi_sum_z(divz(1),i_root,myid,nnz,1)
c
      write(nprt,6000)
 6000 format(' check of divergence for initial state',/,
     +         ' iz ',5x,' divergence')
      write(nprt,6100) (iz,divz(iz),iz=izs,ize)
 6100 format(i5,e15.6)
c
      do iz=izs,ize
         ug(iz)=ugcont
         vg(iz)=vgcont
      enddo
c
      return
      end
      subroutine random
c
c ----------- geostrophic winds designed for comparison case
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use particles
      real :: hsurf_init
      real psi(nnx,iys:iye), psix(nnx,iys:iye),
     +     psiy(nnx,iys:iye,izs:izs), uxx(nnx,iys:iye),
     +     vyy(nnx,iys:iye,izs:izs),mod_Magnus,Ttmp

      !Initialize partcount to 0:
      partcount = 0.0
      partcount_t = 0.0
      vpsum = 0.0
      vpsum_t = 0.0
      vpsqrsum = 0.0
      vpsqrsum_t = 0.0
      upwp_t = 0.0
      upwp = 0.0
      Tpsum = 0.0
      Tpsum_t = 0.0
      Tpsqrsum = 0.0
      Tpsqrsum_t = 0.0
      Tfsum = 0.0
      Tfsum_t = 0.0
      qfsum = 0.0
      qfsum_t = 0.0
      wpTpsum = 0.0
      wpTpsum_t = 0.0
      partsrc = 0.0
      partsrc_t = 0.0
      partTsrc = 0.0
      partTsrc_t = 0.0
      partHsrc = 0.0
      partHsrc_t = 0.0
      partTEsrc = 0.0
      partTEsrc_t = 0.0
      radsum = 0.0 
      radsum_t = 0.0 
      rad2sum = 0.0 
      rad2sum_t = 0.0 
      mpsum = 0.0
      mpsum_t = 0.0
      mwsum = 0.0
      mwsum_t = 0.0
      qstarsum = 0.0
      qstarsum_t = 0.0 

c ------------ note set nmatch in sr. iso so that
c              it is compatible with conditions here
c
      do iz=1,nnz
c        ug(iz)   = ugcont*(zz(iz)/zl)
         ug(iz)   = ugcont
         vg(iz)   = vgcont
         divz(iz) = 0.0
      enddo
c
      izi = (50*nnz)/100
      zi  = z(izi)
      if(myid==0) write(*,*)'The inversion layer initalized at ',zi,' m'
c
      z_lower = zi - 50.0
      t_lower = 300.0
      z_upper = zi + 50.0
      t_upper = 308.0
      slope   = (t_upper - t_lower)/(z_upper - z_lower)
c
c     !--- assume surface layer at 50m 
      hsurf_init=50
       do iy=iys,iye
         do ix=1,nnx
           do iz=izs,ize
            !--- initialize (u,v,w,e) fileds 
            u(ix,iy,iz) = ugcont-ugal
c           u(ix,iy,iz) = ug(iz) - ugal
            if (z(iz).le.hsurf_init.and.zw1/zo.gt.1) then
               u(ix,iy,iz)=u(ix,iy,iz)/log(hsurf_init/zo)*log(z(iz)/zo)
            endif
            v(ix,iy,iz) = vgcont
            w(ix,iy,iz) = 0.0
            e(ix,iy,iz) = 0.0
            !--- initialize Temp fileds with inversion layer specified 
            if(z(iz) .le. z_lower) then
              t(ix,iy,1,iz) = t_lower
              if (itempfld==1) then
                  !Add "bubble" to test cloud model
                  t(ix,iy,1,iz) = t(ix,iy,1,iz) + 0.1*
     +            exp( -( (dx*(ix-1)-500.0)/200.0 )**2 
     +                 -( (dy*(iy-1)-500.0)/200.0 )**2 )*
     +            exp (-zz(iz)/100.0)
              endif
            elseif(z(iz) .ge. z_upper) then
              t(ix,iy,1,iz) = t_upper + (zz(iz+1) - z_upper)*dtdzf(1)
            else
              t(ix,iy,1,iz) = t_lower + slope*(zz(iz+1) - z_lower)
            endif ! for(z(i))
            !--- initialize specific humidity fileds based on Temp.
            !Make a constant q profile at RH based on surface temp:
            t(ix,iy,2,iz) = surf_RH/100.0*
     +      Mw/Ru/t_lower*mod_Magnus(t_lower)/rhoa

         !Make a desired profile of RH:
!            !Convert theta to temperature:
!            Ttmp = 
!     +      t(ix,iy,1,iz)*(psurf/(psurf - zz(iz)*rhoa*grav))**(-Rd/Cpa)
!
!            !Get temperature based on desired RH
!            t(ix,iy,2,iz) = surf_RH/100.0*
!     +      Mw/Ru/Ttmp*mod_Magnus(Ttmp)/rhoa
            r1(ix,iy,iz)  = 0.
            r2(ix,iy,iz)  = 0.
            r3(ix,iy,iz)  = 0.
            r4(ix,iy,1,iz)= 0.
            r4(ix,iy,2,iz)= 0.
            r5(ix,iy,iz)  = 0.
           enddo 
         enddo 
      enddo
c
c ------------- set initial random field to be
c               divergence free
c
      idum = -1 - myid
      do iz=izs,ize
c
c ----------- ampv and ampt are max amplitudes of random 
c             velocity and temperature fields
c             make sure ampv is set if free convection so
c             that we have motions at first time step
c
         ampv = 0.0
         ampv = 0.001
         ampt = 0.10
c  
c ------- simple random field scaled between -0.5 and 0.5
c
         sum_psi = 0.0
         do iy=iys,iye
         do ix=1,nnx
            psi(ix,iy) = ran1(idum)
            sum_psi = sum_psi + psi(ix,iy)
         enddo
         enddo
         sum_psi = sum_psi*fnxy
         call mpi_sum_xy(sum_psi,myid,iss,ise,1)
         do iy=iys,iye
         do ix=1,nnx
            psi(ix,iy) = psi(ix,iy) - sum_psi
            psix(ix,iy)     = psi(ix,iy)
            psiy(ix,iy,izs) = psi(ix,iy)
         enddo
         enddo
         call xderivp(psix(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(psiy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
         vmaxx = 0.0
         do iy=iys,iye
         do ix=1,nnx
            vmag = sqrt(psix(ix,iy)**2 + psiy(ix,iy,izs)**2)
            if(vmag .gt. vmaxx) vmaxx = vmag
         enddo
         enddo
         facv = ampv/vmaxx
c
         if (z(iz) .le. 50.0) then
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz)   = u(ix,iy,iz) - psiy(ix,iy,izs)*facv
            v(ix,iy,iz)   = v(ix,iy,iz) + psix(ix,iy)*facv
            t(ix,iy,1,iz) = t(ix,iy,1,iz) + psi(ix,iy)*ampt
         enddo
         enddo
         endif
c
         if(z(iz) .le. 250.0) then
         do iy=iys,iye
         do ix=1,nnx
            e(ix,iy,iz) = 0.4*(1.0 - z(iz)/250.0)**3
         enddo
         enddo
         endif
c
c ---------- check divergence of initial field
c
         do iy=iys,iye
         do ix=1,nnx
            uxx(ix,iy) = u(ix,iy,iz)
            vyy(ix,iy,izs) = v(ix,iy,iz)
         enddo
         enddo
         call xderivp(uxx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(vyy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
         do iy=iys,iye
         do ix=1,nnx
            divz(iz) = divz(iz) + (uxx(ix,iy) + vyy(ix,iy,izs))**2
         enddo
         enddo
         divz(iz) = divz(iz)*fnxy
c
c -------- end z loop
c
      enddo
c
      call mpi_sum_z(divz(1),i_root,myid,nnz,1)
c
      write(nprt,6000)
 6000 format(' check of divergence for initial state',/,
     +         ' iz ',5x,' divergence')
      write(nprt,6100) (iz,divz(iz),iz=izs,ize)
 6100 format(i5,e15.6)
c
c ------------ fix for baroclinic and subsidence effects !!
c
c     do iz=izs,ize
c        ug(iz)=ugcont
c        vg(iz)=vgcont
c        if (.not.(ibrcl.eq.1)) go to 19988
c        if (.not.(iz.le.izi)) go to 19987
c        ug(iz)=0.
c        vg(iz)=0.
c 19987    continue
c 19988    continue
c        zz2=zz(iz)
c        wls(iz)=-divgls*zz2
c        if (.not.(iz.eq.1)) go to 19986
c        do ix=1,nnx
c        uls(ix)=divgls*(dx*float(ix-1)-xl*.5)
c        enddo
c     enddo
c     write(nprt,9)(uls(ix),ix=1,nnx)
c  9  format(1x,8e12.3)
c 19986 continue
c
      return
      end
      subroutine random_f
c
c ---------- example of using given (sparse) initial 
c            sounding profiles (FIX for ncpu_s).
c            
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real psi(nnx,iys:iye), psix(nnx,iys:iye),
     +     psiy(nnx,iys:iye,izs:izs), uxx(nnx,iys:iye),
     +     vyy(nnx,iys:iye,izs:izs)
c
      parameter (nt=12, nz=11)
      real zg(nz), u_i(nz,nt), v_i(nz,nt), theta_i(nz,nt)
      real ui_temp(nz), vi_temp(nz), ti_temp(nz)
      real time_g(nt)
c
      data time_g /
     +  0.0000E+00,  0.3600E+04,  0.7200E+04,  0.1080E+05,  0.1440E+05,
     +  0.1800E+05,  0.2160E+05,  0.2520E+05,  0.2880E+05,  0.3240E+05,
     +  0.3600E+05,  0.3960E+05
     +/
      data zg /
     +  0.1000E+02,  0.3000E+02,  0.5500E+02,  0.9000E+02,  0.1400E+03,
     +  0.2150E+03,  0.3300E+03,  0.5000E+03,  0.7500E+03,  0.1100E+04,
     +  0.1600E+04
     +/
      data u_i /
     + -0.1510E+01, -0.1560E+01, -0.1580E+01, -0.1580E+01, -0.1560E+01,
     + -0.1530E+01, -0.1510E+01, -0.9000E+00, -0.1390E+01, -0.1220E+01,
     + -0.5100E+00,
     + -0.1090E+01, -0.1110E+01, -0.1120E+01, -0.1120E+01, -0.1030E+01,
     + -0.9900E+00, -0.9500E+00, -0.6200E+00, -0.1230E+01, -0.9400E+00,
     +  0.2800E+00,
     + -0.9100E+00, -0.9200E+00, -0.9100E+00, -0.9000E+00, -0.8800E+00,
     + -0.8400E+00, -0.8000E+00, -0.6500E+00, -0.1510E+01, -0.1070E+01,
     +  0.2400E+00,
     + -0.8900E+00, -0.8900E+00, -0.8900E+00, -0.8800E+00, -0.8700E+00,
     + -0.8500E+00, -0.8100E+00, -0.7000E+00, -0.1830E+01, -0.8400E+00,
     +  0.3500E+00,
     + -0.1250E+01, -0.1260E+01, -0.1260E+01, -0.1250E+01, -0.1240E+01,
     + -0.1220E+01, -0.1160E+01, -0.8800E+00, -0.1980E+01, -0.1900E+00,
     +  0.7500E+00,
     + -0.1800E+01, -0.1810E+01, -0.1820E+01, -0.1820E+01, -0.1800E+01,
     + -0.1780E+01, -0.1710E+01, -0.1150E+01, -0.1960E+01,  0.3900E+00,
     +  0.9200E+00,
     + -0.2110E+01, -0.2130E+01, -0.2140E+01, -0.2140E+01, -0.2130E+01,
     + -0.2110E+01, -0.2050E+01, -0.9300E+00, -0.1400E+01,  0.8800E+00,
     +  0.9600E+00,
     + -0.2250E+01, -0.2280E+01, -0.2290E+01, -0.2300E+01, -0.2290E+01,
     + -0.2260E+01, -0.2070E+01, -0.4000E-01, -0.1600E+00,  0.1440E+01,
     +  0.1190E+01,
     + -0.2160E+01, -0.2200E+01, -0.2220E+01, -0.2220E+01, -0.2220E+01,
     + -0.2190E+01, -0.1610E+01,  0.1470E+01,  0.1420E+01,  0.2050E+01,
     +  0.1610E+01,
     + -0.2230E+01, -0.2270E+01, -0.2290E+01, -0.2300E+01, -0.2300E+01,
     + -0.2260E+01, -0.1350E+01,  0.2480E+01,  0.2380E+01,  0.2320E+01,
     +  0.1740E+01,
     + -0.1890E+01, -0.1930E+01, -0.1950E+01, -0.1950E+01, -0.1940E+01,
     + -0.1890E+01, -0.1120E+01,  0.3010E+01,  0.3030E+01,  0.2800E+01,
     +  0.2000E+01,
     + -0.1210E+01, -0.1230E+01, -0.1240E+01, -0.1230E+01, -0.1210E+01,
     + -0.1140E+01, -0.4600E+00,  0.3320E+01,  0.3510E+01,  0.3420E+01,
     +  0.2340E+01
     +/
      data v_i /
     +  0.4800E+00,  0.5100E+00,  0.5300E+00,  0.5700E+00,  0.6900E+00,
     +  0.7300E+00,  0.7600E+00,  0.1410E+01, -0.4200E+00, -0.3060E+01,
     + -0.3500E+01,
     +  0.7800E+00,  0.8100E+00,  0.8400E+00,  0.8900E+00,  0.1060E+01,
     +  0.1110E+01,  0.1130E+01,  0.1190E+01, -0.1040E+01, -0.2900E+01,
     + -0.3440E+01,
     +  0.3000E+00,  0.3200E+00,  0.3400E+00,  0.3800E+00,  0.4800E+00,
     +  0.5300E+00,  0.5800E+00,  0.5300E+00, -0.1330E+01, -0.2040E+01,
     + -0.2830E+01,
     + -0.2700E+00, -0.2600E+00, -0.2400E+00, -0.2200E+00, -0.1800E+00,
     + -0.1300E+00, -0.5000E-01,  0.1000E+00, -0.1170E+01, -0.1100E+01,
     + -0.2370E+01,
     + -0.5500E+00, -0.5400E+00, -0.5300E+00, -0.5100E+00, -0.4800E+00,
     + -0.4100E+00, -0.2600E+00,  0.1700E+00, -0.4200E+00, -0.2200E+00,
     + -0.2080E+01,
     + -0.2700E+00, -0.2600E+00, -0.2500E+00, -0.2400E+00, -0.2100E+00,
     + -0.1600E+00, -0.1000E-01,  0.8500E+00,  0.9700E+00,  0.3500E+00,
     + -0.2250E+01,
     +  0.5300E+00,  0.5400E+00,  0.5600E+00,  0.5700E+00,  0.6000E+00,
     +  0.6500E+00,  0.7600E+00,  0.1960E+01,  0.2280E+01,  0.3600E+00,
     + -0.2590E+01,
     +  0.1590E+01,  0.1630E+01,  0.1650E+01,  0.1680E+01,  0.1720E+01,
     +  0.1780E+01,  0.2010E+01,  0.3260E+01,  0.3110E+01,  0.1600E+00,
     + -0.2580E+01,
     +  0.2560E+01,  0.2620E+01,  0.2660E+01,  0.2690E+01,  0.2740E+01,
     +  0.2830E+01,  0.3400E+01,  0.4030E+01,  0.3030E+01, -0.7000E-01,
     + -0.2320E+01,
     +  0.3500E+01,  0.3600E+01,  0.3650E+01,  0.3700E+01,  0.3750E+01,
     +  0.3860E+01,  0.4580E+01,  0.4100E+01,  0.2450E+01,  0.6000E-01,
     + -0.1770E+01,
     +  0.4500E+01,  0.4640E+01,  0.4700E+01,  0.4760E+01,  0.4830E+01,
     +  0.4930E+01,  0.5420E+01,  0.3960E+01,  0.2000E+01,  0.5000E+00,
     + -0.1150E+01,
     +  0.5290E+01,  0.5470E+01,  0.5550E+01,  0.5620E+01,  0.5690E+01,
     +  0.5790E+01,  0.6070E+01,  0.4000E+01,  0.1910E+01,  0.9700E+00,
     + -0.5600E+00
     +/
      data theta_i /
     +  0.2936E+03,  0.2936E+03,  0.2937E+03,  0.2937E+03,  0.2938E+03,
     +  0.2942E+03,  0.2948E+03,  0.2980E+03,  0.3027E+03,  0.3092E+03,
     +  0.3186E+03,
     +  0.2937E+03,  0.2937E+03,  0.2937E+03,  0.2938E+03,  0.2939E+03,
     +  0.2942E+03,  0.2946E+03,  0.2978E+03,  0.3024E+03,  0.3090E+03,
     +  0.3184E+03,
     +  0.2939E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,
     +  0.2941E+03,  0.2944E+03,  0.2976E+03,  0.3023E+03,  0.3089E+03,
     +  0.3182E+03,
     +  0.2940E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,
     +  0.2941E+03,  0.2943E+03,  0.2975E+03,  0.3022E+03,  0.3087E+03,
     +  0.3181E+03,
     +  0.2940E+03,  0.2940E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,
     +  0.2940E+03,  0.2942E+03,  0.2974E+03,  0.3021E+03,  0.3086E+03,
     +  0.3180E+03,
     +  0.2941E+03,  0.2940E+03,  0.2940E+03,  0.2940E+03,  0.2941E+03,
     +  0.2941E+03,  0.2941E+03,  0.2973E+03,  0.3019E+03,  0.3085E+03,
     +  0.3179E+03,
     +  0.2942E+03,  0.2942E+03,  0.2942E+03,  0.2942E+03,  0.2941E+03,
     +  0.2941E+03,  0.2941E+03,  0.2973E+03,  0.3020E+03,  0.3086E+03,
     +  0.3179E+03,
     +  0.2943E+03,  0.2943E+03,  0.2943E+03,  0.2943E+03,  0.2943E+03,
     +  0.2943E+03,  0.2943E+03,  0.2975E+03,  0.3022E+03,  0.3087E+03,
     +  0.3181E+03,
     +  0.2945E+03,  0.2945E+03,  0.2945E+03,  0.2945E+03,  0.2945E+03,
     +  0.2944E+03,  0.2946E+03,  0.2978E+03,  0.3025E+03,  0.3090E+03,
     +  0.3184E+03,
     +  0.2947E+03,  0.2947E+03,  0.2947E+03,  0.2947E+03,  0.2946E+03,
     +  0.2946E+03,  0.2949E+03,  0.2980E+03,  0.3027E+03,  0.3093E+03,
     +  0.3187E+03,
     +  0.2949E+03,  0.2949E+03,  0.2949E+03,  0.2948E+03,  0.2948E+03,
     +  0.2948E+03,  0.2950E+03,  0.2982E+03,  0.3028E+03,  0.3094E+03,
     +  0.3188E+03,
     +  0.2950E+03,  0.2950E+03,  0.2950E+03,  0.2950E+03,  0.2950E+03,
     +  0.2950E+03,  0.2950E+03,  0.2982E+03,  0.3029E+03,  0.3095E+03,
     +  0.3188E+03
     +/
c
      save time_g, zg, u_i, v_i, theta_i
c
c --------- find time location of initial profiles 
c
      call lterp(nt,time_g,t_factor,jt,jtp1,t_weit)
c
      do iz=1,nz
         ui_temp(iz) = u_i(iz,jt)*(1.0 - t_weit) +
     +                 u_i(iz,jtp1)*t_weit
         vi_temp(iz) = v_i(iz,jt)*(1.0 - t_weit) +
     +                 v_i(iz,jtp1)*t_weit
         ti_temp(iz) = theta_i(iz,jt)*(1.0 - t_weit) +
     +                 theta_i(iz,jtp1)*t_weit
      enddo
c
c ----------- interpolate vertically
c
      do iz=izs,ize
         call lterp(nz,zg,zz(iz),kk,kkp1,weit)
         u_temp = ui_temp(kk)*(1.0 - weit) +
     +            ui_temp(kkp1)*weit
         v_temp = vi_temp(kk)*(1.0 - weit) +
     +            vi_temp(kkp1)*weit
         theta_temp = ti_temp(kk)*(1.0 - weit) +
     +            ti_temp(kkp1)*weit
c
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz)   = u_temp
            v(ix,iy,iz)   = v_temp
            t(ix,iy,1,iz) = theta_temp
            w(ix,iy,iz)   = 0.
            r1(ix,iy,iz)  = 0.
            r2(ix,iy,iz)  = 0.
            r3(ix,iy,iz)  = 0.
            r4(ix,iy,1,iz)= 0.
            r5(ix,iy,iz)  = 0.
         enddo 
         enddo 
      enddo
c
c ------------- set initial random field to be
c               divergence free
c
      idum = -1
      do iz=izs,ize
         if (iz.le.8) then
c
c ----------- ampv and ampt are max amplitudes of random 
c             velocity and temperature fields
c
         ampv = 0.5
         ampt = 0.1
c  
c ------- simple random field scaled between 0 and 1
c
         do iy=iys,iye
         do ix=1,nnx
            psi(ix,iy) = ran1(idum)
         enddo
         enddo
c
         do iy=iys,iye
         do ix=1,nnx
            psix(ix,iy) = psi(ix,iy)
            psiy(ix,iy,izs) = psi(ix,iy)
         enddo
         enddo
         call xderivp(psix(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(psiy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
c
         vmaxx = 0.0
         do iy=iys,iye
         do ix=1,nnx
            vmag = sqrt(psix(ix,iy)**2 + psiy(ix,iy,izs)**2)
            if(vmag .gt. vmaxx) vmaxx = vmag
         enddo
         enddo
         facv = ampv/vmaxx
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz)   = u(ix,iy,iz) - psiy(ix,iy,izs)*facv
            v(ix,iy,iz)   = v(ix,iy,iz) + psix(ix,iy)*facv
            t(ix,iy,1,iz) = t(ix,iy,1,iz) + psi(ix,iy)*ampt
            e(ix,iy,iz)   = 1.0
         enddo
         enddo
         endif
c
c ---------- check divergence of initial field
c
         do iy=iys,iye
         do ix=1,nnx
            uxx(ix,iy)     = u(ix,iy,iz)
            vyy(ix,iy,izs) = v(ix,iy,iz)
         enddo
         enddo
         call xderivp(uxx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(vyy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
         do iy=iys,iye
         do ix=1,nnx
            divz(iz) = divz(iz) + (uxx(ix,iy) + vyy(ix,iy,izs))**2
         enddo
         enddo
         divz(iz) = divz(iz)*fnxy
c
c -------- end z loop
c
      enddo
c
      call mpi_sum_z(divz(1),i_root,myid,nnz,1)
c
      write(nprt,6000)
 6000 format(' check of divergence for initial state',/,
     +         ' iz ',5x,' divergence')
      write(nprt,6100) (iz,divz(iz),iz=izs,ize)
 6100 format(i5,e15.6)
c
c ------------ fix for baroclinic and subsidence effects !!
c
c     do iz=izs,ize
c        ug(iz)=ugcont
c        vg(iz)=vgcont
c        if (.not.(ibrcl.eq.1)) go to 19988
c        if (.not.(iz.le.izi)) go to 19987
c        ug(iz)=0.
c        vg(iz)=0.
c 19987    continue
c 19988    continue
c        zz2=zz(iz)
c        wls(iz)=-divgls*zz2
c        if (.not.(iz.eq.1)) go to 19986
c        do ix=1,nnx
c        uls(ix)=divgls*(dx*float(ix-1)-xl*.5)
c        enddo
c     enddo
c     write(nprt,9)(uls(ix),ix=1,nnx)
c  9  format(1x,8e12.3)
c 19986 continue
c
      return
      end

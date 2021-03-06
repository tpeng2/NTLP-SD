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
c
      z_lower = zi - 50.0
      t_lower = 300.0
      z_upper = zi + 50.0
      t_upper = 308.0
      slope   = (t_upper - t_lower)/(z_upper - z_lower)
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz) = ugcont-ugal
c           u(ix,iy,iz) = ug(iz) - ugal
            v(ix,iy,iz) = vgcont
            w(ix,iy,iz) = 0.0
            e(ix,iy,iz) = 0.0
         enddo
         enddo
         if(z(iz) .le. z_lower) then
           do iy=iys,iye
           do ix=1,nnx
              t(ix,iy,1,iz) = t_lower

           !Add "bubble" to test cloud model
!           t(ix,iy,1,iz) = t(ix,iy,1,iz) + 0.1*
!     +     exp( -( (dx*(ix-1)-500.0)/200.0 )**2 
!     +          -( (dy*(iy-1)-500.0)/200.0 )**2 )*
!     +     exp (-zz(iz)/100.0)

           enddo
           enddo
         elseif(z(iz) .ge. z_upper) then
           do iy=iys,iye
           do ix=1,nnx
              t(ix,iy,1,iz) = t_upper + (zz(iz+1) - z_upper)*dtdzf(1)
           enddo
           enddo
         else
           do iy=iys,iye
           do ix=1,nnx
              t(ix,iy,1,iz) = t_lower + slope*(zz(iz+1) - z_lower)
           enddo
           enddo
         endif
         do iy=iys,iye
         do ix=1,nnx

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


            w(ix,iy,iz)   = 0.
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

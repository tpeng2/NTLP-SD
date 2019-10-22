      subroutine comp1(istep,it)
c
c ----- 3-order runge-kutta time stepping and monotone scalar fluxes in x,y,z.
c       designed to use mpi in x & y directions.
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      parameter(js = 6, ns = 3, nstat = js + ns*nscl)
      real stat(0:nnz,nstat)
c
c ------ temp arrays to hold rhs from step n-1 and
c        field variables from step n 
c
      real urhs(nnx,iys:iye,izs:ize), 
     +     vrhs(nnx,iys:iye,izs:ize),
     +     wrhs(nnx,iys:iye,izs:ize),
     +     erhs(nnx,iys:iye,izs:ize),
     +     trhs(nnx,iys:iye,nscl,izs:ize)
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            urhs(ix,iy,iz) = u(ix,iy,iz) + dtzeta*r1(ix,iy,iz)
            vrhs(ix,iy,iz) = v(ix,iy,iz) + dtzeta*r2(ix,iy,iz)
            wrhs(ix,iy,iz) = w(ix,iy,iz) + dtzeta*r3(ix,iy,iz)
            erhs(ix,iy,iz) = e(ix,iy,iz) + dtzeta*r5(ix,iy,iz)
         enddo
         enddo
      enddo
      do iz=izs,ize
         do l=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            trhs(ix,iy,l,iz) = t(ix,iy,l,iz) + dtzeta*r4(ix,iy,l,iz)
         enddo
         enddo
         enddo
      enddo
c
c --------- get viscosity and rhs of (e,u,v,w)-equations
c           at next step
c
      if (iDNS .eq. 1) then
         call dns_vis
         call rhs_uvw_DNS(istep)
      else
         call tke_vis(istep)   ! WARNING: this has not been checked for iscl=2
         call rhs_uvw(istep)   ! WARNING: this has not been checked for iscl=2
      end if
c
c -------- evaluate rhs of scalar equations
c
      do l=1,nscl
         if (iDNS .eq. 1) then
            call rhs_scl_dns(istep,l,it) 
         else
            call rhs_scl(istep,l,it) 
         end if
      enddo
c
c ---------- gather stat sums on root processor
c            using mpi_reduction over all processors
c
      if(istep .eq. 1) then
c
        do j=1,nstat
        do iz=0,nnz
           stat(iz,j) = 0.0
        enddo
        enddo
        do iz=izs,ize
           stat(iz,1) = uwsb(iz)
           stat(iz,2) = vwsb(iz)
           stat(iz,3) = wwsb(iz)
!           stat(iz,4) = tr_tau(iz)
           stat(iz,5) = triz(iz)
           stat(iz,6) = shrz(iz)
           !stat(iz,7) = t_diss(iz)
        enddo
        m1 = js
        m2 = js + nscl
        m3 = js + 2*nscl
        do l=1,nscl
           do iz=izs,ize
              stat(iz,m1+l) = utsb(iz,l)
              stat(iz,m2+l) = vtsb(iz,l)
              stat(iz,m3+l) = wtsb(iz,l)
           enddo
        enddo

        !Populate the iz=0 locations if processor lies at the bottom -- only for certain quantities
        if (iss .eq. 0) then
           stat(0,1) = uwsb(0)
           stat(0,2) = vwsb(0)
           do l=1,nscl
              stat(0,m3+l) = wtsb(0,l)
           enddo
        end if

        call mpi_sum_z(stat(0,1),i_root,myid,nstat*(nnz+1),1)
        do iz=0,nnz
           uwsb(iz)   = stat(iz,1)
           vwsb(iz)   = stat(iz,2)
        enddo
        do iz=1,nnz
           wwsb(iz)   = stat(iz,3)
!           tr_tau(iz) = stat(iz,4)
           triz(iz)   = stat(iz,5)
           shrz(iz)   = stat(iz,6)
           !t_diss(iz) = stat(iz,7)
        enddo
        do l=1,nscl
           do iz=1,nnz
              utsb(iz,l) = stat(iz,m1+l)
              vtsb(iz,l) = stat(iz,m2+l)
           enddo
           do iz=0,nnz
              wtsb(iz,l) = stat(iz,m3+l)
           enddo
        enddo
        do iz=1,nnz
           buyz(iz) = batag*wtsb(iz,1)
        enddo

c
c -------- end if block
c
      endif
c
c ------- save old rhs in field variables for RK-advancement
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz) = urhs(ix,iy,iz)
            v(ix,iy,iz) = vrhs(ix,iy,iz)
            w(ix,iy,iz) = wrhs(ix,iy,iz)
            e(ix,iy,iz) = erhs(ix,iy,iz)
         enddo
         enddo
      enddo
      do iz=izs,ize
         do l=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            t(ix,iy,l,iz) = trhs(ix,iy,l,iz)
         enddo
         enddo
         enddo
      enddo
c
      return
      end
      subroutine comp2
c
c ------- add p gradients to rhs. Use already defined p
c         at ize+1 to get w (see sr. pressure).
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real fnt1(nnx,iys:iye,izs:ize), fnt2(nnx,iys:iye)
      real r3_sum(1:nnz)
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      do iz=1,nnz
         r3_sum(iz) = 0.0
      enddo
c
c --------- dp/dy at all z
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            fnt1(ix,iy,iz) = p(ix,iy,iz)
         enddo
         enddo
      enddo
c
      call yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
      do iz=izs,ize
c
         izm1  = iz - 1
         izp1  = iz + 1
c
         do iy=iys,iye
         do ix=1,nnx
            fnt2(ix,iy) = p(ix,iy,iz)
         enddo
         enddo
         call xderivp(fnt2(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         do iy=iys,iye
         do ix=1,nnx
            r1(ix,iy,iz) = r1(ix,iy,iz) - fnt2(ix,iy)
            r2(ix,iy,iz) = r2(ix,iy,iz) - fnt1(ix,iy,iz)
         enddo
         enddo
         if (iz.ne.nnz) then
            do iy=iys,iye
            do ix=1,nnx
               r3(ix,iy,iz) = r3(ix,iy,iz) -
     +            (p(ix,iy,izp1)-p(ix,iy,iz))*dzu_i(izp1)
               r3_sum(iz) = r3_sum(iz) + r3(ix,iy,iz)
            enddo
            enddo
            r3_sum(iz) = r3_sum(iz)*fnxy
         endif
c
c ------------------------ time stepping with 3-order rk method
c                          first w variables
c
      if(iz .ne. nnz) then
         do iy=iys,iye
         do ix=1,nnx
c           w(ix,iy,iz)  = w(ix,iy,iz)+dtgama*r3(ix,iy,iz)
            e(ix,iy,iz)  = e(ix,iy,iz)+dtgama*r5(ix,iy,iz)
         enddo
         enddo
      else
c
c --------- update wout and eout by setting = to bc values
c
         do iy=iys,iye
         do ix=1,nnx
            w(ix,iy,iz)  = wbc(ix,iy,1)
            e(ix,iy,iz)  = ebc(ix,iy,1)
            r3(ix,iy,iz) = 0.0
            r5(ix,iy,iz) = 0.0
         enddo
         enddo
      endif
c
c -------- now all u-variables
c
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz) = u(ix,iy,iz)+dtgama*r1(ix,iy,iz)
            v(ix,iy,iz) = v(ix,iy,iz)+dtgama*r2(ix,iy,iz)
         enddo
         enddo
         do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            t(ix,iy,iscl,iz)  = t(ix,iy,iscl,iz)+
     +                          dtgama*r4(ix,iy,iscl,iz)
         enddo
         enddo
         enddo
c
c -------- end z loop
c
      enddo
c
c ---------- gather partial sums for w computation
c
      call mpi_sum_z(r3_sum,i_root,myid,nnz,1)
c
      do iz=izs,min(ize,nnz-1)
         do iy=iys,iye
         do ix=1,nnx
            r3(ix,iy,iz) = r3(ix,iy,iz) - r3_sum(iz)
            w(ix,iy,iz)  = w(ix,iy,iz) + dtgama*r3(ix,iy,iz)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine comp_p
c
c --------- setup pressure solver
c
      use pars
      use fftwk
      use fields
      use con_data
      use con_stats
      include 'mpif.h'
      real fnt1(nnx,iys:iye,izs:ize)
      real fs(nnx,iys:iye,2), fr(nnx,iys:iye,2)
      integer istatus(mpi_status_size)
c
      gami = 1.0/dtgama
c
      nb = myid - ncpu_s
      nt = myid + ncpu_s
c
c ------------ Send both r3 and updated w (from comp1)
c              to processor above the current myid.
c
      if(iss .eq. 0) then
         nb = mpi_proc_null
      endif
      if(ise .eq. numprocs-1) then
         nt = mpi_proc_null
      endif
      nsend = 2*nnx*(iye + 1 - iys)
      nrecv = nsend
      do iy=iys,iye
      do ix=1,nnx
         fs(ix,iy,1) = r3(ix,iy,ize)
         fs(ix,iy,2) = w(ix,iy,ize)
      enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nt,2,
     +     fr(1,iys,1),nrecv,mpi_real8,nb,2,
     +     mpi_comm_world,istatus,ierr)
      if(iss .ne. 0) then
         do iy=iys,iye
         do ix=1,nnx
            r3(ix,iy,izs-1) = fr(ix,iy,1)
            w(ix,iy,izs-1)  = fr(ix,iy,2)
         enddo
         enddo
      endif
c
c ----------- setup general pressure calculation
c             relies on rhs from step n-1 being included 
c             in velocity-arrays already
c
      do iz=izs,ize
         izm1 = iz -1
         do iy=iys,iye
         do ix=1,nnx
            fnt1(ix,iy,iz) = u(ix,iy,iz)*gami + r1(ix,iy,iz)
         enddo
         enddo
         call xderivp(fnt1(1,iys,iz),trigx(1,1),xk(1),nnx,iys,iye)
c
         if(iz .eq. 1) then
            do iy=iys,iye
            do ix=1,nnx
                p(ix,iy,iz) = fnt1(ix,iy,iz) +  
     +                     ((w(ix,iy,iz) -wbc(ix,iy,2))*gami +
     +                       r3(ix,iy,iz))*dzw_i(iz)

            enddo
            enddo
         else if(iz .eq. nnz) then
            do iy=iys,iye
            do ix=1,nnx
                p(ix,iy,iz) = fnt1(ix,iy,iz) + 
     +                     ((wbc(ix,iy,1) - w(ix,iy,izm1))*gami -
     +                      r3(ix,iy,izm1))*dzw_i(iz)
            enddo
            enddo
         else 
            do iy=iys,iye
            do ix=1,nnx
                p(ix,iy,iz) = fnt1(ix,iy,iz) + 
     +                    ((w(ix,iy,iz)  - w(ix,iy,izm1))*gami +
     +                      r3(ix,iy,iz) - r3(ix,iy,izm1))*dzw_i(iz)
            enddo
            enddo
         endif
c
c --------- end z loop
c
      enddo
c
c ----------- check for radiation boundary condition, all processors
c
      if(ibcu .eq. 1) then
        do iy=iys,iye
        do ix=1,nnx
           ptop(ix,iy,1) = pbc(ix,iy,1)
           ptop(ix,iy,2) = pbc2(ix,iy,1)
        enddo
        enddo
      endif
c
c --------- now y contribution
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            fnt1(ix,iy,iz) = v(ix,iy,iz)*gami + r2(ix,iy,iz)
         enddo
         enddo
      enddo
c
      call yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
             p(ix,iy,iz) = p(ix,iy,iz) + fnt1(ix,iy,iz) 
         enddo
         enddo
      enddo
c
      call pressure
c
      return
      end
      subroutine rhs_scl_dns(istep,iscl)
c
c ------ get right hand side of scalar equation (iscl)
c        monotone scalar fluxes only in z
c        for pencil size (nnx, iys:iye, izs:ize) 
c        care is taken so that if monotone is on then
c        conservative horizontal flux form is used!
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use particles
c
c
      real fnt1(nnx,iys:iye,izs:ize)
      real tx(nnx,iys:iye), ty(nnx,iys:iye,izs:ize)
      real flux_u(nnx,iys:iye), flux_l(nnx,iys:iye)
      real taut3_u(nnx,iys:iye,nscl), taut3_l(nnx,iys:iye,nscl)
      real :: sfc_flx
c
c --------- set sign for ocean simulations that use monotone
c
      sgn = 1.0
      if(iocean .eq. 1) sgn = -1.0
      upwn = 2.0
      if(iupwnd .ne. 1) upwn = 1.0
c
c --------- outer loop over z
c
      do iz=izs,ize
c
      izm2 = iz - 2
      izm1 = iz - 1
      izp1 = iz + 1
      izp2 = iz + 2
      weit  = dzw(iz)/(dzw(iz) + dzw(izp1))
      weit1 = 1.0 - weit
      weit3 = dzw(izm1)/(dzw(iz) + dzw(izm1))
      weit4 = 1.0 - weit3
      dzw2_i = 1.0/(dzw(iz) + dzw(izp1))
      dzw3_i = 2.0*dzw2_i
c
      do iy=iys,iye
      do ix=1,nnx
         tx(ix,iy) = t(ix,iy,iscl,iz)
      enddo
      enddo
      call xderivp(tx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
c
c --------- compute tau_t3 at iz-1 
c

         sfc_flx = 0.0
         do iy=iys,iye
         do ix=1,nnx
            taut3_l(ix,iy,iscl) = -vis_s(ix,iy,iscl,izm1)*
     +              (t(ix,iy,iscl,iz) - t(ix,iy,iscl,izm1))*dzu_i(iz)
         if (iz == 1) then 
             sfc_flx = sfc_flx + taut3_l(ix,iy,iscl)
         end if
         enddo
         enddo


!      Aside - get each wtsfc value:
       if (iz == 1 .and. isfc(iscl) == 1) then
          call mpi_sum_xy(sfc_flx,myid,iss,ise,1)
          wtsfc(iscl) = sfc_flx*fnxy 
       end if
c
c ---------- SGS tau_t1, tau_t3 and resolved u*theta scalar fluxes
c            skew symmetric advective term 0.5(udt/dx + d/dx(ut))
c
      do iy=iys,iye
      do ix=1,nnx
         taut3_u(ix,iy,iscl) = -vis_s(ix,iy,iscl,iz)*
     +   (t(ix,iy,iscl,izp1) - t(ix,iy,iscl,iz))*dzu_i(izp1)
c
         fnt1(ix,iy,iz) = -0.5*((vis_s(ix,iy,iscl,iz)+
     +                  vis_s(ix,iy,iscl,izm1))*
     +                    tx(ix,iy) - upwn*t(ix,iy,iscl,iz)*
     +                      (u(ix,iy,iz)+stokes(iz)))
      enddo
      enddo
      call xderivp(fnt1(1,iys,iz),trigx(1,1),xk(1),nnx,iys,iye)
      do iy=iys,iye
      do ix=1,nnx
         r4(ix,iy,iscl,iz) = - fnt1(ix,iy,iz)
     +           -(taut3_u(ix,iy,iscl)-taut3_l(ix,iy,iscl))*dzw_i(iz)
      enddo
      enddo
c
      if(iupwnd .ne. 1) then
c
c --------- skew symmetric advective form for
c           vertical flux = 0.5(wdt/dz + d/dz(wt))
c
      do iy=iys,iye
      do ix=1,nnx
         theta_u = weit1*t(ix,iy,iscl,iz) +
     +                weit*t(ix,iy,iscl,izp1)
         theta_l = weit3*t(ix,iy,iscl,iz) +
     +                weit4*t(ix,iy,iscl,izm1)
         r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz) 
     +     -0.5*(u(ix,iy,iz)+stokes(iz))*tx(ix,iy)
     +     -0.5*(w(ix,iy,iz)*theta_u - w(ix,iy,izm1)*theta_l)*dzw_i(iz)
c
         r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz)
     +     -0.25*(w(ix,iy,iz)*
     +       (t(ix,iy,iscl,izp1) - t(ix,iy,iscl,iz))*dzu_i(izp1) +
     +            w(ix,iy,izm1)*
     +       (t(ix,iy,iscl,iz) - t(ix,iy,iscl,izm1))*dzu_i(iz))
      enddo
      enddo
c
      else
c
c ----------- z-direction special
c
         if(iz .eq. 1) then
              do iy=iys,iye
              do ix=1,nnx
                 flux_l(ix,iy) = sgn*0.5*w(ix,iy,izm1)*
     +                        (t(ix,iy,iscl,izm1)+t(ix,iy,iscl,iz))
                 flux_u(ix,iy) =
     +           amax1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,iscl,iz) +
     +           rlim(t(ix,iy,iscl,izp1),t(ix,iy,iscl,iz),
     +                t(ix,iy,iscl,izm1))) +
     +           amin1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,iscl,izp1) +
     +           rlim(t(ix,iy,iscl,iz),t(ix,iy,iscl,izp1),
     +                t(ix,iy,iscl,izp2)))
              enddo
              enddo
         else if(iz .eq. nnz) then
              do iy=iys,iye
              do ix=1,nnx
                 flux_u(ix,iy) = sgn*0.5*w(ix,iy,iz)*
     +                        (t(ix,iy,iscl,izp1)+t(ix,iy,iscl,iz))
                 flux_l(ix,iy) =
     +           amax1(sgn*w(ix,iy,izm1),0.)*(t(ix,iy,iscl,izm1) +
     +           rlim(t(ix,iy,iscl,iz),t(ix,iy,iscl,izm1),
     +                t(ix,iy,iscl,izm2))) +
     +           amin1(sgn*w(ix,iy,izm1),0.)*(t(ix,iy,iscl,iz) +
     +           rlim(t(ix,iy,iscl,izm1),t(ix,iy,iscl,iz),
     +                t(ix,iy,iscl,izp1)))
              enddo
              enddo
         else
              do iy=iys,iye
              do ix=1,nnx
                 flux_u(ix,iy) =
     +           amax1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,iscl,iz) +
     +           rlim(t(ix,iy,iscl,izp1),t(ix,iy,iscl,iz),
     +                t(ix,iy,iscl,izm1))) +
     +           amin1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,iscl,izp1) +
     +           rlim(t(ix,iy,iscl,iz),t(ix,iy,iscl,izp1),
     +                t(ix,iy,iscl,izp2)))
                 flux_l(ix,iy) =
     +           amax1(sgn*w(ix,iy,izm1),0.)*(t(ix,iy,iscl,izm1) +
     +           rlim(t(ix,iy,iscl,iz),t(ix,iy,iscl,izm1),
     +                t(ix,iy,iscl,izm2))) +
     +           amin1(sgn*w(ix,iy,izm1),0.)*(t(ix,iy,iscl,iz) +
     +           rlim(t(ix,iy,iscl,izm1),t(ix,iy,iscl,iz),
     +                t(ix,iy,iscl,izp1)))
              enddo
              enddo
         endif
c
c ---------- sum vertical monotone flux
c
         do iy=iys,iye
         do ix=1,nnx
            r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz)
     +          - sgn*(flux_u(ix,iy) - flux_l(ix,iy))*dzw_i(iz)
         enddo
         enddo
c
c -------- end monotone if block
c
      endif
c
c -------- save SGS fluxes for printout, gather sums on exit
c
      if(istep .eq. 1) then
         utsb(iz,iscl) = 0.0
         wtsb(iz,iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            wtsb(iz,iscl) = wtsb(iz,iscl) + taut3_u(ix,iy,iscl)
            utsb(iz,iscl) = utsb(iz,iscl) -
     +            0.5*(vis_s(ix,iy,iscl,iz)+
     +                     vis_s(ix,iy,iscl,izm1))*tx(ix,iy)
         enddo
         enddo
         utsb(iz,iscl) = utsb(iz,iscl)*fnxy
         wtsb(iz,iscl) = wtsb(iz,iscl)*fnxy

         !Do it special for wtsb the lower surface:
         if (iz==1) then
         wtsb(0,iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            wtsb(0,iscl) = wtsb(0,iscl) + taut3_l(ix,iy,iscl)
         enddo
         enddo
         wtsb(0,iscl) = wtsb(0,iscl)*fnxy
         end if

      endif
c
c ---------- end z loop
c
      enddo
c
c --------- outer loop over z for y-depenence
c
      do iz=izs,ize
      do iy=iys,iye
      do ix=1,nnx
         ty(ix,iy,iz)  = t(ix,iy,iscl,iz)
      enddo
      enddo
      enddo
c
c --------- y-derivative of t for [izs:ize]
c
      call yd_mpi(ty(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
c ------------- add skew symmetric advective flux and SGS flux
c               to y-derivative computation. check for monotone
c
      do iz=izs,ize
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
            fnt1(ix,iy,iz) = -0.5*((vis_s(ix,iy,iscl,iz)+
     +              vis_s(ix,iy,iscl,izm1))*ty(ix,iy,iz) -
     +                  upwn*t(ix,iy,iscl,iz)*v(ix,iy,iz))
         enddo
         enddo
         if(iupwnd .ne. 1) then
           do iy=iys,iye
           do ix=1,nnx
              r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz) - 
     +                            0.5*v(ix,iy,iz)*ty(ix,iy,iz)
           enddo
           enddo
         endif
      enddo

!---------add on the thermal coupling from the particles:
      if (iscl == 1) then 
      if (iTcouple == 1) then
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
          r4(ix,iy,1,iz) = r4(ix,iy,1,iz) - partTsrc(ix,iy,iz)
         end do
         end do
      end do
      end if
      end if  

      if (iscl == 2) then 
      if (iHcouple == 1) then
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
          r4(ix,iy,1,iz) = r4(ix,iy,1,iz) - partTEsrc(ix,iy,iz)
          r4(ix,iy,2,iz) = r4(ix,iy,2,iz) - partHsrc(ix,iy,iz)
         end do
         end do
      end do
      end if
      end if 



c
c --------- y-derivatives of scalar fluxes for [izs:ize]
c
      call yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),
     +            nnx,nny,ixs,ixe,ix_s,ix_e,
     +            iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz) - fnt1(ix,iy,iz)
         enddo
         enddo
      enddo
c
c -------- save SGS fluxes for printout
c
      if(istep .eq. 1) then
      do iz=izs,ize
         vtsb(iz,iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            vtsb(iz,iscl) = vtsb(iz,iscl) -
     +            0.5*(vis_s(ix,iy,iscl,iz)+
     +                 vis_s(ix,iy,iscl,izm1))*ty(ix,iy,iz)
         enddo
         enddo
         vtsb(iz,iscl) = vtsb(iz,iscl)*fnxy
      enddo
      endif
c
      return
      end subroutine rhs_scl_dns
      subroutine dns_vis
      use particles
      use pars
      use fields
      implicit none


!     In DNS mode, just set the molecular viscosity (and scalar diffusivities)
!     Also, to make the rest of code work, set the rhs of e equation to 0

      !Both for air at the moment:
      vis_m = nuf 

      !Use Prantdl number for thermal diffusivity:
      vis_s(1:nnx,iys:iye,1,izs-1:ize+1) = nuf/Pra   ! alpha=nu/Prandtl   
      vis_s(1:nnx,iys:iye,2,izs-1:ize+1) = nuf/Sc    !Dv =nu/Sc
      r5 = 0.0
      e = 0.0

      end
      subroutine rhs_uvw_DNS(istep)
c
c ---------- get right hand sides of (u,v,w) equations
c            for pencil size (nnx, iys:iye, izs:ize) 
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use particles
c
      real fntd(nnx,iys:iye,izs:ize)
      real fnt1(nnx,iys:iye), fnt2(nnx,iys:iye) 
      real fnt3(nnx,iys:iye), fnt4(nnx,iys:iye)
      real tau13_u(nnx,iys:iye), tau23_u(nnx,iys:iye)
      real tau13_l(nnx,iys:iye), tau23_l(nnx,iys:iye)
      real r3_sum(1:nnz)
      real sfc_flx(2)
c
      do iz=izs,ize
c
      izm1 = iz - 1
      izp1 = iz + 1
      weit  = dzw(iz)/(dzw(iz) + dzw(izp1))
      weit1 = 1.0 - weit
c
c ---------- dynamics 
c
      do iy=iys,iye
      do ix=1,nnx
         uzm = (u(ix,iy,iz)-u(ix,iy,izm1))*dzu_i(iz)
         vzm = (v(ix,iy,iz)-v(ix,iy,izm1))*dzu_i(iz)
         uz  = (u(ix,iy,izp1)-u(ix,iy,iz))*dzu_i(izp1)
         vz  = (v(ix,iy,izp1)-v(ix,iy,iz))*dzu_i(izp1)
c
         u_avg = u(ix,iy,iz)*weit1 + u(ix,iy,izp1)*weit
         v_avg = v(ix,iy,iz)*weit1 + v(ix,iy,izp1)*weit
c
c ------------ advection
c
         u_adv =  v(ix,iy,iz)*(vx(ix,iy,iz)-uy(ix,iy,iz))-
     +          0.5*(w(ix,iy,iz  )*(uz - wx(ix,iy,iz))+
     +           w(ix,iy,izm1)*(uzm - wx(ix,iy,izm1)))
         v_adv = -u(ix,iy,iz)*(vx(ix,iy,iz)-uy(ix,iy,iz))+
     +          0.5*(w(ix,iy,iz  )*(wy(ix,iy,iz) - vz)+
     +           w(ix,iy,izm1)*(wy(ix,iy,izm1) - vzm))
         w_adv = u_avg*(uz - wx(ix,iy,iz))
     +           - v_avg*(wy(ix,iy,iz) - vz)
c
c ------------ coriolis, vertical and horizontal components
c
         u_cor =  fcor*v(ix,iy,iz) - fcor_h*w(ix,iy,iz)
         v_cor = -fcor*(u(ix,iy,iz) + stokes(iz))
         w_cor =  fcor_h*u(ix,iy,iz)
c
c ------------ buoyancy (with hydrostatic part)
c
         w_buy = batag*(t(ix,iy,1,iz)*weit1 +
     +                  t(ix,iy,1,izp1)*weit)
c
c ------------ geostrophic wind
c
         !u_geo = -fcor*vg(iz)
         !v_geo =  fcor*(ug(iz)-ugal)
         !Instead of geostrophic wind (which is a pressure gradient)
         !make u_geo and v_geo equal to my pressure gradient:

         u_geo = -dpdx
         v_geo = 0.0

c
c ------------ totals
c
         r1(ix,iy,iz) = u_adv + u_cor + u_geo
         r2(ix,iy,iz) = v_adv + v_cor + v_geo
         r3(ix,iy,iz) = w_adv + w_cor + w_buy


         !Add particle momentum coupling
         if (icouple == 1) then
         r1(ix,iy,iz) = r1(ix,iy,iz) - partsrc(ix,iy,iz,1)
         r2(ix,iy,iz) = r2(ix,iy,iz) - partsrc(ix,iy,iz,2)
         !Note: partsrc(3,ix,iy,iz) is located at u,v-locations
         !Interpolate to w-location:
         r3(ix,iy,iz) = r3(ix,iy,iz) - (weit*partsrc(ix,iy,izp1,3)+
     +                  weit1*partsrc(ix,iy,iz,3))
         end if

      enddo
      enddo
c
c ---------------- stokes term in ocean cases
c
      if(iocean .eq. 1) then
        stokavg = stokes(iz)*weit1 + stokes(izp1)*weit
        do iy=iys,iye
        do ix=1,nnx
            r2(ix,iy,iz) = r2(ix,iy,iz) + stokes(iz)*
     +                    (uy(ix,iy,iz) - vx(ix,iy,iz))
            uz = (u(ix,iy,izp1) - u(ix,iy,iz))*dzu_i(izp1)
            r3(ix,iy,iz) = r3(ix,iy,iz) + stokavg* 
     +                    (uz - wx(ix,iy,iz))
        enddo
        enddo
      endif
c
c --------- get tau_13,_23 at iz-1 
c
!      Have it compute t13,t23 like normal, even at bottom
!      REQUIRES ghost points to be set correctly for no-slip (done in lower,upper)
!      Also, get rid of the mean correction for that 2-part model
         sfc_flx = 0.0
         do iy=iys,iye
         do ix=1,nnx
            uzm = (u(ix,iy,iz)-u(ix,iy,izm1))*dzu_i(iz)
            vzm = (v(ix,iy,iz)-v(ix,iy,izm1))*dzu_i(iz)
            tau13_l(ix,iy) = -vis_m(ix,iy,izm1)*(uzm + wx(ix,iy,izm1))
            tau23_l(ix,iy) = -vis_m(ix,iy,izm1)*(vzm + wy(ix,iy,izm1))
            if (iz == 1) then
                sfc_flx(1) = sfc_flx(1) + tau13_l(ix,iy)
                sfc_flx(2) = sfc_flx(2) + tau23_l(ix,iy)
            end if
         enddo
         enddo
!      else
!         do iy=iys,iye
!         do ix=1,nnx
!            tau13_l(ix,iy) = tau13m(ix,iy)
!            tau23_l(ix,iy) = tau23m(ix,iy)
!         enddo
!         enddo
!      endif

!      As an aside, compute the mean sfc_flx to put into uwsfc and vwsfc
       if (iz == 1) then 
          call mpi_sum_xy(sfc_flx,myid,iss,ise,2)
          uwsfc = sfc_flx(1)*fnxy
          vwsfc = sfc_flx(2)*fnxy
       end if
c
c ----------- x and z horizontal SGS fluxes for u, v, w
c             tau_11, tau_12, tau_13, tau_23 at iz
c
      do iy=iys,iye
      do ix=1,nnx
         fnt1(ix,iy) = -(vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*
     +                    ux(ix,iy,iz)
         fnt2(ix,iy) = -.5*(vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*
     +                    (uy(ix,iy,iz)+vx(ix,iy,iz))
         uz = (u(ix,iy,izp1)-u(ix,iy,iz))*dzu_i(izp1)
         vz = (v(ix,iy,izp1)-v(ix,iy,iz))*dzu_i(izp1)
         tau13_u(ix,iy) = -vis_m(ix,iy,iz)*(uz+wx(ix,iy,iz))
         tau23_u(ix,iy) = -vis_m(ix,iy,iz)*(vz + wy(ix,iy,iz))
         fnt3(ix,iy) = tau13_u(ix,iy)
      enddo
      enddo
      call xderivp(fnt1(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
      call xderivp(fnt2(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
      call xderivp(fnt3(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
      do iy=iys,iye
      do ix=1,nnx
         r1(ix,iy,iz) = r1(ix,iy,iz) - fnt1(ix,iy)
     +           -(tau13_u(ix,iy)-tau13_l(ix,iy))*dzw_i(iz)
         r2(ix,iy,iz) = r2(ix,iy,iz) - fnt2(ix,iy)
     +           -(tau23_u(ix,iy)-tau23_l(ix,iy))*dzw_i(iz)
         fnt4(ix,iy) = -(vis_m(ix,iy,izm1)+vis_m(ix,iy,iz))*
     +                (w(ix,iy,iz)-w(ix,iy,izm1))*dzw_i(iz)
         fnt2(ix,iy) = -(vis_m(ix,iy,izp1)+vis_m(ix,iy,iz))*
     +                (w(ix,iy,izp1)-w(ix,iy,iz))*dzw_i(izp1)
         r3(ix,iy,iz) = r3(ix,iy,iz) - fnt3(ix,iy) -
     +                   (fnt2(ix,iy)-fnt4(ix,iy))*dzu_i(izp1)
      enddo
      enddo
c
c -------- save SGS fluxes for printout
!          NOTE: now uwsb is t13_viscous,vwsb is t23_viscous, wwsb is t33_viscous
c
      if(istep .eq. 1) then
         uwsb(iz)   = 0.0
         vwsb(iz)   = 0.0
         wwsb(iz)   = 0.0
!         tr_tau(iz) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            uwsb(iz) = uwsb(iz) + tau13_u(ix,iy)
            vwsb(iz) = vwsb(iz) + tau23_u(ix,iy)
            wwsb(iz) = wwsb(iz) + fnt4(ix,iy)
            ufluc    = (u(ix,iy,izp1) - uxym(izp1))*weit +
     +                 (u(ix,iy,iz) - uxym(iz))*weit1
            vfluc    = (v(ix,iy,izp1) - vxym(izp1))*weit +
     +                 (v(ix,iy,iz) - vxym(iz))*weit1
!            tr_tau(iz) = tr_tau(iz) +
!     +                 tau13_u(ix,iy)*ufluc + tau23_u(ix,iy)*vfluc
         enddo
         enddo
         uwsb(iz)   = uwsb(iz)*fnxy
         vwsb(iz)   = vwsb(iz)*fnxy
         wwsb(iz)   = wwsb(iz)*fnxy
!         tr_tau(iz) = tr_tau(iz)*fnxy

         !Save the surface viscous stresses:
         if (iz==1) then
            uwsb(0) = 0.0
            vwsb(0) = 0.0
            do iy=iys,iye
            do ix=1,nnx
               uwsb(0) = uwsb(0) + tau13_l(ix,iy)
               vwsb(0) = vwsb(0) + tau23_l(ix,iy)
            end do
            end do
            uwsb(0) = uwsb(0)*fnxy
            vwsb(0) = vwsb(0)*fnxy
         end if
         

      endif
c
c ---------- end z loop
c
      enddo
c
c ---------- SGS fluxes tau_12, tau_22, tau_23 that depend on 
c            y-derivatives 
c
      do iz=izs,ize
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
            fntd(ix,iy,iz) = -.5*(vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*
     +                   (uy(ix,iy,iz)+vx(ix,iy,iz))
         enddo
         enddo
      enddo
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do iz=izs,ize
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
            r1(ix,iy,iz)   = r1(ix,iy,iz) - fntd(ix,iy,iz)
            fntd(ix,iy,iz) = -(vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*
     +                    vy(ix,iy,iz)
         enddo
         enddo
      enddo
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do iz=izs,ize
         izp1 = iz + 1
         do iy=iys,iye
         do ix=1,nnx
            r2(ix,iy,iz)   = r2(ix,iy,iz) - fntd(ix,iy,iz)
            vz             = (v(ix,iy,izp1)-v(ix,iy,iz))*dzu_i(izp1)
            fntd(ix,iy,iz) = -vis_m(ix,iy,iz)*(vz + wy(ix,iy,iz))
         enddo
         enddo
      enddo
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do iz=1,nnz
         r3_sum(iz) = 0.0
      enddo
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            r3(ix,iy,iz) = r3(ix,iy,iz) - fntd(ix,iy,iz)
            r3_sum(iz)   = r3_sum(iz) + r3(ix,iy,iz)
         enddo
         enddo
         r3_sum(iz) = r3_sum(iz)*fnxy
      enddo
c
      call mpi_sum_z(r3_sum,i_root,myid,nnz,1)
c
c ------- make sure <r3> = 0 and set r3 = 0 at top
c
      do iz=izs,ize
         if(iz .eq. nnz) then
            do iy=iys,iye
            do ix=1,nnx
               r3(ix,iy,iz) = 0.0
            enddo
            enddo
         else
            do iy=iys,iye
            do ix=1,nnx
               r3(ix,iy,iz) = r3(ix,iy,iz) - r3_sum(iz)
            enddo
            enddo
         endif
      enddo
c
      return
      end subroutine rhs_uvw_DNS
      subroutine rhs_scl(istep,iscl)
c
c ------ get right hand side of scalar equation (iscl)
c        monotone scalar fluxes only in z
c        for pencil size (nnx, iys:iye, izs:ize) 
c        care is taken so that if monotone is on then
c        conservative horizontal flux form is used!
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use particles
c
c
      real fnt1(nnx,iys:iye,izs:ize)
      real tx(nnx,iys:iye), ty(nnx,iys:iye,izs:ize)
      real flux_u(nnx,iys:iye), flux_l(nnx,iys:iye)
      real taut3_u(nnx,iys:iye,nscl), taut3_l(nnx,iys:iye,nscl)
c
c --------- set sign for ocean simulations that use monotone
c
      sgn = 1.0
      if(iocean .eq. 1) sgn = -1.0
      upwn = 2.0
      if(iupwnd .ne. 1) upwn = 1.0
c
c --------- outer loop over z
c
      do iz=izs,ize
c
      izm2 = iz - 2
      izm1 = iz - 1
      izp1 = iz + 1
      izp2 = iz + 2
      weit  = dzw(iz)/(dzw(iz) + dzw(izp1))
      weit1 = 1.0 - weit
      weit3 = dzw(izm1)/(dzw(iz) + dzw(izm1))
      weit4 = 1.0 - weit3
      dzw2_i = 1.0/(dzw(iz) + dzw(izp1))
      dzw3_i = 2.0*dzw2_i
c
      do iy=iys,iye
      do ix=1,nnx
         tx(ix,iy) = t(ix,iy,iscl,iz)
      enddo
      enddo
      call xderivp(tx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
c
c --------- compute tau_t3 at iz-1 
c
      if (iz.ne.1 .or. ibcl.ne.0) then
         do iy=iys,iye
         do ix=1,nnx
            taut3_l(ix,iy,iscl) = -vis_s(ix,iy,iscl,izm1)*
     +              (t(ix,iy,iscl,iz) - t(ix,iy,iscl,izm1))*dzu_i(iz)
         enddo
         enddo
      else
         do iy=iys,iye
         do ix=1,nnx
            taut3_l(ix,iy,iscl) = taut3m(ix,iy,iscl)
         enddo
         enddo
      endif
c
c ---------- SGS tau_t1, tau_t3 and resolved u*theta scalar fluxes
c            skew symmetric advective term 0.5(udt/dx + d/dx(ut))
c
      do iy=iys,iye
      do ix=1,nnx
         taut3_u(ix,iy,iscl) = -vis_s(ix,iy,iscl,iz)*
     +      (t(ix,iy,iscl,izp1) - t(ix,iy,iscl,iz))*dzu_i(izp1)
         fnt1(ix,iy,iz) = -0.5*((vis_s(ix,iy,iscl,iz)+
     +           vis_s(ix,iy,iscl,izm1))*tx(ix,iy) - 
     +           upwn*t(ix,iy,iscl,iz)*
     +                      (u(ix,iy,iz)+stokes(iz)))
      enddo
      enddo
      call xderivp(fnt1(1,iys,iz),trigx(1,1),xk(1),nnx,iys,iye)
      do iy=iys,iye
      do ix=1,nnx
         r4(ix,iy,iscl,iz) = - fnt1(ix,iy,iz)
     +           -(taut3_u(ix,iy,iscl)-taut3_l(ix,iy,iscl))*dzw_i(iz)
      enddo
      enddo
c
      if(iupwnd .ne. 1) then
c
c --------- skew symmetric advective form for
c           vertical flux = 0.5(wdt/dz + d/dz(wt))
c
      do iy=iys,iye
      do ix=1,nnx
         theta_u = weit1*t(ix,iy,iscl,iz) +
     +                weit*t(ix,iy,iscl,izp1)
         theta_l = weit3*t(ix,iy,iscl,iz) +
     +                weit4*t(ix,iy,iscl,izm1)
         r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz) 
     +     -0.5*(u(ix,iy,iz)+stokes(iz))*tx(ix,iy)
     +     -0.5*(w(ix,iy,iz)*theta_u - w(ix,iy,izm1)*theta_l)*dzw_i(iz)
c
         r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz)
     +     -0.25*(w(ix,iy,iz)*
     +       (t(ix,iy,iscl,izp1) - t(ix,iy,iscl,iz))*dzu_i(izp1) +
     +            w(ix,iy,izm1)*
     +       (t(ix,iy,iscl,iz) - t(ix,iy,iscl,izm1))*dzu_i(iz))
      enddo
      enddo
c
      else
c
c ----------- z-direction special
c
         if(iz .eq. 1) then
              do iy=iys,iye
              do ix=1,nnx
                 flux_l(ix,iy) = sgn*0.5*w(ix,iy,izm1)*
     +                        (t(ix,iy,iscl,izm1)+t(ix,iy,iscl,iz))
                 flux_u(ix,iy) =
     +           amax1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,iscl,iz) +
     +           rlim(t(ix,iy,iscl,izp1),t(ix,iy,iscl,iz),
     +                t(ix,iy,iscl,izm1))) +
     +           amin1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,iscl,izp1) +
     +           rlim(t(ix,iy,iscl,iz),t(ix,iy,iscl,izp1),
     +                t(ix,iy,iscl,izp2)))
              enddo
              enddo
         else if(iz .eq. nnz) then
              do iy=iys,iye
              do ix=1,nnx
                 flux_u(ix,iy) = sgn*0.5*w(ix,iy,iz)*
     +                        (t(ix,iy,iscl,izp1)+t(ix,iy,iscl,iz))
                 flux_l(ix,iy) =
     +           amax1(sgn*w(ix,iy,izm1),0.)*(t(ix,iy,iscl,izm1) +
     +           rlim(t(ix,iy,iscl,iz),t(ix,iy,iscl,izm1),
     +                t(ix,iy,iscl,izm2))) +
     +           amin1(sgn*w(ix,iy,izm1),0.)*(t(ix,iy,iscl,iz) +
     +           rlim(t(ix,iy,iscl,izm1),t(ix,iy,iscl,iz),
     +                t(ix,iy,iscl,izp1)))
              enddo
              enddo
         else
              do iy=iys,iye
              do ix=1,nnx
                 flux_u(ix,iy) =
     +           amax1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,iscl,iz) +
     +           rlim(t(ix,iy,iscl,izp1),t(ix,iy,iscl,iz),
     +                t(ix,iy,iscl,izm1))) +
     +           amin1(sgn*w(ix,iy,iz),0.)*(t(ix,iy,iscl,izp1) +
     +           rlim(t(ix,iy,iscl,iz),t(ix,iy,iscl,izp1),
     +                t(ix,iy,iscl,izp2)))
                 flux_l(ix,iy) =
     +           amax1(sgn*w(ix,iy,izm1),0.)*(t(ix,iy,iscl,izm1) +
     +           rlim(t(ix,iy,iscl,iz),t(ix,iy,iscl,izm1),
     +                t(ix,iy,iscl,izm2))) +
     +           amin1(sgn*w(ix,iy,izm1),0.)*(t(ix,iy,iscl,iz) +
     +           rlim(t(ix,iy,iscl,izm1),t(ix,iy,iscl,iz),
     +                t(ix,iy,iscl,izp1)))
              enddo
              enddo
         endif
c
c ---------- sum vertical monotone flux
c
         do iy=iys,iye
         do ix=1,nnx
            r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz)
     +          - sgn*(flux_u(ix,iy) - flux_l(ix,iy))*dzw_i(iz)
         enddo
         enddo
c
c -------- end monotone if block
c
      endif
c
c -------- save SGS fluxes for printout, gather sums on exit
c
      if(istep .eq. 1) then
         utsb(iz,iscl) = 0.0
         wtsb(iz,iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            wtsb(iz,iscl) = wtsb(iz,iscl) + taut3_u(ix,iy,iscl)
            utsb(iz,iscl) = utsb(iz,iscl) -
     +            0.5*(vis_s(ix,iy,iscl,iz)+
     +                    vis_s(ix,iy,iscl,izm1))*tx(ix,iy)
         enddo
         enddo
         utsb(iz,iscl) = utsb(iz,iscl)*fnxy
         wtsb(iz,iscl) = wtsb(iz,iscl)*fnxy
      endif
c
c ---------- end z loop
c
      enddo
c
c --------- outer loop over z for y-depenence
c
      do iz=izs,ize
      do iy=iys,iye
      do ix=1,nnx
         ty(ix,iy,iz)  = t(ix,iy,iscl,iz)
      enddo
      enddo
      enddo
c
c --------- y-derivative of t for [izs:ize]
c
      call yd_mpi(ty(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
c ------------- add skew symmetric advective flux and SGS flux
c               to y-derivative computation. check for monotone
c
      do iz=izs,ize
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
            fnt1(ix,iy,iz) = -0.5*(
     +           (vis_s(ix,iy,iscl,iz)+vis_s(ix,iy,iscl,izm1))*
     +                 ty(ix,iy,iz) - upwn*t(ix,iy,iscl,iz)*v(ix,iy,iz))
         enddo
         enddo
         if(iupwnd .ne. 1) then
           do iy=iys,iye
           do ix=1,nnx
              r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz) - 
     +                            0.5*v(ix,iy,iz)*ty(ix,iy,iz)
           enddo
           enddo
         endif
      enddo
c
c --------- y-derivatives of scalar fluxes for [izs:ize]
c
      call yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),
     +            nnx,nny,ixs,ixe,ix_s,ix_e,
     +            iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            r4(ix,iy,iscl,iz) = r4(ix,iy,iscl,iz) - fnt1(ix,iy,iz)
         enddo
         enddo
      enddo

!---------add on the thermal coupling from the particles:
      if (iscl == 1) then
      if (iTcouple == 1) then
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
          r4(ix,iy,1,iz) = r4(ix,iy,1,iz) - partTsrc(ix,iy,iz)
         end do
         end do
      end do
      end if
      end if

      if (iscl == 2) then
      if (iHcouple == 1) then
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
          r4(ix,iy,1,iz) = r4(ix,iy,1,iz) - partTEsrc(ix,iy,iz)
          r4(ix,iy,2,iz) = r4(ix,iy,2,iz) - partHsrc(ix,iy,iz)
         end do
         end do
      end do
      end if
      end if

c
c -------- save SGS fluxes for printout
c
      if(istep .eq. 1) then
      do iz=izs,ize
         vtsb(iz,iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            vtsb(iz,iscl) = vtsb(iz,iscl) -
     +            0.5*(vis_s(ix,iy,iscl,iz)+
     +              vis_s(ix,iy,iscl,izm1))*ty(ix,iy,iz)
         enddo
         enddo
         vtsb(iz,iscl) = vtsb(iz,iscl)*fnxy
      enddo
      endif
c
      return
      end
      subroutine rhs_uvw(istep)
c
c ---------- get right hand sides of (u,v,w) equations
c            for pencil size (nnx, iys:iye, izs:ize) 
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use particles
c
      real fntd(nnx,iys:iye,izs:ize)
      real fnt1(nnx,iys:iye), fnt2(nnx,iys:iye) 
      real fnt3(nnx,iys:iye), fnt4(nnx,iys:iye)
      real tau13_u(nnx,iys:iye), tau23_u(nnx,iys:iye)
      real tau13_l(nnx,iys:iye), tau23_l(nnx,iys:iye)
      real r3_sum(1:nnz)
c
      do iz=izs,ize
c
      izm1 = iz - 1
      izp1 = iz + 1
      weit  = dzw(iz)/(dzw(iz) + dzw(izp1))
      weit1 = 1.0 - weit
c
c ---------- dynamics 
c
      do iy=iys,iye
      do ix=1,nnx
         uzm = (u(ix,iy,iz)-u(ix,iy,izm1))*dzu_i(iz)
         vzm = (v(ix,iy,iz)-v(ix,iy,izm1))*dzu_i(iz)
         uz  = (u(ix,iy,izp1)-u(ix,iy,iz))*dzu_i(izp1)
         vz  = (v(ix,iy,izp1)-v(ix,iy,iz))*dzu_i(izp1)
c
         u_avg = u(ix,iy,iz)*weit1 + u(ix,iy,izp1)*weit
         v_avg = v(ix,iy,iz)*weit1 + v(ix,iy,izp1)*weit
c
c ------------ advection
c
         u_adv =  v(ix,iy,iz)*(vx(ix,iy,iz)-uy(ix,iy,iz))-
     +          0.5*(w(ix,iy,iz  )*(uz - wx(ix,iy,iz))+
     +           w(ix,iy,izm1)*(uzm - wx(ix,iy,izm1)))
         v_adv = -u(ix,iy,iz)*(vx(ix,iy,iz)-uy(ix,iy,iz))+
     +          0.5*(w(ix,iy,iz  )*(wy(ix,iy,iz) - vz)+
     +           w(ix,iy,izm1)*(wy(ix,iy,izm1) - vzm))
         w_adv = u_avg*(uz - wx(ix,iy,iz))
     +           - v_avg*(wy(ix,iy,iz) - vz)
c
c ------------ coriolis, vertical and horizontal components
c
         u_cor =  fcor*v(ix,iy,iz) - fcor_h*w(ix,iy,iz)
         v_cor = -fcor*(u(ix,iy,iz) + stokes(iz))
         w_cor =  fcor_h*u(ix,iy,iz)
c
c ------------ buoyancy (with hydrostatic part)
c

!MODIFY BUOYANCY TERM TO INCLUDE HUMIDITY
!         w_buy = batag*(t(ix,iy,1,iz)*weit1 +
!     +                  t(ix,iy,1,izp1)*weit)

         w_buy = batag*(t(ix,iy,1,iz)*weit1 +
     +                  t(ix,iy,1,izp1)*weit)*
     +   (1.0+0.61*(weit1*t(ix,iy,2,iz) + 
     +    weit*t(ix,iy,2,izp1)))


c
c ------------ geostrophic wind
c
         u_geo = -fcor*vg(iz)
         v_geo =  fcor*(ug(iz)-ugal)
c
c ------------ totals
c
         r1(ix,iy,iz) = u_adv + u_cor + u_geo
         r2(ix,iy,iz) = v_adv + v_cor + v_geo
         r3(ix,iy,iz) = w_adv + w_cor + w_buy


         !Add particle momentum coupling
         if (icouple == 1) then
         r1(ix,iy,iz) = r1(ix,iy,iz) - partsrc(ix,iy,iz,1)
         r2(ix,iy,iz) = r2(ix,iy,iz) - partsrc(ix,iy,iz,2)
         !Note: partsrc(3,ix,iy,iz) is located at u,v-locations
         !Interpolate to w-location:
         r3(ix,iy,iz) = r3(ix,iy,iz) - (weit*partsrc(ix,iy,izp1,3)+
     +                  weit1*partsrc(ix,iy,iz,3))
         end if

      enddo
      enddo
c
c ---------------- stokes term in ocean cases
c
      if(iocean .eq. 1) then
        stokavg = stokes(iz)*weit1 + stokes(izp1)*weit
        do iy=iys,iye
        do ix=1,nnx
            r2(ix,iy,iz) = r2(ix,iy,iz) + stokes(iz)*
     +                    (uy(ix,iy,iz) - vx(ix,iy,iz))
            uz = (u(ix,iy,izp1) - u(ix,iy,iz))*dzu_i(izp1)
            r3(ix,iy,iz) = r3(ix,iy,iz) + stokavg* 
     +                    (uz - wx(ix,iy,iz))
        enddo
        enddo
      endif
c
c --------- get tau_13,_23 at iz-1 
c
      if (iz.ne.1 .or. ibcl.ne.0) then
         do iy=iys,iye
         do ix=1,nnx
            uzm = (u(ix,iy,iz)-u(ix,iy,izm1))*dzu_i(iz)
            vzm = (v(ix,iy,iz)-v(ix,iy,izm1))*dzu_i(iz)
            tau13_l(ix,iy) = -vis_m(ix,iy,izm1)*(uzm + wx(ix,iy,izm1)) -
     +             vis_mean(izm1)*(u_mn(iz)-u_mn(izm1))*dzu_i(iz)
            tau23_l(ix,iy) = -vis_m(ix,iy,izm1)*(vzm + wy(ix,iy,izm1)) -
     +             vis_mean(izm1)*(v_mn(iz)-v_mn(izm1))*dzu_i(iz)
         enddo
         enddo
      else
         do iy=iys,iye
         do ix=1,nnx
            tau13_l(ix,iy) = tau13m(ix,iy)
            tau23_l(ix,iy) = tau23m(ix,iy)
         enddo
         enddo
      endif
c
c ----------- x and z horizontal SGS fluxes for u, v, w
c             tau_11, tau_12, tau_13, tau_23 at iz
c
      do iy=iys,iye
      do ix=1,nnx
         fnt1(ix,iy) = -(vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*
     +                    ux(ix,iy,iz)
         fnt2(ix,iy) = -.5*(vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*
     +                    (uy(ix,iy,iz)+vx(ix,iy,iz))
         uz = (u(ix,iy,izp1)-u(ix,iy,iz))*dzu_i(izp1)
         vz = (v(ix,iy,izp1)-v(ix,iy,iz))*dzu_i(izp1)
         tau13_u(ix,iy) = -vis_m(ix,iy,iz)*(uz+wx(ix,iy,iz)) -
     +            vis_mean(iz)*(u_mn(izp1)-u_mn(iz))*dzu_i(izp1)
         tau23_u(ix,iy) = -vis_m(ix,iy,iz)*(vz + wy(ix,iy,iz)) -
     +            vis_mean(iz)*(v_mn(izp1)-v_mn(iz))*dzu_i(izp1)
         fnt3(ix,iy) = tau13_u(ix,iy)
      enddo
      enddo
      call xderivp(fnt1(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
      call xderivp(fnt2(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
      call xderivp(fnt3(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
      do iy=iys,iye
      do ix=1,nnx
         r1(ix,iy,iz) = r1(ix,iy,iz) - fnt1(ix,iy)
     +           -(tau13_u(ix,iy)-tau13_l(ix,iy))*dzw_i(iz)
         r2(ix,iy,iz) = r2(ix,iy,iz) - fnt2(ix,iy)
     +           -(tau23_u(ix,iy)-tau23_l(ix,iy))*dzw_i(iz)
         fnt4(ix,iy) = -(vis_m(ix,iy,izm1)+vis_m(ix,iy,iz))*
     +                (w(ix,iy,iz)-w(ix,iy,izm1))*dzw_i(iz)
         fnt2(ix,iy) = -(vis_m(ix,iy,izp1)+vis_m(ix,iy,iz))*
     +                (w(ix,iy,izp1)-w(ix,iy,iz))*dzw_i(izp1)
         r3(ix,iy,iz) = r3(ix,iy,iz) - fnt3(ix,iy) -
     +                   (fnt2(ix,iy)-fnt4(ix,iy))*dzu_i(izp1)
      enddo
      enddo
c
c -------- save SGS fluxes for printout
c
      if(istep .eq. 1) then
         uwsb(iz)   = 0.0
         vwsb(iz)   = 0.0
         wwsb(iz)   = 0.0
         tr_tau(iz) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            uwsb(iz) = uwsb(iz) + tau13_u(ix,iy)
            vwsb(iz) = vwsb(iz) + tau23_u(ix,iy)
            wwsb(iz) = wwsb(iz) + fnt4(ix,iy)
            ufluc    = (u(ix,iy,izp1) - uxym(izp1))*weit +
     +                 (u(ix,iy,iz) - uxym(iz))*weit1
            vfluc    = (v(ix,iy,izp1) - vxym(izp1))*weit +
     +                 (v(ix,iy,iz) - vxym(iz))*weit1
            tr_tau(iz) = tr_tau(iz) +
     +                 tau13_u(ix,iy)*ufluc + tau23_u(ix,iy)*vfluc
         enddo
         enddo
         uwsb(iz)   = uwsb(iz)*fnxy
         vwsb(iz)   = vwsb(iz)*fnxy
         wwsb(iz)   = wwsb(iz)*fnxy
         tr_tau(iz) = tr_tau(iz)*fnxy
      endif
c
c ---------- end z loop
c
      enddo
c
c ---------- SGS fluxes tau_12, tau_22, tau_23 that depend on 
c            y-derivatives 
c
      do iz=izs,ize
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
            fntd(ix,iy,iz) = -.5*(vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*
     +                   (uy(ix,iy,iz)+vx(ix,iy,iz))
         enddo
         enddo
      enddo
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do iz=izs,ize
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
            r1(ix,iy,iz)   = r1(ix,iy,iz) - fntd(ix,iy,iz)
            fntd(ix,iy,iz) = -(vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*
     +                    vy(ix,iy,iz)
         enddo
         enddo
      enddo
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do iz=izs,ize
         izp1 = iz + 1
         do iy=iys,iye
         do ix=1,nnx
            r2(ix,iy,iz)   = r2(ix,iy,iz) - fntd(ix,iy,iz)
            vz             = (v(ix,iy,izp1)-v(ix,iy,iz))*dzu_i(izp1)
            fntd(ix,iy,iz) = -vis_m(ix,iy,iz)*(vz + wy(ix,iy,iz)) -
     +            vis_mean(iz)*(v_mn(izp1)-v_mn(iz))*dzu_i(izp1)
         enddo
         enddo
      enddo
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do iz=1,nnz
         r3_sum(iz) = 0.0
      enddo
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            r3(ix,iy,iz) = r3(ix,iy,iz) - fntd(ix,iy,iz)
            r3_sum(iz)   = r3_sum(iz) + r3(ix,iy,iz)
         enddo
         enddo
         r3_sum(iz) = r3_sum(iz)*fnxy
      enddo
c
      call mpi_sum_z(r3_sum,i_root,myid,nnz,1)
c
c ------- make sure <r3> = 0 and set r3 = 0 at top
c
      do iz=izs,ize
         if(iz .eq. nnz) then
            do iy=iys,iye
            do ix=1,nnx
               r3(ix,iy,iz) = 0.0
            enddo
            enddo
         else
            do iy=iys,iye
            do ix=1,nnx
               r3(ix,iy,iz) = r3(ix,iy,iz) - r3_sum(iz)
            enddo
            enddo
         endif
      enddo
c
      return
      end
      subroutine tke_budget

!NOTE::: THIS HAS NOT BEEN UPDATED TO REFLECT HUMIDITY IN THE BUOYANCY
!EQUATION   DHR 7/5/16

c
c -------- get terms in resolved scale tke budget
c          as in gabls writeup at w-points
c          at istage = 1. 
c
      use pars
      use particles
      use fields
      use con_data
      use con_stats
c
      real :: stat(1:nnz,5)
      real :: s11s,s22s,s33s,s12s,s13s,s23s,wz,wzp
      real :: s13,s23,s33
      real :: ufluc,ufluc_t,ufluc_b,vfluc,vfluc_t,vfluc_b,wfluc
      real :: ffluc1,ffluc2,ffluc3
      real :: ffluc1p,ffluc2p,ffluc3p
      real :: weit,weit1
      integer :: iz,i,j,izp1,izm1
c
c -------- stat(.,1) = tke transport  = wq
c          stat(.,2) = pressure transport  = wp
c          stat(.,3) = tke dissipation
c          stat(.,4) = tke dissipation
c          stat(.,5) = particle force correlation
c
      do iz=1,nnz
         stat(iz,1) = 0.0
         stat(iz,2) = 0.0
         stat(iz,3) = 0.0
         stat(iz,4) = 0.0
         stat(iz,5) = 0.0
      enddo

!Compute DNS dissipation, since there is no subgrid dissipation now:
      do iz=izs,ize
         izp1 = iz + 1
         izm1 = iz - 1
         weit = dzw(iz)/(dzw(iz) + dzw(izp1))
         weit1 = 1.0 - weit
!
! ---- get fluctuating strain rates:
!      here, sij = 1/2*(duidxj + dujdxi)
!      then t_diss = 2*nu*<sij sij>
! ---- NOTE: these are computed at the w-locations!  (not u,v locations)
!
         t_diss(iz) = 0.0
         do j=iys,iye
         do i=1,nnx

            !Things for dissipation - these are computed at w-locations since
            !there is no z-derivative after
            s11s = weit1*ux(i,j,iz)**2 + weit*ux(i,j,izp1)**2
            s22s = weit1*vy(i,j,iz)**2 + weit*vy(i,j,izp1)**2
            wz  = (w(i,j,iz)-w(i,j,izm1))*dzw_i(iz)
            wzp = (w(i,j,izp1)-w(i,j,iz))*dzw_i(izp1)
            s33s = weit*wzp**2 + weit1*wz**2
            s12s = weit1*(0.5*(uy(i,j,iz) + vx(i,j,iz)))**2 +
     +            weit*(0.5*(uy(i,j,izp1) + vx(i,j,izp1)))**2
            s13s = (0.5*((u(i,j,izp1) - u(i,j,iz) +
     +            uxym(iz) - uxym(izp1))*dzu_i(izp1) +
     +            wx(i,j,iz)))**2
            s23s = (0.5*((v(i,j,izp1) - v(i,j,iz) +
     +          vxym(iz) - vxym(izp1))*dzu_i(izp1) +
     +          wy(i,j,iz)))**2

         stat(iz,3) = stat(iz,3) + 2.0*vis_m(i,j,iz)*(
     +               s11s+s22s+s33s+2.0*(s12s+s13s+s23s))

            !Things for viscous transport - these are computed at u-locations since
            !a z-derivative is done to the average, which lands t_tau on w-locations
            
            ufluc_t   = u(i,j,izp1) - uxym(izp1)
            ufluc   = u(i,j,iz) - uxym(iz)
            ufluc_b   = u(i,j,izm1) - uxym(izm1)
            vfluc_t   = v(i,j,izp1) - vxym(izp1)
            vfluc   = v(i,j,iz) - vxym(iz)
            vfluc_b   = v(i,j,izm1) - vxym(izm1)
            wfluc = 0.5*( (w(i,j,iz)-wxym(iz)) 
     +                  + (w(i,j,izm1)-wxym(izm1)) )

            uz_t = (ufluc_t-ufluc)*dzu_i(izp1)
            uz_b = (ufluc-ufluc_b)*dzu_i(iz)
            vz_t = (vfluc_t-vfluc)*dzu_i(izp1)
            vz_b = (vfluc-vfluc_b)*dzu_i(iz)
            
            uz = 0.5*(uz_t+uz_b)
            vz = 0.5*(vz_t+vz_b)

            s13 = 0.5*(uz + 0.5*(wx(i,j,iz)+wx(i,j,izm1)))
            s23 = 0.5*(vz + 0.5*(wy(i,j,iz)+wy(i,j,izm1)))
            s33 = wz

         !Note: just uses vis_m(1,1,iz) since it's equal everywhere:
         stat(iz,4) = stat(iz,4) + 2.0*vis_m(i,j,iz)*(
     +               ufluc*s13 + vfluc*s23 + wfluc*s33)
                    
         !Finally get the particle force correlation term:
         ffluc1 = partsrc(i,j,iz,1)-m1src(iz)
         ffluc2 = partsrc(i,j,iz,2)-m2src(iz)
         ffluc3 = partsrc(i,j,iz,3)-m3src(iz)
         if (iz==nnz) then
         ffluc1p = 0.0
         ffluc2p = 0.0
         ffluc3p = 0.0
         else
         ffluc1p = partsrc(i,j,izp1,1)-m1src(izp1)
         ffluc2p = partsrc(i,j,izp1,2)-m2src(izp1)
         ffluc3p = partsrc(i,j,izp1,3)-m3src(izp1)
         end if
         stat(iz,5) = stat(iz,5) +
     +       weit*(ufluc_t*ffluc1p)+weit1*(ufluc*ffluc1)+
     +       weit*(vfluc_t*ffluc2p) + weit1*(vfluc*ffluc2) + 
     +       (w(i,j,iz)-wxym(iz))*(weit*ffluc3p+weit1*ffluc3)
         enddo
         enddo
         stat(iz,3) = stat(iz,3)*fnxy
         stat(iz,4) = stat(iz,4)*fnxy
         stat(iz,5) = stat(iz,5)*fnxy
      enddo
c
c --------------- get transport terms as vertical arrays
c
      do iz=izs,ize
c
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
c
c --------- get estimate of turbulent transport term
c
            ufluc   = u(ix,iy,iz) - uxym(iz)
            vfluc   = v(ix,iy,iz) - vxym(iz)
            wfluc   = w(ix,iy,iz) - wxym(iz)
            wfluc_l = w(ix,iy,izm1) - wxym(izm1)
            stat(iz,1)  = stat(iz,1) + 0.25*(wfluc + wfluc_l)*
     +             (ufluc**2 + vfluc**2 + 0.5*(wfluc_l**2 + wfluc**2))
c
c --------- get estimate of pressure transport term
c
            pfluc = p(ix,iy,iz) - pxym(iz)
     +           -0.5*((u(ix,iy,iz)+stokes(iz))**2 + 
     +                 v(ix,iy,iz)*v(ix,iy,iz) + 
     +      0.5*(w(ix,iy,iz)*w(ix,iy,iz)+w(ix,iy,izm1)*w(ix,iy,izm1)))
            stat(iz,2) = stat(iz,2) + pfluc*0.5*(wfluc_l + wfluc)
         enddo
         enddo
         stat(iz,1) = stat(iz,1)*fnxy
         stat(iz,2) = stat(iz,2)*fnxy
      enddo
      call mpi_sum_z(stat(1,1),i_root,myid,nnz*5,1)
c
c ------ we have all terms on all processors for all z, add them up
c        treat tr_tau at bottom special, tr_tau = 0
c
!      tr_tau(0) = 0.0
      do iz=1,nnz
c
         izp1 = iz + 1
         izm1 = iz - 1
         if(iz .eq. nnz) then
            t_tau(iz) = 0.0
            t_wp(iz)  = 0.0
            t_wq(iz)  = 0.0
         else
            t_tau(iz) = (stat(izp1,4) - stat(iz,4))*dzu_i(izp1) 
            t_wq(iz)  = -(stat(izp1,1) - stat(iz,1))*dzu_i(izp1)
            t_wp(iz)  = -(stat(izp1,2) - stat(iz,2))*dzu_i(izp1)
         endif
         dudz = (uxym(izp1) - uxym(iz))*dzu_i(izp1)
         dvdz = (vxym(izp1) - vxym(iz))*dzu_i(izp1)
c
c ------------- gather all the budget terms
c
         t_tran(iz)  = t_wq(iz) + t_wp(iz) + t_tau(iz)
         t_rprod(iz) = -(dudz*uwle(iz) + dvdz*vwle(iz))
         !Old t_sprod had subgrid stuff
         !t_sprod(iz) =  (dudz*uwsb(iz) + dvdz*vwsb(iz))
         !Now make t_sprod the spray tke term to reduce new variables:
         t_sprod(iz) = -stat(iz,5)
         t_buoy(iz)  =  batag*wtle(iz,1)
         t_diss(iz) = stat(iz,3)
c
      enddo
c
      return
      end
      subroutine tke_vis(istep)
c
c -------------- get viscosity using deardorff tke model with
c                stability correction. fixes for surface layer. 
c                 get rhs of e-equation
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      real fnt1(nnx,iys:iye), fnt2(nnx,iys:iye,izs:ize)
      real fnt3(nnx,iys:iye)
      real ex(nnx,iys:iye), ey(nnx,iys:iye,izs:ize)
      real u_avg(nnx,iys:iye), v_avg(nnx,iys:iye), dissp(nnx,iys:iye)
      real alk(nnx,iys:iye,izs-1:ize+1)
c
      do iz=izs-1,ize+1
c
      izp1 = iz + 1
      dslk  = dsl_z(iz)
      if(iz .gt. 0) dslk  = amin1(dsl_z(iz),vk*abs(z(iz))/csmag)
      almin = almin_c*dsl_z(iz)
      if(iz .eq. 0 .or. iz .eq. nnz+1) then
         dfack = 1.0
      else
         dfack = dfac(iz)
      endif
      if(ivis .eq. 1 .and. iz .le. nmatch) then
c
c --------------- no stability corrected length scales
c
         do j=iys,iye
         do i=1,nnx
            alk(i,j,iz) = dslk
         end do
         end do
      else
         do j=iys,iye
         do i=1,nnx
            alk(i,j,iz) = dslk
            stab = batag*(t(i,j,1,izp1) - t(i,j,1,iz))*dzu_i(izp1)
            if(stab.gt.stabmin) then
              als = stab_c*sqrt(e(i,j,iz)/stab)
              alk(i,j,iz) = amin1(dslk,als)
            endif
            alk(i,j,iz)  = amax1(almin,alk(i,j,iz))
         enddo
         enddo
      endif
      do j=iys,iye
      do i=1,nnx
         vis_m(i,j,iz) = ck*alk(i,j,iz)*sqrt(e(i,j,iz))*dfack
         vis_s(i,j,1:nscl,iz) = (1.+2.*alk(i,j,iz)/dslk)*vis_m(i,j,iz) 
      enddo
      enddo
c
c -------------- special case for iz = 1
c
      if(iz.eq.1 .and. ibcl .eq. 0) then
         do iy=iys,iye
         do ix=1,nnx
            vis_m(ix,iy,iz-1) = vis_m(ix,iy,iz)
            vis_s(ix,iy,1:nscl,iz-1) = vis_s(ix,iy,1:nscl,iz)
         enddo
         enddo
      endif
c
c -------------- end z loop
c
      enddo
c
c -------------- if special 2 part surface layer model is on
c                get "mean" viscosity
c
      do iz=izs-1,ize
         izm1         = iz - 1
         izp1         = iz + 1
         vis_mean(iz) = 0.0
         if(ivis .eq. 1 .and. iz .le. nmatch) then
            if(iz .le. 1) then
              vis_mean(iz) = xksurf
            else
              stravg = sqrt((u_mn(izp1)-u_mn(iz))**2 + 
     +              (v_mn(izp1)-v_mn(iz))**2)*abs(dzu_i(izp1))
              vis_mean(iz) = xksurf*viscon*stravg
            endif
         endif
      enddo
c
c --------- update rhs of sgs e from x and z pieces
c           cube of size (nnx, iys,iye, izs:ize)
c
      do iz=izs,ize
c
      izm1   = iz - 1
      izp1   = iz + 1
      weit   = dzw(iz)/(dzw(iz) + dzw(izp1))
      weit1  = 1.0 - weit
      dzw2_i = 1.0/(dzw(iz) + dzw(izp1))
      dzw3_i = 2.0*dzw2_i
      dslk   = dsl_z(iz)
c
      do iy=iys,iye
      do ix=1,nnx
         ex(ix,iy) = e(ix,iy,iz)
      enddo
      enddo
      call xderivp(ex(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
c
c ------------ include stokes contribution in advection
c              and horizontal x-diffusion
c
      do iy=iys,iye
      do ix=1,nnx
         u_avg(ix,iy)   = (stokes(iz) + u(ix,iy,iz))*weit1 +
     +                    (stokes(izp1) + u(ix,iy,izp1))*weit
         fnt1(ix,iy)    = e(ix,iy,iz)*u_avg(ix,iy) - 
     +                    4.0*vis_m(ix,iy,iz)*ex(ix,iy)
      enddo
      enddo
      call xderivp(fnt1(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
      do iy=iys,iye
      do ix=1,nnx
         r5(ix,iy,iz) = -fnt1(ix,iy) - 
     +         (w(ix,iy,izp1)*e(ix,iy,izp1) -
     +          w(ix,iy,izm1)*e(ix,iy,izm1))*dzw2_i
c
	r5(ix,iy,iz)=0.25*((r5(ix,iy,iz) - u_avg(ix,iy)*ex(ix,iy))*2.0
     +        - w(ix,iy,iz)*(e(ix,iy,izp1)-e(ix,iy,izm1))*dzw3_i)
      enddo
      enddo
c
c ------------- 9/1989 add ihflt=1 option--mean shear does not generate sgs tke
c
      uxymm=0.
      uxymp=0.
      vxymm=0.
      vxymp=0.
      if(ivis .eq. 1 .and. iz .le. nmatch) then
         uxymm = u_mn(iz)
         uxymp = u_mn(izp1)
         vxymm = v_mn(iz)
         vxymp = v_mn(izp1)
      endif
c
      do iy=iys,iye
      do ix=1,nnx
c
c ----------------- dissipation 
c
         dissp(ix,iy) =  (0.19+0.74*alk(ix,iy,iz)/dslk)*
     +            e(ix,iy,iz)*sqrt(e(ix,iy,iz))/alk(ix,iy,iz)
         r5(ix,iy,iz)=r5(ix,iy,iz) - dissp(ix,iy)
c
c ----------------- vertical diffusion
c
         fnt3(ix,iy) = 
     +      ((vis_m(ix,iy,izp1)+vis_m(ix,iy,iz))*
     +       (e(ix,iy,izp1)-e(ix,iy,iz))*dzw_i(izp1) -
     +       (vis_m(ix,iy,iz)+vis_m(ix,iy,izm1))*
     +       (e(ix,iy,iz  )-e(ix,iy,izm1))*dzw_i(iz))*dzu_i(izp1)
         r5(ix,iy,iz) = r5(ix,iy,iz) + fnt3(ix,iy)
c
c ----------------- shear production
c
         s11 = weit1*ux(ix,iy,iz)**2 + weit*ux(ix,iy,izp1)**2
         s22 = weit1*vy(ix,iy,iz)**2 + weit*vy(ix,iy,izp1)**2
         wz  = (w(ix,iy,iz)-w(ix,iy,izm1))*dzw_i(iz)
         wzp = (w(ix,iy,izp1)-w(ix,iy,iz))*dzw_i(izp1)
         s33 = weit*wzp**2 + weit1*wz**2
         s12 = weit1*(uy(ix,iy,iz) + vx(ix,iy,iz))**2 +
     +         weit*(uy(ix,iy,izp1) + vx(ix,iy,izp1))**2
         uzmn=(u(ix,iy,izp1)-uxymp-u(ix,iy,iz)+uxymm)*dzu_i(izp1) 
         vzmn=(v(ix,iy,izp1)-vxymp-v(ix,iy,iz)+vxymm)*dzu_i(izp1)
         s13 = (uzmn + wx(ix,iy,iz))**2
         s23 = (vzmn + wy(ix,iy,iz))**2
c
         fnt1(ix,iy) = vis_m(ix,iy,iz)*(2.0*(s11 + s22 + s33) +
     +                                   s13 + s23 + s12)
         r5(ix,iy,iz) = r5(ix,iy,iz) + fnt1(ix,iy)
c
c ----------------- buoyancy, get tau_w*theta
c
         buoy_sgs = -vis_s(ix,iy,1,iz)*(t(ix,iy,1,izp1) -
     +                      t(ix,iy,1,iz))*dzu_i(izp1)
                     !check if iscl is needed
         r5(ix,iy,iz) = r5(ix,iy,iz) + batag*buoy_sgs
c
         enddo
         enddo
c
c ---------------- compute shear, buoyancy, diffusion
c                  terms in SGS e eqn for printout
c            **** triz is only vertical diffusion ****
c
      if(istep .eq. 1) then
         shrz(iz)   = 0.0
         triz(iz)   = 0.0
         t_diss(iz) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            shrz(iz)   = shrz(iz) + fnt1(ix,iy)
            t_diss(iz) = t_diss(iz) + dissp(ix,iy)
            triz(iz)   = triz(iz) + fnt3(ix,iy)
         enddo
         enddo
         shrz(iz)   = shrz(iz)*fnxy
         t_diss(iz) = t_diss(iz)*fnxy
         triz(iz)   = triz(iz)*fnxy
      endif
c
c -------------- end z loop
c
      enddo
c
c --------- update tendency of sgs e from y contributions
c           pencil size (nnx,iys:iye,izs:ize)
c
      do iz=izs,ize
      do iy=iys,iye
      do ix=1,nnx
         ey(ix,iy,iz) = e(ix,iy,iz)
      enddo
      enddo
      enddo
c
      call yd_mpi(ey(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
c ------ skew symmetic advection [vde/dy + d/dy(ve)]/2
c        plus SGS diffusion contribution
c
      do iz=izs,ize
      izm1   = iz - 1
      izp1   = iz + 1
      weit   = dzw(iz)/(dzw(iz) + dzw(izp1))
      weit1  = 1.0 - weit
      do iy=iys,iye
      do ix=1,nnx
         v_avg(ix,iy)   = v(ix,iy,iz)*weit1 + v(ix,iy,izp1)*weit
         fnt2(ix,iy,iz) = e(ix,iy,iz)*v_avg(ix,iy) -
     +                    4.0*vis_m(ix,iy,iz)*ey(ix,iy,iz)
         r5(ix,iy,iz)   = r5(ix,iy,iz) - 0.5*(v_avg(ix,iy)*ey(ix,iy,iz)) 
      enddo
      enddo
      enddo
c
      call yd_mpi(fnt2(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
      do iz=izs,ize
      do iy=iys,iye
      do ix=1,nnx
         r5(ix,iy,iz) = r5(ix,iy,iz) - 0.5*fnt2(ix,iy,iz)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine Tvar_budget(iscl)
      use pars
      use particles
      use fields
      use con_data
      use con_stats
      use fftwk
      implicit none

      real :: stat(1:nnz,4)
      real :: weit,weit1
      real :: tx_tmp(nnx,iys:iye), ty(nnx,iys:iye,izs-1:ize+1)
      real :: tx(nnx,iys:iye,izs-1:ize+1)
      real :: trans(izs:ize)
      real :: gradTp(3),Tfluc,dTpdz1,dTpdz,dTdz
      real :: Tflucp1,Tflucm1,qfluc,qflucp1,Tmean
      integer :: iz,i,j,izp1,izm1,iscl

c -------- stat(.,1) = Transport: -del.[U<T'2> + <u T'2> - alpha*grad(T'2)]
c          stat(.,2) = Dissipation: -2*alpha <grad(T').grad(T')> 
c          stat(.,3) = Particle: <T' Q'>

      !Need the y gradient of temp:
      do iz=izs-1,ize+1
      do j=iys,iye
      do i=1,nnx
         ty(i,j,iz)  = t(i,j,iscl,iz)
         tx_tmp(i,j)  = t(i,j,iscl,iz)
      enddo
      enddo
      call xderivp(tx_tmp(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
      tx(1:nnx,iys:iye,iz) = tx_tmp(1:nnx,iys:iye)
      enddo

      call yd_mpi(ty(1,iys,izs-1),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
       

      stat = 0.0
      Tv_part2 = 0.0
      do iz=izs,ize
         izp1 = iz + 1
         izm1 = iz - 1
         weit = dzw(iz)/(dzw(iz) + dzw(izp1))
         weit1 = 1.0 - weit

       do j=iys,iye
       do i=1,nnx

         if (iz==1)  then
         Tmean = 2.0*Tbot(iscl) - txym(iz,iscl)
         Tflucm1 = t(i,j,iscl,izm1)-Tmean
         Tflucp1 = t(i,j,iscl,izp1)-txym(izp1,iscl)
         elseif (iz==nnz) then
         Tmean = 2.0*Ttop(iscl) - txym(iz,iscl)
         Tflucp1 = t(i,j,iscl,izp1)-Tmean
         Tflucm1 = t(i,j,iscl,izm1)-txym(izm1,iscl)
         else
         Tflucp1 = t(i,j,iscl,izp1)-txym(izp1,iscl)
         Tflucm1 = t(i,j,iscl,izm1)-txym(izm1,iscl)
         end if
         Tfluc = t(i,j,iscl,iz)-txym(iz,iscl)

         !First dissipation: 
         !Note that gradients of total T and T' in x,y directions are equal since d<T>/dx = d<T>/dy = 0
         gradTp(1) = weit1*tx(i,j,iz) + weit*tx(i,j,izp1)
         gradTp(2) = weit1*ty(i,j,iz) + weit*ty(i,j,izp1)

         !Now get dT'/dz:
         gradTp(3) = (Tflucp1 - Tfluc)*dzu_i(izp1)
         
         stat(iz,2) = stat(iz,2)  - 2.0*vis_s(i,j,iscl,iz)*
     +                (gradTp(1)**2+gradTp(2)**2+gradTp(3)**2)

         !Next transport

         !Store the transport sum at u-locations since z-derivative at the end

         stat(iz,1) = stat(iz,1) + w(i,j,iz)*Tfluc**2

         dTpdz1 = (Tflucp1**2-Tfluc**2)*dzu_i(izp1)
         dTpdz = (Tfluc**2-Tflucm1**2)*dzu_i(iz)
         stat(iz,1) = stat(iz,1) - vis_s(i,j,iscl,iz)*0.5*(dTpdz1+dTpdz)

          !Particle source:
          if (iscl == 1) then
          if (iTcouple ==1) then
            qfluc = partTsrc(i,j,iz)-Tpsrc(iz)
            if (iz==nnz) then
            qflucp1 = 0.0
            else
            qflucp1 = partTsrc(i,j,izp1)-Tpsrc(izp1)
            end if
            stat(iz,3) = stat(iz,3) + 
     +         weit*(qflucp1*Tflucp1) + weit1*(qfluc*Tfluc)
          endif
          endif


          if (iscl == 2) then
          if (iHcouple ==1) then
            qfluc = partHsrc(i,j,iz)-Hpsrc(iz)
            if (iz==nnz) then
            qflucp1 = 0.0
            else
            qflucp1 = partHsrc(i,j,izp1)-Hpsrc(izp1)
            end if
            stat(iz,3) = stat(iz,3) + 
     +         weit*(qflucp1*Tflucp1) + weit1*(qfluc*Tfluc)
          endif
          endif

          if (iscl .eq. 2 .and. iHcouple .eq. 1 ) then
              qfluc = partTEsrc(i,j,iz)-TEpsrc(iz)
              if (iz==nnz) then
              qflucp1 = 0.0
              else
              qflucp1 = partTEsrc(i,j,izp1)-TEpsrc(izp1)
              end if
              Tv_part2(iz) = Tv_part2(iz) + 
     +          weit*(qflucp1*Tflucp1) + weit1*(qfluc*Tfluc)
  
         endif
     
       end do
       end do

         stat(iz,1) = stat(iz,1)*fnxy
         stat(iz,2) = stat(iz,2)*fnxy
         stat(iz,3) = stat(iz,3)*fnxy
         if (iscl .eq. 2 .and. iHcouple .eq. 1 ) then
         Tv_part2(iz) = Tv_part2(iz)*fnxy
         endif
      end do
         

      call mpi_sum_z(stat(1,1),i_root,myid,nnz*3,1)
      if (iscl .eq. 2 .and. iHcouple .eq. 1 ) then
      call mpi_sum_z(Tv_part2(1),i_root,myid,nnz,1)
      endif


c
c ------ we have all terms on all processors for all z, add them up
c
      do iz=1,nnz
         izp1 = iz + 1
         izm1 = iz - 1
         if(iz .eq. nnz) then
            Tv_tran(iz,iscl) = 0.0
         else
            Tv_tran(iz,iscl) = -(stat(izp1,1)-stat(iz,1))*dzu_i(izp1)
         endif
         dTdz = (txym(izp1,iscl)-txym(iz,iscl))*dzu_i(izp1)
c
c ------------- gather all the budget terms
c
         Tv_prod(iz,iscl) = -2.0*wtle(iz,iscl)*dTdz
         Tv_diss(iz,iscl) = stat(iz,2)
         Tv_part1(iz,iscl) = -stat(iz,3)
         !Tv_part2 is already gathered
      enddo
c
      return
      end
      subroutine dealias
c
c --------- wave cutoff filter using 2d fft
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real wve(nny,jxs:jxe,izs:ize)
      real wves(nnxp2,iys:iye,izs:ize)
c
c --------- sharp spectral cutoff, specific to current 2dfft
c
      ix_cut   = 2*int(float(nnx)/3.) + 3
      iy_cut_l = int(float(nny)/3.) + 2
      iy_cut_u = nnyp2 - iy_cut_l
c
c ---------- u-equation
c
      call fft2d_mpi(u(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
      call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
      call fft2d_mpi(u(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
c
c ---------- v-equation
c
      call fft2d_mpi(v(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
      call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
      call fft2d_mpi(v(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
c
c ---------- w-equation
c
      call fft2d_mpi(w(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
      call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
      call fft2d_mpi(w(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
c
c ---------- e-equation
c
      call fft2d_mpi(e(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
      call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
      call fft2d_mpi(e(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
c
c ------------- scalars, not stored in correct order
c
      do iscl=1,nscl
         do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            wves(ix,iy,iz) = t(ix,iy,iscl,iz)
         enddo
         enddo
         enddo
         call fft2d_mpi(wves(1,iys,izs),wve(1,jxs,izs),trigx(1,1),
     +           trigc,nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
         call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
         call fft2d_mpi(wves(1,iys,izs),wve(1,jxs,izs),trigx(1,1),
     +           trigc,nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
         do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            t(ix,iy,iscl,iz) = wves(ix,iy,iz)
         enddo
         enddo
         enddo
      enddo
c
      return
      end
      subroutine pressure
c
c -------- solve for pressure using a matrix transpose
c          across mpi tasks and tridiagonal solver. 
c          The transposed array
c          is dimensioned (0:nnz+1). Values 
c          (0 & nnz+1) are not needed but are useful in the 
c          matrix transpose when we return (see send_ztox).
c          On exit p is defined at all [izs-1:ize+1].
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real pfft(nny,jxs:jxe,izs-1:ize+1)
      real pt(0:nnz+1,jxs:jxe,jys:jye)
      real ptopfft(nny,jxs:jxe,1:2)
      real psum(1:nnz)
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
c ------------ Fourier analyze the right hand side
c              at all iz = izs,ize. results are in pfft
c
c


      call fft2d_mpi(p(1,iys,izs),pfft(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)

c
c ------------ Fourier analyze the radiation bc arrays
c
      if(ibcu .eq. 1) then
        call fft2d_mpi(ptop(1,iys,1),ptopfft(1,jxs,1),
     +           trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           1,2,myid,ncpu_s,numprocs,-2)
      endif
c
c ---------- transpose first and last index of array
c            the order of pfft is (y,x,z)
c
      call xtoz_trans(pfft,pt,nny,nnz,jys,jye,jy_s,jy_e,
     +                jxs,jxe,izs,ize,iz_s,iz_e,myid,ncpu_s,
     +                numprocs)


      call solve_trid(pt, ptopfft)

c
c ------------- transpose back
c
      call ztox_trans(pt,pfft,nny,nnz,jys,jye,jy_s,jy_e,
     +                jxs,jxe,izs,ize,iz_s,iz_e,myid,ncpu_s,
     +                numprocs)
c
      iz_ee = ize+1
      if(ise .eq. numprocs-1) then
         iz_ee = ize
      endif
c
c --------- inverse fft at all iz=izs,iz_ee to get p
c           see z indices
c
      call fft2d_mpi(p(1,iys,izs),pfft(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,iz_ee,myid,ncpu_s,numprocs,2)
c
c -------- partial sums for pressure
c
      do iz=1,nnz
         psum(iz) = 0.0
      enddo
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            psum(iz) = psum(iz) + p(ix,iy,iz)
         enddo
         enddo
         psum(iz) = psum(iz)*fnxy
      enddo
      call mpi_sum_z(psum,i_root,myid,nnz,1)
c
      do iz=izs,iz_ee
c        psum(iz) = -psum(iz) + engz(iz) + c23*engsbz(iz)
         do iy=iys,iye
         do ix=1,nnx
            p(ix,iy,iz) = p(ix,iy,iz) - psum(iz)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine extra_flux_terms
      use pars
      use particles
      use fields
      use con_data
      use con_stats
      use fftwk
      implicit none

      real :: stat(1:nnz,3)
      real :: weit,weit1
      real :: Tflucp1,Tflucm1,Tfluc,Tmean
      real :: qflucp1,qflucm1,qfluc,qmean,dqpdz
      real :: Sq,Sqp1
      real :: gradTp(3),dTpdz1,dTpdz,dTdz
      integer :: iz,i,j,izp1,izm1,iscl

!     Compute the "extra" enthalpy budget terms, located at w locations:
c -------- stat(.,1) = < q' T' w' > 
c          stat(.,2) = < T' Sq > where Sq is q source
c          stat(.,3) = Dv*< T' dq'/dz > 

       

      stat = 0.0
      do iz=izs,ize
         izp1 = iz + 1
         izm1 = iz - 1
         weit = dzw(iz)/(dzw(iz) + dzw(izp1))
         weit1 = 1.0 - weit

       do j=iys,iye
       do i=1,nnx

         if (iz==1)  then
           Tmean = 2.0*Tbot(1) - txym(iz,1)
           Tflucm1 = t(i,j,1,izm1)-Tmean
           Tflucp1 = t(i,j,1,izp1)-txym(izp1,1)
         elseif (iz==nnz) then
           Tmean = 2.0*Ttop(1) - txym(iz,1)
           Tflucp1 = t(i,j,1,izp1)-Tmean
           Tflucm1 = t(i,j,1,izm1)-txym(izm1,1)
         else
           Tflucp1 = t(i,j,1,izp1)-txym(izp1,1)
           Tflucm1 = t(i,j,1,izm1)-txym(izm1,1)
         end if
         Tfluc = t(i,j,1,iz)-txym(iz,1)

         if (iz==1)  then
           qmean = 2.0*Tbot(2) - txym(iz,2)
           qflucm1 = t(i,j,2,izm1)-qmean
           qflucp1 = t(i,j,2,izp1)-txym(izp1,2)
           Sqp1 = -partHsrc(i,j,izp1)
         elseif (iz==nnz) then
           qmean = 2.0*Ttop(2) - txym(iz,2)
           qflucp1 = t(i,j,2,izp1)-qmean
           qflucm1 = t(i,j,2,izm1)-txym(izm1,2)
           Sqp1 = 0.0
         else
           qflucp1 = t(i,j,2,izp1)-txym(izp1,2)
           qflucm1 = t(i,j,2,izm1)-txym(izm1,2)
         end if
         qfluc = t(i,j,2,iz)-txym(iz,2)
         Sq = -partHsrc(i,j,iz)

         !Get dq'/dz at w-locations:
         dqpdz = (qflucp1 - qfluc)*dzu_i(izp1)

         !Then can get source #3:
         stat(iz,3) = stat(iz,3) + vis_s(i,j,2,iz)*
     +   (weit1*Tfluc + weit*Tflucp1)*dqpdz


         !Now source #1, which must be evaluated at w-location:

         stat(iz,1) = stat(iz,1) + w(i,j,iz)*
     +   (weit1*Tfluc + weit*Tflucp1)*
     +   (weit1*qfluc + weit*qflucp1)


         !Now source #2:
         stat(iz,2) = stat(iz,2) + 
     +   (weit1*Tfluc + weit*Tflucp1)*
     +   (weit1*Sq + weit*Sqp1)

     
       end do
       end do

         stat(iz,1) = stat(iz,1)*fnxy
         stat(iz,2) = stat(iz,2)*fnxy
         stat(iz,3) = stat(iz,3)*fnxy
      end do
         

      call mpi_sum_z(stat(1,1),i_root,myid,nnz*3,1)


c
c ------ we have all terms on all processors for all z, add them up
c
      do iz=1,nnz
c
c ------------- gather all the budget terms
c
         trip(iz) = stat(iz,1)
         TpSq(iz) = stat(iz,2)
         Tpdqp(iz) = stat(iz,3)
      enddo
c
      return
      end
      subroutine forcing
c
c ----------- update surface temperature based on a 
c             constant cooling rate
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
      tsfcc(1) = t_surf_i - c_rate*time
c
      return
      end
